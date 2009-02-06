// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "bpmatrix.h"
#include "../common/rna.h"
#include <cmath>
#include <stack>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
#include <boost/thread.hpp>
#endif

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/PS_dot.h>
#if 0
  extern int pfl_fold(char *sequence, int winSize, int pairdist,
		      float cutoff, struct plist **pl);
  extern void init_pf_foldLP(int length);
  extern void free_pf_arraysLP(void);
#endif
  extern char* pbacktrack(char *sequence);
  extern int   st_back;
};
};

#ifdef HAVE_LIBCONTRAFOLD
#include <contrafold.h>
#endif

typedef boost::shared_ptr<BPMatrix> BPMatrixPtr;
namespace po = boost::program_options;

// options
void
BPMatrix::Options::
add_options(po::options_description& desc)
{
  desc.add_options()
    ("noGU",
     po::value<bool>(&no_GU)->zero_tokens()->default_value(false),
     "disallow GU wobble base-pairs")
    ("noClosingGU",
     po::value<bool>(&no_closingGU)->zero_tokens()->default_value(false),
     "disallow closing GU base-pairs")
    ("noLonelyPairs",
     po::value<bool>(&no_LonelyPairs)->zero_tokens()->default_value(false),
     "disallow lonely base-pairs")
    ("use-alifold",
     po::value<bool>(&alifold)->zero_tokens()->default_value(false),
     "use pf_alifold for producing base pairing probability matrices")
#ifdef HAVE_LIBCONTRAFOLD
    ("use-contrafold",
     po::value<bool>(&contrafold)->zero_tokens()->default_value(false),
      "use the CONTRAfold model for producing base pairing probability matrices")
#endif
    ("pf-scale",
     po::value<bool>(&use_pf_scale_mfe)->zero_tokens()->default_value(false),
     "calculate appropriciate pf_scales using MFE")
#if 0
    ("sampling",
     po::value<uint>(&n_samples)->default_value(0),
     "use stochastic sampling for producing base pairing probability matrices")
#endif
#if 0
    ("window-size,w", po::value<uint>(&win_sz)->default_value(0),
     "set the window size for folding RNAs")
    ("pair-width", po::value<uint>(&pair_sz)->default_value(0),
     "set the pair width for pairing bases")
#endif
    ;
}

uint
BPMatrix::Options::
method() const
{
  uint m=FOLD;
  if (alifold) m=ALIFOLD;
  if (n_samples>0) m=SFOLD;
  if (contrafold) m=CONTRAFOLD;
  return m;
}

// helpers
static
bool
make_bp_matrix(BPMatrix& bp, const std::string &x, float pf_scale, 
	       const BPMatrix::Options& opts);

static
bool
make_bp_matrix(BPMatrix& bp, const std::list<std::string>& x, float pf_scale, 
	       const BPMatrix::Options& opts);

BPMatrix::
BPMatrix(const std::string& s, const Options& opts)
  : sz_(s.size()), table_(sz_+1)
{
  table_.fill(0.0);
  make_bp_matrix(*this, s, -1, opts);
}

BPMatrix::
BPMatrix(const std::string& s, float pf_scale, const Options& opts)
  : sz_(s.size()), table_(sz_+1)
{
  table_.fill(0.0);
  make_bp_matrix(*this, s, pf_scale, opts);
}

BPMatrix::
BPMatrix(const std::list<std::string>& ma, float pf_scale, const Options& opts)
  : sz_(ma.begin()->size()), table_(sz_+1)
{
  table_.fill(0.0);
  make_bp_matrix(*this, ma, pf_scale, opts);
}

BPMatrix::
BPMatrix(const std::list<std::string>& ma, const Options& opts)
  : sz_(ma.begin()->size()), table_(sz_+1)
{
  table_.fill(0.0);
  make_bp_matrix(*this, ma, -1, opts);
}


// folding for sigle sequences
static
bool
make_bp_matrix(BPMatrix& bp, const std::string &x, float pf_scale,
	       const BPMatrix::Options& opts)
{
  std::string s(x);
  boost::algorithm::to_lower(s);
  Vienna::noGU = opts.no_GU ? 1 : 0;
  Vienna::no_closingGU = opts.no_closingGU ? 1 : 0;
  Vienna::noLonelyPairs = opts.no_LonelyPairs ? 1 : 0;
  switch (opts.method()) {
    case BPMatrix::FOLD:
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      Vienna::pf_scale = pf_scale;
      if (opts.use_pf_scale_mfe) {
        int length = s.size();
        char* str = new char[length+1];
        double min_en = Vienna::fold(s.c_str(), str);
        double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
        Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
        delete[] str;
      }
      Vienna::init_pf_fold(s.size());
      Vienna::pf_fold(const_cast<char*>(s.c_str()), NULL);
      for (uint j=2; j!=bp.size()+1; ++j) {
        for (uint i=j-1; ; --i) {
          bp(i,j) = Vienna::pr[Vienna::iindx[i]-j];
          if (i==1) break;
        }
      }
      Vienna::free_pf_arrays();
      return true;
    }
    break;

    case BPMatrix::SFOLD:
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      int bk_st_back=Vienna::st_back;
      Vienna::st_back=1;

      Vienna::pf_scale = pf_scale;
      if (opts.use_pf_scale_mfe) {
        int length = s.size();
        char* str = new char[length+1];
        double min_en = Vienna::fold(s.c_str(), str);
        double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
        Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
        delete[] str;
      }
      Vienna::init_pf_fold(s.size());
      Vienna::pf_fold(const_cast<char*>(s.c_str()), NULL);

      // stochastic sampling
      for (uint n=0; n!=opts.n_samples; ++n) {
        char *str = Vienna::pbacktrack(const_cast<char*>(s.c_str()));
        //std::cout << s << std::endl << str << std::endl;
        assert(s.size()==strlen(str));
        // count basepairs
        std::stack<uint> st;
        for (uint i=0; i!=s.size(); ++i) {
          if (str[i]=='(') {
            st.push(i);
          } else if (str[i]==')') {
            bp(st.top()+1, i+1) += 1;
            st.pop();
          }
        }
        //std::cout << str << std:endl;
	
        free(str);
      }

      // normalize
      for (uint j=1; j!=bp.size(); ++j) {
        for (uint i=j-1; /*i>=0*/; --i) {
          bp(i+1,j+1) = bp(i+1,j+1) / opts.n_samples;
          if (i==0) break;
        }
      }

      Vienna::st_back=bk_st_back;
      Vienna::free_pf_arrays();
      return true;
    }
    break;
    
    case BPMatrix::LFOLD:
#if 0
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      Vienna::plist *pl;
      Vienna::pf_scale = pf_scale;
      if (opts.use_pf_scale_mfe) {
        int length = s.size();
        char* str = new char[length+1];
        double min_en = Vienna::fold(s.c_str(), str);
        double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
        Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
        delete[] str;
      }
      uint wsz = s.size()<win_sz_ ? s.size() : win_sz_;
      uint psz = wsz<pair_sz_ ? wsz : pair_sz_;
      Vienna::init_pf_foldLP(s.size());
      Vienna::pfl_fold(const_cast<char*>(s.c_str()), wsz, psz, th_, &pl);
      for (uint k=0; pl[k].i!=0; ++k)
        (*this)(pl[k].i, pl[k].j) = pl[k].p;
      free(pl);
      Vienna::free_pf_arraysLP();
      return true;
    }
    break;
#endif

#ifdef HAVE_LIBCONTRAFOLD
    case BPMatrix::CONTRAFOLD:
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      CONTRAfold<float> cf;
      std::vector<float> posterior;
      cf.ComputePosterior(s, posterior);
      uint k=0;
      for (uint i=0; i!=s.size()+1; ++i) {
        for (uint j=i; j!=s.size()+1; ++j) {
          if (i!=0) bp(i-1, j-1) = posterior[k];
          ++k;
        }
      }
    }
    break;
#endif

    default:
      assert(!"unsupported folding method");
      break;
  }
  return false;
}

template <class Seq>
static
void
make_index_map(const Seq& seq, std::vector<uint>& idxmap)
{
  typedef typename Seq::value_type rna_type;
  const rna_type GAP = RNASymbol<rna_type>::GAP;

  idxmap.resize(seq.size(), static_cast<uint>(-1));
  for (uint i=0, j=0; i!=seq.size(); ++i) {
    if (seq[i]!=GAP) idxmap[i]=j++;
  }
}

template <class Seq>
static
void
average_matrix(BPMatrix& bp,
	       const std::list<Seq>& ali,
	       const std::list< boost::shared_ptr<BPMatrix> >& bps)
{
  assert(ali.size()==bps.size());

  // align bp matrices according to the given alignment
  uint n_seq=ali.size();
  typename std::list<Seq>::const_iterator a;
  std::list< boost::shared_ptr<BPMatrix> >::const_iterator b;
  for (a=ali.begin(), b=bps.begin(); a!=ali.end() && b!=bps.end(); ++a, ++b) {
    std::vector<uint> idxmap;
    make_index_map(*a, idxmap);
    for (uint j=1; j!=bp.size(); ++j) {
      if (idxmap[j]!=static_cast<uint>(-1)) {
	for (uint i=j-1; /*i>=0*/; --i) {
	  if (idxmap[i]!=static_cast<uint>(-1)) {
	    bp(i+1,j+1) += (**b)(idxmap[i]+1,idxmap[j]+1);
	  }
	  if (i==0) break;
	}
      }
    }
    bp.add_matrix(*b, idxmap);
  }

  // averaging
  for (uint j=1; j!=bp.size(); ++j) {
    for (uint i=j-1; /*i>=0*/; --i) {
      bp(i+1,j+1) = bp(i+1,j+1) / n_seq;
      if (i==0) break;
    }
  }
}

// folding for aligned sequences
static
bool
make_bp_matrix(BPMatrix& bp, const std::list<std::string>& ma, float pf_scale,
	       const BPMatrix::Options& opts)
{
  Vienna::noGU = opts.no_GU ? 1 :0;
  Vienna::no_closingGU = opts.no_closingGU ? 1 : 0;
  Vienna::noLonelyPairs = opts.no_LonelyPairs ? 1 : 0;

  switch (opts.method()) {
    case BPMatrix::ALIFOLD:
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      // prepare an alignment
      uint length = ma.begin()->size();
      char** seqs = new char*[ma.size()+1];
      seqs[ma.size()]=NULL;
      std::list<std::string>::const_iterator x;
      uint i=0;
      for (x=ma.begin(); x!=ma.end(); ++x) {
	assert(x->size()==length);
	seqs[i] = new char[length+1];
	strcpy(seqs[i++], x->c_str());
      }
      if (pf_scale<0.0 || opts.use_pf_scale_mfe) {
	// scaling parameters to avoid overflow
	char* str = new char[length+1];
	double min_en = Vienna::alifold(seqs, str);
	delete[] str;
	Vienna::free_alifold_arrays();
	double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
	Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
      } else {
	Vienna::pf_scale = pf_scale;
      }
      // build a base pair probablity matrix
#ifdef HAVE_VIENNA18
      Vienna::plist* pi;
#else
      Vienna::pair_info* pi;
#endif
      Vienna::alipf_fold(seqs, NULL, &pi);
      for (uint k=0; pi[k].i!=0; ++k)
	bp(pi[k].i, pi[k].j) = pi[k].p;
      free(pi);
      for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
      delete[] seqs;
      return true;
    }
    break;

    case BPMatrix::FOLD:
    case BPMatrix::LFOLD:
    case BPMatrix::SFOLD:
    case BPMatrix::CONTRAFOLD:
    {
      std::list<std::string> ali;
      std::list<BPMatrixPtr> bps;
      std::list<std::string>::const_iterator x;
      for (x=ma.begin(); x!=ma.end(); ++x) {
	std::string s(*x);
	boost::algorithm::to_lower(s);
	ali.push_back(s);
	s=erase_gap(s);
	BPMatrixPtr b(new BPMatrix(s, pf_scale, opts));
	bps.push_back(b);
      }
      average_matrix(bp, ali, bps);
      return true;
    }
    break;

    default:
      assert(!"unsupported folding method");
      break;
  }

  return false;
}
