// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cmath>
#include <vector>
#include <list>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
#include <boost/thread.hpp>
#endif
#include "data.h"
#include "bpmatrix.h"
#include "../common/cyktable.h"
#include "../common/glob.h"
#include "fa.h"
#include "maf.h"
#include "aln.h"

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
};
};

typedef boost::shared_ptr<BPMatrix> BPMatrixPtr;

static uint folding_method_;
static float th_;
static uint win_sz_;
static uint pair_sz_;
static bool use_pf_scale_;

void set_folding_method(uint method)
{
  folding_method_ = method;
}

void set_bp_threshold(float th)
{
  th_ = th;
}

void set_window_size(uint win_sz, uint pair_sz)
{
  win_sz_ = win_sz;
  pair_sz_ = pair_sz;
  if (pair_sz_==0 || win_sz_<pair_sz_)
    pair_sz_ = win_sz_;
}

void set_use_pf_scale(bool use_pf_scale)
{
  use_pf_scale_ = use_pf_scale;
}

template < class Seq, class BPM, class Node >
static
void
make_tree(std::vector<Node>& tree, const Seq& seq,
	  const BPM& pf, double th);

template < class Node >
static
void
find_root(std::vector<uint>& root, const std::vector<Node>& tree);

template <class BPF>
static
void
fill_weight(const BPF& pf, std::vector<float>& weight);

template < class Node >
static
void
find_max_parent(std::vector<uint>& max_pa, const std::vector<Node>& x);

template <class IS>
static
BPMatrixPtr
make_bp_matrix(const IS& s)
{
  return BPMatrixPtr();
}

template <class IS, class OS>
static
void
convert_seq(const IS& in, OS& out)
{
  char2rna(out, in);
}

template < class S, class IS, class N >
Data<S,IS,N>::
Data(const IS& s)
  : tree(), seq(), root(), weight(), max_pa()
{
  if (folding_method_!=NO_BPMATRIX) {
    BPMatrixPtr bp = make_bp_matrix(s);
    convert_seq(s, seq);
    make_tree(tree, s, *bp, th_);
    find_root(root, tree);
    weight.resize(seq.size());
    fill_weight(*bp, weight);
    find_max_parent(max_pa, tree);
  } else {
    convert_seq(s, seq);
  }
}
 
// helper functions 

template < class BPM, class Seq, class Node >
static
uint
make_tree_helper(std::vector<Node>& tree,
		 const Seq& seq,
		 CYKTable<uint>& vt,
		 const BPM& pf,
		 const CYKTable< std::list<Pos> >& table,
		 const Pos& pos)
{
  typedef typename Node::Edge Edge;
  if (vt(pos.first, pos.second)==static_cast<uint>(-1)) {
    if (pos.first==pos.second) { // leaf node
      Node leaf(pos, 1.0);
      tree.push_back(leaf);
      vt(pos)=tree.size()-1;
    } else {
      const std::list<Pos>& cur = table(pos);
      if (cur.empty()) {
	Node node(pos, seq[pos.first], seq[pos.second],
		  pf(pos.first+1, pos.second+1), 1);
	uint ret=make_tree_helper(tree, seq, vt, pf, table,
				  Pos(pos.first, pos.first));
	node[0] = Edge(ret, pos);
	tree.push_back(node);
	vt(pos)=tree.size()-1;
      } else {
	Node node(pos, seq[pos.first], seq[pos.second],
		  pf(pos.first+1, pos.second+1), cur.size());
	std::list<Pos>::const_iterator x;
	uint i;
	for (x=cur.begin(), i=0; x!=cur.end(); ++x, ++i) {
	  uint ret=make_tree_helper(tree, seq, vt, pf, table, *x);
	  node[i] = Edge(ret, pos, *x);
	}
	tree.push_back(node);
	vt(pos)=tree.size()-1;
      }
    }
  }
  return vt(pos.first, pos.second);
}

struct is_child : public std::binary_function<Pos,Pos,bool> {
  bool operator()(const Pos& pa, const Pos& ch) const
  {
    return /*pa.first < ch.first &&*/ pa.second > ch.second;
  }
};

template < class Seq, class BPM, class Node >
static
void
make_tree(std::vector<Node>& tree, const Seq& seq,
	  const BPM& pf, double th)
{
  // scan the matrix in bottom up order
  CYKTable< std::list<Pos> > bp(seq.size());
  CYKTable< std::list<Pos> > ch(seq.size());
  std::vector< std::list<Pos> > head(seq.size());
  for (uint j=1; j!=seq.size(); ++j)  {
    for (uint i=j-1; ; --i) {
      if (pf(i+1,j+1)>=th) {
#if 0
	std::back_insert_iterator< std::list<Pos> > ii(bp(i,j));
	std::copy(ch(i+1,j-1).begin(), ch(i+1,j-1).end(), ii);
#else
	std::swap(bp(i,j),ch(i+1,j-1));
#endif
	ch(i,j).push_back(Pos(i,j));
	head[i].push_back(Pos(i,j));
      } else {
	std::back_insert_iterator< std::list<Pos> > ii(ch(i,j));
	std::remove_copy_if(ch(i+1,j).begin(), ch(i+1,j).end(), ii,
			    std::bind1st(is_child(), head[i].back()));
	std::copy(head[i].begin(), head[i].end(), ii);
      }
      if (i==0) break;
    }
  }

  // trace back from root nodes
  CYKTable<uint> vt(seq.size());
  vt.fill(static_cast<uint>(-1));
  for (uint i=0; i!=seq.size(); ++i) {
    std::list<Pos>::const_reverse_iterator x;
    for (x=head[i].rbegin(); x!=head[i].rend(); ++x) {
      make_tree_helper(tree, seq, vt, pf, bp, *x);
    }
  }
}

template < class Node >
static
void
find_root(std::vector<uint>& root, const std::vector<Node>& tree)
{
  std::vector<bool> is_root(tree.size(), true);
  for (uint i=0; i!=tree.size(); ++i) {
    typename Node::const_iterator ix;
    for (ix=tree[i].begin(); ix!=tree[i].end(); ++ix) {
      is_root[ix->to()]=false;
      //std::cout << i << "->" << ix->to() << std::endl;
    }
  }
  uint n=0;
  for (uint i=0; i!=is_root.size(); ++i) {
    if (is_root[i]) n++;
  }
  root.resize(n);
  uint c=0;
  for (uint i=0; i!=is_root.size(); ++i) {
    if (is_root[i]) root[c++]=i;
  }
}

template <class BPF>
static
void
fill_weight(const BPF& pf, std::vector<float>& weight)
{
  for (uint i=0; i!=weight.size(); ++i) {
    weight[i] = 1.0;
    for (uint j=0; j!=i; ++j)
      weight[i] -= pf(j+1, i+1);
    for (uint j=i+1; j!=weight.size(); ++j)
      weight[i] -= pf(i+1, j+1);
    if (weight[i]<0.0) weight[i] = 0.0;
  }
}

template < class Node >
static
void
find_max_parent(std::vector<uint>& max_pa, const std::vector<Node>& x)
{
  max_pa.resize(x.size());
  std::fill(max_pa.begin(), max_pa.end(), static_cast<uint>(-1));
  typename Node::const_iterator ix;
  for (uint i=0; i!=x.size(); ++i) {
    for (ix=x[i].begin(); ix!=x[i].end(); ++ix) {
      if (max_pa[ix->to()]==static_cast<uint>(-1) || max_pa[ix->to()]<i) {
	max_pa[ix->to()] = i;
      }
    }
  }
}

static
BPMatrixPtr
make_bp_matrix_helper(const std::string &s, uint method)
{
  switch (method) {
  case FOLD:
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      if (use_pf_scale_) {
	int length = s.size();
	char* str = new char[length+1];
	double min_en = Vienna::fold(s.c_str(), str);
	double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
	Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
	delete[] str;
      } else {
	Vienna::pf_scale = -1;
      }
      Vienna::init_pf_fold(s.size());
      Vienna::pf_fold(const_cast<char*>(s.c_str()), NULL);
      BPMatrixPtr bp(new BPMatrix(s.size(), Vienna::pr, Vienna::iindx));
      Vienna::free_pf_arrays();
      return bp;
    }
    break;
  case LFOLD:
#if 0
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      Vienna::plist *pl;
      Vienna::pf_scale = -1;
      uint wsz = s.size()<win_sz_ ? s.size() : win_sz_;
      uint psz = wsz<pair_sz_ ? wsz : pair_sz_;
      Vienna::init_pf_foldLP(s.size());
      Vienna::pfl_fold(const_cast<char*>(s.c_str()), wsz, psz, th_, &pl);
      BPMatrixPtr bp(new BPMatrix(s.size()));
      for (uint k=0; pl[k].i!=0; ++k)
	(*bp)(pl[k].i, pl[k].j) = pl[k].p;
      free(pl);
      Vienna::free_pf_arraysLP();
      return bp;
    }
    break;
#endif
  default:
    assert(!"unsupported folding method");
    break;
  }
  return BPMatrixPtr();
}

template < >
//static
BPMatrixPtr
make_bp_matrix(const std::string& x)
{
  std::string s(x);
  boost::algorithm::to_lower(s);
  switch (folding_method_) {
  case FOLD:
  case LFOLD:
    return make_bp_matrix_helper(s, folding_method_);
    break;
  default:
    assert(!"unsupported folding method");
    break;
  }
  return BPMatrixPtr();
}

template < >
//static
BPMatrixPtr
make_bp_matrix(const MASequence<std::string>& x)
{
  switch (folding_method_) {
  case ALIFOLD:
    {
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
      static boost::mutex mtx;
      boost::mutex::scoped_lock lock(mtx);
#endif
      // prepare an alignment
      uint length = x.get_seq(0).size();
      char** seqs = new char*[x.n_seqs()+1];
      seqs[x.n_seqs()]=NULL;
      for (uint i=0; i!=x.n_seqs(); ++i) {
	std::string s(x.get_seq(i));
	assert(s.size()==length);
	seqs[i] = new char[length+1];
	strcpy(seqs[i], s.c_str());
      }
      {
	// scaling parameters to avoid overflow
	char* str = new char[length+1];
	double min_en = Vienna::alifold(seqs, str);
	delete[] str;
	Vienna::free_alifold_arrays();
	double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
	Vienna::pf_scale = exp(-(1.07*min_en)/kT/length);
      }
      // build a base pair probablity matrix
      Vienna::pair_info* pi;
      Vienna::alipf_fold(seqs, NULL, &pi);
      BPMatrixPtr ret(new BPMatrix(length));
      for (uint k=0; pi[k].i!=0; ++k)
	(*ret)(pi[k].i, pi[k].j) = pi[k].p;
      free(pi);
      for (uint i=0; i!=x.n_seqs(); ++i) delete[] seqs[i];
      delete[] seqs;
      return ret;
    }
    break;

  case FOLD:
  case LFOLD:
    {
      std::list<std::string> ali;
      std::list<BPMatrixPtr> bps;
      for (uint i=0; i!=x.n_seqs(); ++i) {
	std::string s(x.get_seq(i));
	boost::algorithm::to_lower(s);
	ali.push_back(s);
	s=erase_gap(s);
	bps.push_back(make_bp_matrix_helper(s, folding_method_));
      }
      BPMatrixPtr ret(new BPMatrix(ali, bps));
      return ret;
    }
    break;
  default:
    assert(!"unsupported folding method");
    break;
  }

  return BPMatrixPtr();
}

template <>
//static
void
convert_seq(const MASequence<std::string>& in, std::vector<Col>& out)
{
  for (uint i=0; i!=in.n_seqs(); ++i) {
    RNASequence a;
    char2rna(a, in.get_seq(i));
    out.resize(a.size());
    for (uint j=0; j!=out.size(); ++j) {
      if (a[j] != RNASymbol<RNASequence::value_type>::GAP) {
	out[j].cnt(index(a[j]))++;
      }
    }
  }
}

template <>
//static
void
convert_seq(const MASequence<std::string>& in, MASequence<RNASequence>& out)
{
  for (uint i=0; i!=in.n_seqs(); ++i) {
    RNASequence a;
    char2rna(a, in.get_seq(i));
    out.add_seq(a);
  }
}

enum { TP_UNKNOWN, TP_FA, TP_ALN, TP_MAF };

static
uint
check_filetype(const char* f)
{
#if 0  
  boost::spirit::file_iterator<> fi(f);
  std::string s;
  MASequence<std::string> ma;
  if (load_fa(s, fi)) return TP_FA;
  if (load_aln(ma, fi)) return TP_ALN;
  if (load_maf(ma, fi)) return TP_MAF;
  return TP_UNKNOWN;
#else
  std::ifstream in(f);
  std::string l;
  while (std::getline(in, l)) {
    if (l[0]=='>') return TP_FA;
    else if (l.compare(0, 7, "CLUSTAL")==0) return TP_ALN;
    else if (l.compare(0, 2, "a ")==0) return TP_MAF;
  }
  return TP_UNKNOWN;
#endif
}

MakeData<SData>::
MakeData(const char* f)
  : filename(f), type(check_filetype(f)), fi(f)
{ }

bool
MakeData<SData>::
operator()(SData& data, bool skip)
{
  if (!fi) {
    std::ostringstream os;
    os << filename << ": no such file";
    throw os.str().c_str();
  }
  if (type==TP_FA) {
    std::string s;
    if (load_fa(s, fi)) {
      if (!skip) data=SData(s);
      return true;
    } else {
      return false;
    }
#if 0
  } else {
    std::ostringstream os;
    os << filename << ": bad format";
    throw os.str().c_str();
#endif
  }
  return false;
}

MakeData<MData>::
MakeData(const char* f)
  : filename(f), type(check_filetype(f))
{
  switch (type) {
  case TP_ALN:
    fi = boost::spirit::file_iterator<>(f);
    available = fi;
    break;
  case TP_MAF:
    available = maf.open(f);
    break;
  default:
    break;
  }
}
  
bool
MakeData<MData>::
operator()(MData& data, bool skip)
{
  if (!available) {
    std::ostringstream os;
    os << filename << ": no such file";
    throw os.str().c_str();
    //return false;
  }
  bool ret=false;
  MASequence<std::string> ma;
  switch (type) {
  case TP_ALN:
    ret=load_aln(ma, fi);
    break;
  case TP_MAF:
    ret=maf.get(ma);
    break;
  default:
#if 0
    {
      std::ostringstream os;
      os << filename << ": bad format";
      throw os.str().c_str();
    }
#endif
    break;
  }
  if (ret && !skip) data=MData(ma);
  return ret;
}


template <>
bool
load_examples(const std::string& label,
	      const char* filename,
	      std::vector< std::pair<std::string, SData> >& ex,
	      const std::vector< uint >& sv_index,
	      uint n_th /*=1*/, uint th_no /*=0*/)
{
  std::list<std::string> fa;
  Glob glob(filename);
  if (glob.empty()) {
    std::ostringstream os;
    os << filename << ": no matches found";
    throw os.str().c_str();
    //return false;
  }
  Glob::const_iterator p;
  for (p=glob.begin(); p!=glob.end(); ++p) {
    if (!load_fa(fa, p->c_str())) {
#if 0
      std::ostringstream os;
      os << p->c_str() << ": bad format";
      throw os.str().c_str();
#endif
      return false;
    }
  }

  uint cnt=0;
  for (cnt=0; cnt!=sv_index.size(); ++cnt) {
    if (sv_index[cnt]>=ex.size()) break;
  }

  std::list<std::string>::const_iterator x;
  for (x=fa.begin(); x!=fa.end(); ++x) {
    if (sv_index[cnt]==ex.size()) {
      if (cnt++%n_th==th_no)
	ex.push_back(std::make_pair(label,SData(*x)));
      else
	ex.push_back(std::make_pair(label,SData()));
    } else {
      ex.push_back(std::make_pair(label,SData()));
    }
  }
  return true;
}

template <>
bool
load_examples(const std::string& label,
	      const char* filename,
	      std::vector< std::pair<std::string, SData> >& ex,
	      uint n_th /*=1*/, uint th_no /*=0*/)
{
  std::list<std::string> fa;
  Glob glob(filename);
  if (glob.empty()) {
    std::ostringstream os;
    os << filename << ": no matches found";
    throw os.str().c_str();
    //return false;
  } 
  Glob::const_iterator p;
  for (p=glob.begin(); p!=glob.end(); ++p) {
    if (!load_fa(fa, p->c_str())) {
#if 0
      std::ostringstream os;
      os << p->c_str() << ": bad format";
      throw os.str().c_str();
#endif
      return false;
    }
  }

  std::list<std::string>::const_iterator x;
  for (x=fa.begin(); x!=fa.end(); ++x) {
    if (ex.size()%n_th==th_no) {
      SData e(*x);
      ex.push_back(std::make_pair(label,e));
    } else {
      ex.push_back(std::make_pair(label,SData()));
    }
  }
  return true;
}

template <>
bool
load_examples(const std::string& label,
	      const char* filename,
	      std::vector< std::pair<std::string, MData> >& ex,
	      uint n_th /*=1*/, uint th_no /*=0*/)
{
  std::list< MASequence<std::string> > ma;
  Glob glob(filename);
  if (glob.empty()) {
    std::ostringstream os;
    os << filename << ": no matches found";
    throw os.str().c_str();
    //return false;
  }
  Glob::const_iterator p;
  for (p=glob.begin(); p!=glob.end(); ++p) {
    if (!load_maf(ma, p->c_str()) &&
	!load_aln(ma, p->c_str()) &&
	!load_fa(ma, p->c_str())) {
#if 0
      std::ostringstream os;
      os << p->c_str() << ": bad format";
      throw os.str().c_str();
#endif
      return false;
    }
  }

  std::list< MASequence<std::string> >::const_iterator x;
  for (x=ma.begin(); x!=ma.end(); ++x) {
    if (ex.size()%n_th==th_no)
      ex.push_back(std::make_pair(label,MData(*x)));
    else
      ex.push_back(std::make_pair(label,MData()));
  }
  return true;
}

template <>
bool
load_examples(const std::string& label,
	      const char* filename,
	      std::vector< std::pair<std::string, MData> >& ex,
	      const std::vector< uint >& sv_index,
	      uint n_th /*=1*/, uint th_no /*=0*/)
{
  std::list< MASequence<std::string> > ma;
  Glob glob(filename);
  if (glob.empty()) {
    std::ostringstream os;
    os << filename << ": no matches found";
    throw os.str().c_str();
    //return false;
  }
  Glob::const_iterator p;
  for (p=glob.begin(); p!=glob.end(); ++p) {
    if (!load_maf(ma, p->c_str()) &&
	!load_aln(ma, p->c_str()) &&
	!load_fa(ma, p->c_str())) {
#if 0
      std::ostringstream os;
      os << p->c_str() << ": bad format";
      throw os.str().c_str();
#endif
      return false;
    }
  }

  uint cnt=0;
  for (cnt=0; cnt!=sv_index.size(); ++cnt) {
    if (sv_index[cnt]>=ex.size()) break;
  }

  std::list< MASequence<std::string> >::const_iterator x;
  for (x=ma.begin(); x!=ma.end(); ++x) {
    if (sv_index[cnt]==ex.size()) {
      if (cnt++%n_th==th_no)
	ex.push_back(std::make_pair(label,MData(*x)));
      else
	ex.push_back(std::make_pair(label,MData()));
    } else {
      ex.push_back(std::make_pair(label,MData()));
    }
  }
  return true;
}

// instantiation

#include <string>
#include "bpmatrix.h"
#include "../common/rna.h"

template
class Data<RNASequence,std::string>;

template
class Data<std::vector<Col>, MASequence<std::string> >;
