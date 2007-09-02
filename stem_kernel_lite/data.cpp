// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cmath>
#include <vector>
#include <list>
#include <sstream>
#include "data.h"
#include "bpmatrix.h"
#include "../common/profile.h"
#include "../common/cyktable.h"
#include "../common/fa.h"
#include "../common/maf.h"
#include "../common/aln.h"

template < class Seq >
static
uint
seq_size(const Seq& seq)
{
  return seq.size();
}

template < >
static
uint
seq_size(const std::list<std::string>& seq)
{
  return seq.begin()->size();
}

class Profiler
{
public:
  Profiler(const std::string& seq, const BPMatrix& bpm, float w=1.0)
    : seq_(seq), bpm_(bpm), w_(w), pr_(seq),
      idx_(seq_.size(), static_cast<uint>(-1)),
      nbp_(seq_.size(), 1.0)
  {
    make_idxmap();
    non_bp_profile();
  }

  Profiler(const Profiler& p)
    : seq_(p.seq_), bpm_(p.bpm_), w_(p.w_), pr_(p.pr_),
      idx_(p.idx_), nbp_(p.nbp_)
  {
  }

  float loop_profile(uint i) const
  {
    assert(idx_[i]!=static_cast<uint>(-1));
    return w_*nbp_[i];
  }

  uint index(uint i) const { return idx_[i]; }
  float weight() const { return w_; }

  void bp_profile(uint i, uint j, std::map<DAG::bp_t, float>& v) const
  {
    if (idx_[i]!=static_cast<uint>(-1) &&
	idx_[j]!=static_cast<uint>(-1)) {
      float p = bpm_.size()!=seq_.size() ?
	bpm_(idx_[i]+1,idx_[j]+1) : bpm_(i+1,j+1);

      for (uint a=0; a!=N_RNA; ++a) {
	if (pr_[i][a]==0.0) continue;
	for (uint b=0; b!=N_RNA; ++b) {
	  if (pr_[j][b]==0.0) continue;
	  DAG::bp_t k = std::make_pair(a,b);
	  std::map<DAG::bp_t, float>::iterator x = v.find(k);
	  if (x==v.end()) {
	    v.insert(std::make_pair(k, w_*p*pr_[i][a]*pr_[j][b]));
	  } else {
	    x->second += w_*p*pr_[i][a]*pr_[j][b];
	  }
	}
      }
    }
  }

private:
  void make_idxmap()
  {
    typedef std::string::value_type rna_type;
    const rna_type GAP = RNASymbol<rna_type>::GAP;
    uint j=0;
    for (uint i=0; i!=seq_.size(); ++i) {
      if (seq_[i]!=GAP) idx_[i]=j++;
    }
  }

  void non_bp_profile()
  {
    if (bpm_.size()!=seq_.size()) {
      for (uint i=0; i!=nbp_.size(); ++i) {
	//nbp_[i] = 1.0;
	if (idx_[i]!=static_cast<uint>(-1)) {
	  for (uint j=0; j!=i; ++j)
	    if (idx_[j]!=static_cast<uint>(-1))
	      nbp_[i] -= bpm_(idx_[j]+1, idx_[i]+1);
	  for (uint j=i+1; j!=nbp_.size(); ++j) 
	    if (idx_[j]!=static_cast<uint>(-1))
	      nbp_[i] -= bpm_(idx_[i]+1, idx_[j]+1);
	  if (nbp_[i]<0.0) nbp_[i] = 0.0;
	}
      }
    } else {			// the case of using pf_alifold
      for (uint i=0; i!=nbp_.size(); ++i) {
	//nbp_[i] = 1.0;
	if (idx_[i]!=static_cast<uint>(-1)) {
	  for (uint j=0; j!=i; ++j)
	    if (idx_[j]!=static_cast<uint>(-1))
	      nbp_[i] -= bpm_(j+1, i+1);
	  for (uint j=i+1; j!=nbp_.size(); ++j)
	    if (idx_[j]!=static_cast<uint>(-1))
	      nbp_[i] -= bpm_(i+1, j+1);
	  if (nbp_[i]<0.0) nbp_[i] = 0.0;
	}
      }
    }
  }

private:
  const std::string& seq_;
  const BPMatrix& bpm_;
  float w_;
  ProfileSequence pr_;
  std::vector<uint> idx_;
  std::vector<float> nbp_;
};

struct is_child : public std::binary_function<Pos,Pos,bool> {
  bool operator()(const Pos& pa, const Pos& ch) const
  {
    return /*pa.first < ch.first &&*/ pa.second > ch.second;
  }
};

template < class Node >
class DAGBuilder
{
public:
  DAGBuilder(const std::list<Profiler>& prof, const BPMatrix& bpm, float th)
    : prof_(prof), bpm_(bpm), th_(th), sz_(/*seq_size(seq_)*/ bpm.size()),
      bp_(sz_), head_(sz_), vt_(sz_)
  {
  }
  
  void build(std::vector<Node>& tree)
  {
    initialize();
    for (uint i=0; i!=sz_; ++i) {
      std::list<Pos>::const_reverse_iterator x;
      for (x=head_[i].rbegin(); x!=head_[i].rend(); ++x) {
	build_helper(tree, *x);
      }
    }
  }

private:
  typedef typename Node::Edge Edge;

  void initialize()
  {
    // scan the matrix in bottom up order
    CYKTable< std::list<Pos> > ch(sz_);
    for (uint j=1; j!=sz_; ++j)  {
      for (uint i=j-1; ; --i) {
	//std::cout <<  pf(i+1,j+1) << std::endl;
	if (bpm_(i+1,j+1)>=th_) {
#if 0
	  std::back_insert_iterator< std::list<Pos> > ii(bp_(i,j));
	  std::copy(ch(i+1,j-1).begin(), ch(i+1,j-1).end(), ii);
#else
	  std::swap(bp_(i,j),ch(i+1,j-1));
#endif
	  ch(i,j).push_back(Pos(i,j));
	  head_[i].push_back(Pos(i,j));
	} else {
	  std::back_insert_iterator< std::list<Pos> > ii(ch(i,j));
	  std::remove_copy_if(ch(i+1,j).begin(), ch(i+1,j).end(), ii,
			      std::bind1st(is_child(), head_[i].back()));
	  std::copy(head_[i].begin(), head_[i].end(), ii);
	}
	if (i==0) break;
      }
    }
    vt_.fill(static_cast<uint>(-1));
  }

  void make_leaf(std::vector<Node>& tree, const Pos& pos)
  {
    Node leaf(pos);
    tree.push_back(leaf);
    vt_(pos)=tree.size()-1;
  }

  void make_loop(std::vector<Node>& tree, const Pos& pos)
  {
    std::list<DAG::bp_freq_t> bp_freq;
    bp_profile(pos.first, pos.second, bp_freq);
    float n_w =  loop_profile(pos.first) * loop_profile(pos.second);
    Node node(pos, n_w, bp_freq, 1);
    uint ret = build_helper(tree, Pos(pos.first, pos.first));
    float e_w = edge_score(pos);
    node[0] = Edge(ret, pos, e_w);
    tree.push_back(node);
    vt_(pos)=tree.size()-1;
  }

  void make_stem(std::vector<Node>& tree, const Pos& pos)
  {
    const std::list<Pos>& cur = bp_(pos);
    std::list<DAG::bp_freq_t> bp_freq;
    bp_profile(pos.first, pos.second, bp_freq);
    float n_w = loop_profile(pos.first) * loop_profile(pos.second);
    Node node(pos, n_w, bp_freq, cur.size());
    std::list<Pos>::const_iterator x;
    uint i;
    for (x=cur.begin(), i=0; x!=cur.end(); ++x, ++i) {
      uint ret = build_helper(tree, *x);
      float e_w = edge_score(pos, *x);
      node[i] = Edge(ret, pos, *x, e_w);
    }
    tree.push_back(node);
    vt_(pos)=tree.size()-1;
  }

  uint build_helper(std::vector<Node>& tree, const Pos& pos)
  {
    // parse nodes in the depth first order
    if (vt_(pos.first, pos.second)==static_cast<uint>(-1)) {
      if (pos.first==pos.second) {
	make_leaf(tree, pos);
      } else if (bp_(pos).empty()) {
	make_loop(tree, pos);
      } else {
	make_stem(tree, pos);
      }
    }
    return vt_(pos.first, pos.second);
  }

  void bp_profile(uint i, uint j, std::list<DAG::bp_freq_t>& bp_freq) const
  {
    std::map<DAG::bp_t, float> v;
    std::list<Profiler>::const_iterator x;
    float t = 0.0;
    for (x=prof_.begin(); x!=prof_.end(); ++x) {
      if (x->index(i)!=static_cast<uint>(-1) &&
	  x->index(j)!=static_cast<uint>(-1)) {
	x->bp_profile(i, j, v);
      }
      t += x->weight();
    }
    std::map<DAG::bp_t, float>::const_iterator y;
    for (y=v.begin(); y!=v.end(); ++y) {
      bp_freq.push_back(std::make_pair(y->first, y->second/t));
    }
  }

  float loop_profile(uint i) const
  {
    float v = 0.0;
    float t = 0.0;
    std::list<Profiler>::const_iterator x;
    for (x=prof_.begin(); x!=prof_.end(); ++x) {
      if (x->index(i)!=static_cast<uint>(-1)) {
	v += x->loop_profile(i);
      }
      t += x->weight();
    }
    return v/t;
  }

  float edge_score(const Pos& p_pos)
  {
    assert(p_pos.first<p_pos.second);
    float ret = 1.0;
    for (uint i=p_pos.first+1; i!=p_pos.second; ++i)
      ret *= loop_profile(i);
    return ret;
  }

  float edge_score(const Pos& p_pos, const Pos& c_pos) const
  {
    assert(p_pos.first<c_pos.first);
    assert(c_pos.second<p_pos.second);
    float ret = 1.0;
    for (uint i=p_pos.first+1; i!=c_pos.first; ++i)
      ret *= loop_profile(i);
    for (uint i=c_pos.second+1; i!=p_pos.second; ++i)
      ret *= loop_profile(i);
    return ret;
  }

private:
  const std::list<Profiler>& prof_;
  const BPMatrix& bpm_;
  const float th_;
  uint sz_;
  CYKTable< std::list<Pos> > bp_;
  std::vector< std::list<Pos> > head_;
  CYKTable<uint> vt_;
};


template < class Node >
static
void
find_root(std::vector<uint>& root, const std::vector<Node>& tree);

template < class Node >
static
void
find_max_parent(std::vector<uint>& max_pa, const std::vector<Node>& x);


template < class S, class IS, class N >
Data<S,IS,N>::
Data(const IS& s, const BPMatrix::Options& opts)
  : tree(), seq(s), root(), max_pa()
{
  BPMatrix bp(s, opts);
  std::list<Profiler> prof;
  typename IS::const_iterator x;
  if (bp.n_matrices()>1) {
    BPMatrix::matrix_iterator bp_it = bp.matrix_begin();
    for (x=s.begin(); x!=s.end(); ++x, ++bp_it)
      prof.push_back(Profiler(*x, **bp_it, 1.0));
  } else {
    for (x=s.begin(); x!=s.end(); ++x)
      prof.push_back(Profiler(*x, bp, 1.0));
  }
  DAGBuilder<N> builder(prof, bp, opts.th);
  builder.build(tree);
  find_root(root, tree);
  find_max_parent(max_pa, tree);
}
 
template < class S, class IS, class N >
Data<S,IS,N>::
Data(const IS& s)
  : tree(), seq(s), root(), max_pa()
{
}

// helper functions
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

#if 0
DataLoader<SData>::
DataLoader(const char* filename,
	   const BPMatrix::Options& bp_opts, bool use_bp)
  : bp_opts_(bp_opts),
    use_bp_(use_bp),
    filename_(filename),
    type_(check_filetype(filename)),
    fi_(filename)
{
  if (!fi_) {
    std::ostringstream os;
    os << filename_ << ": no such file";
    throw os.str().c_str();
  }
}

SData*
DataLoader<SData>::
get()
{
  std::string s;
  switch (type_) {
  case TP_FA:
    if (load_fa(s, fi_)) {
      if (use_bp_)
	return new SData(s, bp_opts_);
      else
	return new SData(s);
    } else {
      return NULL;
    }
    break;
  default:
    return NULL;
    break;
  }
  return NULL;
}
#endif

DataLoader<MData>::
DataLoader(const char* filename,
	   const BPMatrix::Options& bp_opts, bool use_bp)
  : bp_opts_(bp_opts),
    use_bp_(use_bp),
    filename_(filename),
    type_(check_filetype(filename)),
    fi_()
{
  switch (type_) {
  case TP_FA:
  case TP_ALN:
  case TP_MAF:
    fi_ = boost::spirit::file_iterator<>(filename_);
    break;
  default:
    break;
  }
  if (!fi_) {
    std::ostringstream os;
    os << filename_ << ": no such file";
    throw os.str().c_str();
    //return false;
  }
}

MData*
DataLoader<MData>::
get()
{
  bool ret=false;
  std::list<std::string> ma;
  switch (type_) {
  case TP_FA:
    ret=load_fa(ma, fi_);
    break;
  case TP_ALN:
    ret=load_aln(ma, fi_);
    break;
  case TP_MAF:
    ret=load_maf(ma, fi_);
    break;
  default:
    break;
  }
  if (ret) {
    if (use_bp_)
      return new MData(ma, bp_opts_);
    else
      return new MData(ma);
  }
  return NULL;
}

#if 0
DataLoader<SData>*
DataLoaderFactory< DataLoader<SData> >::
get_loader(const char* filename) const
{
  switch (check_filetype(filename)) {
  case TP_FA:
    return new DataLoader<SData>(filename, bp_opts_, use_bp_);
    break;
  default:
    return NULL;
    break;
  }
  return NULL;
}
#endif

// instantiation

#include <string>
#include "bpmatrix.h"
#include "../common/rna.h"

#if 0
template
class Data<ProfileSequence, std::string>;
#endif

template
class Data<ProfileSequence, std::list<std::string> >;

template
class DAGBuilder<DAG::Node<> >;
