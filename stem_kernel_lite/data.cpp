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

// helpers
template < class Seq, class BPM, class Node >
static
void
make_tree(std::vector<Node>& tree, const Seq& seq,
	  const BPM& pf, const std::vector<float>& weight, double th);

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

template < class S, class IS, class N >
Data<S,IS,N>::
Data(const IS& s, const BPMatrix::Options& opts)
  : tree(), seq(s), root(), weight(seq.size()), max_pa()
{
  BPMatrix bp(s, opts);
  fill_weight(bp, weight);
  make_tree(tree, s, bp, weight, opts.th);
  find_root(root, tree);
  find_max_parent(max_pa, tree);
}
 
template < class S, class IS, class N >
Data<S,IS,N>::
Data(const IS& s)
  : tree(), seq(s), root(), weight(), max_pa()
{
}

// helper functions
static
void
make_bp_profile(std::list<DAG::bp_freq_t>& bp_freq, uint i, uint j,
		const BPProfileMaker& maker, const BPMatrix& bpm)
{
  std::map<DAG::bp_t, float> c;
  BPMatrix::matrix_iterator mtx = bpm.matrix_begin();
  BPMatrix::idx_map_iterator idx = bpm.idx_map_begin();
  if (maker.n_seqs() == bpm.n_matrices()) {
    // weight the bp profile on (i,j) by bp prob for each seq
    std::vector<float> p(maker.n_seqs(), 0.0);
    for (uint x=0; x!=maker.n_seqs(); ++x, ++mtx, ++idx) {
      if ((*idx)[i]!=static_cast<uint>(-1) &&
	  (*idx)[j]!=static_cast<uint>(-1)) { // not GAP
	p[x] = (**mtx)((*idx)[i]+1, (*idx)[j]+1);
      }
    }
    maker.make_profile(p, i, j, c);
  } else {
    maker.make_profile(bpm(i+1, j+1), i, j, c);
  }

  std::map<DAG::bp_t, float>::iterator y;
  for (y=c.begin(); y!=c.end(); ++y) {
    bp_freq.push_back(*y);
  }
}

static
float
calc_edge_weight(const std::vector<float>& weight, const Pos& p_pos)
{
  assert(p_pos.first<p_pos.second);

  float ret = 1.0;
  for (uint i=p_pos.first+1; i!=p_pos.second; ++i)
    ret *= weight[i];

  return ret;
}

static
float
calc_edge_weight(const std::vector<float>& weight,
		 const Pos& p_pos, const Pos& c_pos)
{
  assert(p_pos.first<c_pos.first);
  assert(c_pos.second<p_pos.second);

  float ret = 1.0;
  for (uint i=p_pos.first+1; i!=c_pos.first; ++i)
    ret *= weight[i];
  for (uint i=c_pos.second+1; i!=p_pos.second; ++i)
    ret *= weight[i];

  return ret;
}

template < class BPM, class Node >
static
uint
make_tree_helper(std::vector<Node>& tree,
		 const BPProfileMaker& maker,
		 CYKTable<uint>& vt,
		 const BPM& pf,
		 const CYKTable< std::list<Pos> >& table,
		 const std::vector<float>& weight,
		 const Pos& pos)
{
  typedef typename Node::Edge Edge;
  if (vt(pos.first, pos.second)==static_cast<uint>(-1)) {
    if (pos.first==pos.second) { // leaf node
      Node leaf(pos);
      tree.push_back(leaf);
      vt(pos)=tree.size()-1;
    } else {
      const std::list<Pos>& cur = table(pos);
      if (cur.empty()) {
	std::list<DAG::bp_freq_t> bp_freq;
	make_bp_profile(bp_freq, pos.first, pos.second, maker, pf);
	Node node(pos, bp_freq, 1);
	uint ret=make_tree_helper(tree, maker, vt, pf, table,
				  weight, Pos(pos.first, pos.first));
	float e_w = calc_edge_weight(weight, pos);
	node[0] = Edge(ret, pos, e_w);
	tree.push_back(node);
	vt(pos)=tree.size()-1;
      } else {
	std::list<DAG::bp_freq_t> bp_freq;
	make_bp_profile(bp_freq, pos.first, pos.second, maker, pf);
	Node node(pos, bp_freq, cur.size());
	std::list<Pos>::const_iterator x;
	uint i;
	for (x=cur.begin(), i=0; x!=cur.end(); ++x, ++i) {
	  float e_w = calc_edge_weight(weight, pos, *x);
	  uint ret=make_tree_helper(tree, maker, vt, pf, table, weight, *x);
	  node[i] = Edge(ret, pos, *x, e_w);
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

template < class Seq, class BPM, class Node >
static
void
make_tree(std::vector<Node>& tree, const Seq& seq,
	  const BPM& pf, const std::vector<float>& weight, double th)
{
  // scan the matrix in bottom up order
  uint sz = seq_size(seq);
  CYKTable< std::list<Pos> > bp(sz);
  CYKTable< std::list<Pos> > ch(sz);
  std::vector< std::list<Pos> > head(sz);
  for (uint j=1; j!=sz; ++j)  {
    for (uint i=j-1; ; --i) {
      //std::cout <<  pf(i+1,j+1) << std::endl;
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
  BPProfileMaker maker(seq);
  CYKTable<uint> vt(sz);
  vt.fill(static_cast<uint>(-1));
  for (uint i=0; i!=sz; ++i) {
    std::list<Pos>::const_reverse_iterator x;
    for (x=head[i].rbegin(); x!=head[i].rend(); ++x) {
      make_tree_helper(tree, maker, vt, pf, bp, weight, *x);
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
