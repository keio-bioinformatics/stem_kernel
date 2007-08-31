// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cmath>
#include <cassert>
#include "score_table.h"
#include "../common/rna.h"
#include "ribosum.h"

template < class V, class N >
static
inline
V
simple_node_score(const V& covar, const N& x, const N& y)
{
  typedef V value_type;
  value_type v_c = 0.0;
  float n=0;
  DAG::bp_freq_iterator ix,iy;
  for (ix=x.bp_freq_begin(); ix!=x.bp_freq_end(); ++ix) {
    rna_t i = ix->first.first;
    rna_t j = ix->first.second;
    float cx = ix->second;
    for (iy=y.bp_freq_begin(); iy!=y.bp_freq_end(); ++iy) {
      rna_t k = iy->first.first;
      rna_t l = iy->first.second;
      float cy = iy->second;
      value_type v = (i!=k || j!=l) ? covar : 1.0;
      v_c += v*cx*cy;
      n += cx*cy;
    }
  }
  return n==0 ? 1.0 : v_c / n;
}

// node score for MATCH
template < class V, class D >
typename SimpleNodeScore<V,D>::value_type
SimpleNodeScore<V,D>::
node_score(const Data& xx, const Data& yy,
	   uint i, uint j) const
{
  const std::vector<Node>& x(xx.tree);
  const std::vector<Node>& y(yy.tree);
  const std::vector<float>& x_w(xx.weight);
  const std::vector<float>& y_w(yy.weight);
  float x_seqs = xx.seq.n_seqs();
  float y_seqs = yy.seq.n_seqs();

  //value_type v_s = x[i].weight()*y[j].weight();
  
  value_type v_c = 0.0;
  DAG::bp_freq_iterator ix,iy;
  float x_n = x_seqs;
  for (ix=x.bp_freq_begin(); ix!=x.bp_freq_end(); ++ix) {
    rna_t i = ix->first.first;
    rna_t j = ix->first.second;
    float cx = ix->second;
    x_n -= cx;

    float y_n = y_seqs;
    for (iy=y.bp_freq_begin(); iy!=y.bp_freq_end(); ++iy) {
      rna_t k = iy->first.first;
      rna_t l = iy->first.second;
      float cy = iy->second;
      y_n -= cy;

      value_type v = (i!=k || j!=l) ? mismatch_ : match_;
      v_c += v*cx*cy;
      n += cx*cy;
    }

    // gap for x
    
  }

  v_s *= simple_node_score(covar_, x[i], y[j]);
  return v_s;
}

// static private
template < class V, class D >
std::vector<typename SimpleEdgeScore<V,D>::value_type>
SimpleEdgeScore<V,D>::gap_vec_;

template < class V, class D >
void
SimpleEdgeScore<V,D>::
initialize(WorkArea& wa, const Data& xx, const Data& yy) const
{
  const Seq& seq_x(xx.seq);
  const Seq& seq_y(yy.seq);
  uint cur_sz=(std::max(seq_x.size(),seq_y.size()))*2;
  if (gap_vec_.empty()) {
    gap_vec_.resize(cur_sz);
    gap_vec_[0]=1.0;
    for (uint i=1; i!=cur_sz; ++i) gap_vec_[i]=gap_vec_[i-1]*gap_;
  } else if (cur_sz>gap_vec_.size()) {
    uint old_sz = gap_vec_.size();
    gap_vec_.resize(cur_sz);
    for (uint i=old_sz; i!=cur_sz; ++i) gap_vec_[i]=gap_vec_[i-1]*gap_;
  }    
}

// edge score for MATCH
template < class V, class D >
typename SimpleEdgeScore<V,D>::value_type
SimpleEdgeScore<V,D>::
edge_score(const WorkArea& wa,
	   const Data& xx, const Data& yy,
	   const Edge& ix, const Edge& iy,
	   uint i, uint j) const
{
  const std::vector<value_type>& g(gap_vec_);
  return g[ix.gaps()]*g[iy.gaps()]*ix.weight()*iy.weight();
}

// edge score for GAP
template < class V, class D >
typename SimpleEdgeScore<V,D>::value_type
SimpleEdgeScore<V,D>::
edge_score(const WorkArea& wa,
	   const Data& xx, const Edge& ix, uint i) const
{
  const std::vector<value_type>& g(gap_vec_);
  return g[ix.gaps()]*ix.weight();
}

#if 0
// edge and extented node score for MATCH where node does not use for STEM
template < class V, class D >
typename SimpleEdgeScore<V,D>::value_type
SimpleEdgeScore<V,D>::
edge_ext_score(const WorkArea& wa,
	       const Data& xx, const Data& yy,
	       const Edge& ix, const Edge& iy,
	       uint i, uint j) const
{
  const std::vector<value_type>& g(gap_vec_);
  return g[4]*g[ix.gaps()]*g[iy.gaps()];
}
#endif

template < class V, class D >
SubstNodeScore<V,D>::
SubstNodeScore(value_type gap, value_type beta)
  : gap_(gap),
    co_subst_(boost::extents[N_RNA][N_RNA][N_RNA][N_RNA]),
    beta_(beta)
{
  for (uint i=0; i!=N_RNA; ++i) {
    for (uint j=0; j!=N_RNA; ++j) {
      for (uint k=0; k!=N_RNA; ++k) {
	for (uint l=0; l!=N_RNA; ++l) {
	  co_subst_[i][j][k][l] = exp(ribosum_p[i][j][k][l] * beta_);
	}
      }
    }
  }
}

template < class ST, class N >
static
inline
typename ST::element
subst_node_score(const ST& st, const N& x, const N& y)
{
  typedef typename ST::element value_type;
  value_type v_c = 0.0;
  float n=0;
  DAG::bp_freq_iterator ix,iy;
  for (ix=x.bp_freq_begin(); ix!=x.bp_freq_end(); ++ix) {
    rna_t i = ix->first.first;
    rna_t j = ix->first.second;
    float cx = ix->second;
    for (iy=y.bp_freq_begin(); iy!=y.bp_freq_end(); ++iy) {
      rna_t k = iy->first.first;
      rna_t l = iy->first.second;
      float cy = iy->second;
      v_c += st[i][j][k][l]*cx*cy;
      n += cx*cy;
    }
  }
  return n==0 ? 1.0 : v_c / n;
}

// node score for GAP
template < class V, class D >
typename SubstNodeScore<V,D>::value_type
SubstNodeScore<V,D>::
node_score(const Data& xx, const Data& yy,
	   uint i, uint j) const
{
  const std::vector<Node>& x(xx.tree);
  const std::vector<Node>& y(yy.tree);
  value_type v_s = x[i].weight()*y[j].weight();
  v_s *= subst_node_score(co_subst_, x[i], y[j]);
  return v_s;
}


#if 0
template < class V, class D >
SubstEdgeScore<V,D>::
SubstEdgeScore(value_type gap, value_type beta)
  : gap_(gap), beta_(beta),
    si_subst_(boost::extents[N_RNA][N_RNA])
{
  for (uint i=0; i!=N_RNA; ++i) {
    for (uint k=0; k!=N_RNA; ++k) {
      si_subst_[i][k] = exp(ribosum_s[i][k] * beta_);
    }
  }
}

template < class V, class D >
void
SubstEdgeScore<V,D>::
initialize(WorkArea& wa, const Data& xx, const Data& yy) const
{
  const Seq& x(xx.seq);
  const Seq& y(yy.seq);
  std::vector<value_type>& g(wa.gap);
  g.resize(std::max(x.size(),y.size()));
  g[0]=1.0;
  for (uint i=1; i!=g.size(); ++i) g[i]=g[i-1]*gap_;
}

template < class ST, class C >
static
inline
typename ST::element
subst_score(const ST& st, const C& x, const C& y)
{
  return st[index(x)][index(y)];
}

template < class ST, class T >
static
inline
typename ST::element
subst_score(const ST& st, const Column<T>& x, const Column<T>& y)
{
  typedef typename ST::element value_type;
  const T GAP = RNASymbol<T>::GAP;
  value_type v_c = 0.0;
  uint n=0;
  for (uint s=0; s!=x.n_seqs(); ++s) {
    if (x[s]==GAP) continue;
    for (uint t=0; t!=y.n_seqs(); ++t) {
      if (y[t]==GAP) continue;
      v_c += st[index(x[s])][index(y[t])];
      ++n;
    }
  }
  return n==0 ? 1.0 : v_c / n;
}

template <class Data, class V>
static
V
calc_edge_score(const boost::multi_array<V,2>& si_subst,
		const std::vector<V>& gap,
		const Data& xx, const Data& yy,
		uint x_i, uint x_j, uint y_i, uint y_j)
{
  typedef V value_type;
  typedef typename Data::Seq Seq;
  typedef boost::multi_array<value_type,2> dp_type;
  const Seq& x(xx.seq);
  const Seq& y(yy.seq);
  const std::vector<float>& w_x(xx.weight);
  const std::vector<float>& w_y(yy.weight);
  uint sz_x = x_i<=x_j ? x_j-x_i+1 : 0;
  uint sz_y = y_i<=y_j ? y_j-y_i+1 : 0;

  if (sz_x==0) { return gap[sz_y]; }
  if (sz_y==0) { return gap[sz_x]; }

  dp_type K0(boost::extents[sz_x+1][sz_y+1]);
  dp_type G0(boost::extents[sz_x+1][sz_y+1]);
  std::vector<value_type> K1(sz_y+1);
  std::vector<value_type> G1(sz_y+1);

  K0[0][0] = G0[0][0] = 1.0;
  for (uint i=1; i!=sz_x+1; ++i) {
    K0[i][0] = 1.0;
    G0[i][0] = G0[i-1][0]*gap[1];
  }
  for (uint j=1; j!=sz_y+1; ++j) {
    K0[0][j] = 1.0;
    G0[0][j] = G0[0][j-1]*gap[1];
  }

  for (uint i=1; i!=sz_x+1; ++i) {
    K1[0] = G1[0] = 0.0;
    uint ii=x_i+i;
    for (uint j=1; j!=sz_y+1; ++j) {
      uint jj=y_i+j;
      value_type v = G0[i-1][j-1];
      v *= subst_score(si_subst,x[ii-1],y[jj-1]);
      v *= w_x[ii-1]*w_y[jj-1];
      K1[j] = v + K1[j-1];
      G1[j] = v + G1[j-1]*gap[1];
      K0[i][j] = K1[j] + K0[i-1][j];
      G0[i][j] = G1[j] + G0[i-1][j]*gap[1];
    }
  }

  return K0[sz_x][sz_y];
}

// edge score for MATCH
template < class V, class D >
typename SubstEdgeScore<V,D>::value_type
SubstEdgeScore<V,D>::
edge_score(const WorkArea& wa,
	   const Data& xx, const Data& yy,
	   const Edge& ix, const Edge& iy,
	   uint i, uint j) const
{
  value_type s=0.0;
  if (ix.p_pos().first==ix.c_pos().first ||
      iy.p_pos().first==iy.c_pos().first) { // loop region
    s = calc_edge_score(si_subst_, wa.gap, xx, yy,
	 		ix.p_pos().first+1, ix.p_pos().second-1,
                        iy.p_pos().first+1, iy.p_pos().second-1);
  } else {			// buldge region
    s  =  calc_edge_score(si_subst_, wa.gap, xx, yy,
			  ix.p_pos().first+1, ix.c_pos().first-1,
			  iy.p_pos().first+1, iy.c_pos().first-1);
    s += calc_edge_score(si_subst_, wa.gap, xx, yy,
			 ix.c_pos().second+1, ix.p_pos().second-1,
			 iy.c_pos().second+1, iy.p_pos().second-1);
  }
  return s;
}

// edge score for GAP
template < class V, class D >
typename SubstEdgeScore<V,D>::value_type
SubstEdgeScore<V,D>::
edge_score(const WorkArea& wa,
	   const Data& xx, const Edge& ix, uint i) const
{
  const std::vector<value_type>& g(wa.gap);
  return g[ix.gaps()];
}

// edge and extented node score for MATCH where node does not use for STEM
template < class V, class D >
typename SubstEdgeScore<V,D>::value_type
SubstEdgeScore<V,D>::
edge_ext_score(const WorkArea& wa,
	       const Data& xx, const Data& yy,
	       const Edge& ix, const Edge& iy,
	       uint i, uint j) const
{
  value_type s=0.0;
  if (ix.p_pos().first==ix.c_pos().first ||
      iy.p_pos().first==iy.c_pos().first) { // loop region
    s = calc_edge_score(si_subst_, wa.gap, xx, yy,
	 		ix.p_pos().first, ix.p_pos().second,
                        iy.p_pos().first, iy.p_pos().second);
  } else {				    // buldge
    s  =  calc_edge_score(si_subst_, wa.gap, xx, yy,
			  ix.p_pos().first, ix.c_pos().first-1,
			  iy.p_pos().first, iy.c_pos().first-1);
    s += calc_edge_score(si_subst_, wa.gap, xx, yy,
			 ix.c_pos().second+1, ix.p_pos().second,
			 iy.c_pos().second+1, iy.p_pos().second);
  }
  return s;
}
#endif

// instantiation
#include "data.h"

typedef double ValueType;

template
class SimpleNodeScore<ValueType,SData>;

template
class SimpleEdgeScore<ValueType,SData>;

template
class SubstNodeScore<ValueType,SData>;

#if 0
template
class SubstEdgeScore<ValueType,SData>;
#endif

template
class SimpleNodeScore<ValueType,MData>;

template
class SimpleEdgeScore<ValueType,MData>;

template
class SubstNodeScore<ValueType,MData>;

#if 0
template
class SubstEdgeScore<ValueType,MData>;
#endif
