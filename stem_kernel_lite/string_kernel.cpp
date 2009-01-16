// $Id$

#include <boost/multi_array.hpp>
#include <cmath>
#include "string_kernel.h"
#include "ribosum.h"
#include "data.h"
#include "dptable.h"

template < class V, class D >
StringKernel<V,D>::
StringKernel(value_type gap, value_type alpha)
  : gap_(gap), 
    subst_(boost::extents[N_RNA][N_RNA])
{
  for (uint i=0; i!=N_RNA; ++i) {
    for (uint k=0; k!=N_RNA; ++k) {
      subst_[i][k] = exp(ribosum_s[i][k] * alpha);
    }
  }
}

template < class V, class D >
StringKernel<V,D>::
StringKernel(value_type gap, value_type match, value_type mismatch)
  : gap_(gap), 
    subst_(boost::extents[N_RNA][N_RNA])
{
  for (uint i=0; i!=N_RNA; ++i) {
    for (uint k=0; k!=N_RNA; ++k) {
      subst_[i][k] = i==k ? match : mismatch;
    }
  }
}

template < class ST, class C >
static
inline
typename ST::element
subst_score(const ST& st, const C& x, const C& y)
{
  return st[index(x)][index(y)];
}

template < class ST >
static
inline
typename ST::element
subst_score(const ST& st,
	    const ProfileSequence::Column& x, const ProfileSequence::Column& y)
{
  typedef typename ST::element value_type;
  value_type v_c = 0.0;
  float n = 0;
  for (uint i=0; i!=N_RNA; ++i) {
    if (x[i]==0) continue;
    for (uint j=0; j!=N_RNA; ++j) {
      if (y[j]==0) continue;
      n += x[i]*y[j];
      v_c += st[i][j]*x[i]*y[j];
    }
  }
  return n==0 ? 1.0 : v_c / n;
}

template < class V, class D >
typename StringKernel<V,D>::value_type
StringKernel<V,D>::
operator()(const Data& xx, const Data& yy) const
{
  typedef typename Data::Seq Seq;

  const Seq& x(xx.seq);
  const Seq& y(yy.seq);
  const std::vector<float>& w_x(xx.weight);
  const std::vector<float>& w_y(yy.weight);
  bool use_weight = !w_x.empty() && !w_y.empty();
  uint sz_x = x.size();
  uint sz_y = y.size();
  DPTable<value_type> K0(sz_x+1, sz_y+1);
  DPTable<value_type> G0(sz_x+1, sz_y+1);
  std::vector<value_type> K1(sz_y+1);
  std::vector<value_type> G1(sz_y+1);

  K0.allocate(0);
  G0.allocate(0);
  K0[0][0] = G0[0][0] = 1.0;
  for (uint j=1; j!=sz_y+1; ++j) {
    K0[0][j] = 1.0;
    G0[0][j] = G0[0][j-1]*gap_;
  }

  if (use_weight) {
    for (uint i=1; i!=sz_x+1; ++i) {
      K0.allocate(i);
      G0.allocate(i);
      K0[i][0] = 1.0;
      G0[i][0] = G0[i-1][0]*gap_;
      K1[0] = G1[0] = 0.0;
      for (uint j=1; j!=sz_y+1; ++j) {
	value_type v = G0[i-1][j-1]*w_x[i-1]*w_y[j-1];
	v *= subst_score(subst_, x[i-1], y[j-1]);
	K1[j] = v + K1[j-1];
	G1[j] = v + G1[j-1]*gap_;
	K0[i][j] = K1[j] + K0[i-1][j];
	G0[i][j] = G1[j] + G0[i-1][j]*gap_;
      }
      K0.deallocate(i-1);
      G0.deallocate(i-1);
    }
  } else {
    for (uint i=1; i!=sz_x+1; ++i) {
      K0.allocate(i);
      G0.allocate(i);
      K0[i][0] = 1.0;
      G0[i][0] = G0[i-1][0]*gap_;
      K1[0] = G1[0] = 0.0;
      for (uint j=1; j!=sz_y+1; ++j) {
	value_type v = G0[i-1][j-1];
	v *= subst_score(subst_, x[i-1], y[j-1]);
	K1[j] = v + K1[j-1];
	G1[j] = v + G1[j-1]*gap_;
	K0[i][j] = K1[j] + K0[i-1][j];
	G0[i][j] = G1[j] + G0[i-1][j]*gap_;
      }
      K0.deallocate(i-1);
      G0.deallocate(i-1);
    }
  }

  return K0[sz_x][sz_y];
}

// instantiation
#include <string>
#include "data.h"
#include "../common/rna.h"

typedef double ValueType;

#if 0
template
class StringKernel<ValueType,SData>;
#endif

template
class StringKernel<ValueType,MData>;
