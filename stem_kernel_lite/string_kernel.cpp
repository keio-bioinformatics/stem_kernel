// $Id$

#include <boost/multi_array.hpp>
#include <cmath>
#include "string_kernel.h"
#include "ribosum.h"
#include "data.h"

template < class V, class D >
StringKernel<V,D>::
StringKernel(const value_type& gap, const value_type& alpha)
  : gap_(gap), alpha_(alpha),
    si_subst_(boost::extents[N_IUPAC][N_IUPAC])
{
  for (uint i=0; i!=N_IUPAC; ++i) {
    for (uint k=0; k!=N_IUPAC; ++k) {
      uint cnt=0;
      value_type val=0.0;
      for (uint a=0; a!=N_RNA; ++a) {
	if (!iupac_symbol[i][a]) continue;
	for (uint c=0; c!=N_RNA; ++c) {
	  if (!iupac_symbol[k][c]) continue;
	  val += ribosum_s[a][c];
	  cnt++;
	}
      }
      if (cnt==0)
	si_subst_[i][k] = gap_;
      else
	si_subst_[i][k] = exp(val/cnt * alpha_);
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
subst_score(const ST& st, const Col& x, const Col& y)
{
  typedef typename ST::element value_type;
  value_type v_c = 0.0;
  uint n = 0;
  for (uint i=0; i!=N_RNA; ++i) {
    if (x.cnt(i)==0) continue;
    for (uint j=0; j!=N_RNA; ++j) {
      if (y.cnt(j)==0) continue;
      n += x.cnt(i)*y.cnt(j);
      v_c += st[i][j]*x.cnt(i)*y.cnt(j);
    }
  }
  return n==0 ? 1.0 : v_c / n;
}

template < class ST, class T >
static
inline
typename ST::element
subst_score(const ST& st, const Column<T>& x, const Column<T>& y)
{
  typedef typename ST::element value_type;
  value_type v_c = 0.0;
  uint n = 0;
  for (uint i=0; i!=N_RNA; ++i) {
    if (x.cnt(i)==0) continue;
    for (uint j=0; j!=N_RNA; ++j) {
      if (y.cnt(j)==0) continue;
      n += x.cnt(i)*y.cnt(j);
      v_c += st[i][j]*x.cnt(i)*y.cnt(j);
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
  typedef boost::multi_array<value_type,2> dp_type;

  const Seq& x(xx.seq);
  const Seq& y(yy.seq);
  const std::vector<float>& w_x(xx.weight);
  const std::vector<float>& w_y(yy.weight);
  uint sz_x = x.size();
  uint sz_y = y.size();
  dp_type K0(boost::extents[sz_x+1][sz_y+1]);
  dp_type G0(boost::extents[sz_x+1][sz_y+1]);
  std::vector<value_type> K1(sz_y+1);
  std::vector<value_type> G1(sz_y+1);

  K0[0][0] = G0[0][0] = 1.0;
  for (uint i=1; i!=sz_x+1; ++i) {
    K0[i][0] = 1.0;
    G0[i][0] = G0[i-1][0]*gap_;
  }
  for (uint j=1; j!=sz_y+1; ++j) {
    K0[0][j] = 1.0;
    G0[0][j] = G0[0][j-1]*gap_;
  }

  for (uint i=1; i!=sz_x+1; ++i) {
    K1[0] = G1[0] = 0.0;
    for (uint j=1; j!=sz_y+1; ++j) {
      value_type v = G0[i-1][j-1]*w_x[i-1]*w_y[j-1];
      v *= subst_score(si_subst_, x[i-1], y[j-1]);
      K1[j] = v + K1[j-1];
      G1[j] = v + G1[j-1]*gap_;
      K0[i][j] = K1[j] + K0[i-1][j];
      G0[i][j] = G1[j] + G0[i-1][j]*gap_;
    }
  }

  return K0[sz_x][sz_y];
}

// instantiation
#include <string>
#include "data.h"
#include "../common/rna.h"

typedef double ValueType;

template
class StringKernel<ValueType,SData>;

template
class StringKernel<ValueType,MData>;
