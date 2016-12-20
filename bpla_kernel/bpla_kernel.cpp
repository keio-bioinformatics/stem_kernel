// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <boost/multi_array.hpp>
#include "bpla_kernel.h"

template < class ValueType, class Data >
struct LAScore
{
  LAScore(const boost::multi_array<ValueType,2>& score_table)
    : score_table_(score_table)
  { }

  inline
  ValueType
  operator()(const Data& x, const Data& y, uint i, uint j) const
  {
    assert(x.seq[0].size()==y.seq[0].size());
    assert(score_table_.size()==x.seq[0].size()-1);
    assert(score_table_[0].size()==x.seq[0].size()-1);

    ValueType v = 0.0;
    float n = 0;
    for (uint k=0; k!=x.seq[i].size()-1; ++k) {
      if (x.seq[i][k]==0) continue;
      for (uint l=0; l!=y.seq[j].size()-1; ++l) {
	if (y.seq[j][l]==0) continue;
	n += x.seq[i][k]*y.seq[j][l];
	v += score_table_[k][l]*x.seq[i][k]*y.seq[j][l];
      }
    }
    return n==0 ? 0.0 : v/n;
  }

  const boost::multi_array<ValueType,2>& score_table_;
};

template < class ValueType, class Data >
struct BPLAScore : public LAScore<ValueType,Data>
{
  BPLAScore(const boost::multi_array<ValueType,2>& score_table, ValueType alpha)
    : LAScore<ValueType,Data>(score_table), alpha_(alpha) { }

  inline
  ValueType
  operator()(const Data& x, const Data& y, uint i, uint j) const
  {
    return alpha_ * (x.p_right[i]*y.p_right[j] + x.p_left[i] * y.p_left[j])
      + x.p_unpair[i]*y.p_unpair[j] * LAScore<ValueType,Data>::operator()(x,y,i,j);
  }

  ValueType alpha_;
};

template < class Data, class ValueType, class Score >
inline
ValueType
local_alignment_exp(const Data& x, const Data& y,
                    ValueType beta, ValueType gap, ValueType ext,
                    Score score)
{
  ValueType beta_gap = exp(beta*gap);
  ValueType beta_ext = exp(beta*ext);
  typedef boost::multi_array<ValueType,2> dp_type;
  dp_type M(boost::extents[x.size()+1][y.size()+1]);
  dp_type X(boost::extents[x.size()+1][y.size()+1]);
  dp_type Y(boost::extents[x.size()+1][y.size()+1]);
  dp_type X2(boost::extents[x.size()+1][y.size()+1]);
  dp_type Y2(boost::extents[x.size()+1][y.size()+1]);

  //initialize  
  for (uint i=0; i!=x.size()+1; ++i) {
    M[i][0]=0;
    X[i][0]=0;
    Y[i][0]=0;
    X2[i][0]=0;
    Y2[i][0]=0;
  }

  for (uint j=0; j!=y.size()+1; ++j) {
    M[0][j]=0;
    X[0][j]=0;
    Y[0][j]=0;
    X2[0][j]=0;
    Y2[0][j]=0;
  }

  //calculate
  for (uint i=1; i != x.size()+1; ++i) {
    for (uint j=1; j != y.size()+1; ++j) {
      M[i][j] = exp(beta * score(x, y, i-1, j-1))
	* (1 + X[i-1][j-1] + Y[i-1][j-1] + M[i-1][j-1]);
      X[i][j] = beta_gap * M[i-1][j] 
	+ beta_ext * X[i-1][j] ;
      
      Y[i][j] = beta_gap * (M[i][j-1] + X[i][j-1]) 
	+ beta_ext * Y[i][j-1] ;

      X2[i][j]= M[i-1][j] + X2[i-1][j];
      Y2[i][j]= M[i][j-1] + X2[i][j-1] + Y2[i][j-1];
    }
  }

  return 1+X2[x.size()][y.size()]
    + Y2[x.size()][y.size()] + M[x.size()][y.size()];
}

template < class Data, class ValueType, class Score >
inline
ValueType
local_alignment_max(const Data& x, const Data& y,
                    ValueType gap, ValueType ext,
                    Score score)
{
  typedef boost::multi_array<ValueType,2> dp_type;
  dp_type M(boost::extents[x.size()+1][y.size()+1]);
  dp_type X(boost::extents[x.size()+1][y.size()+1]);
  dp_type Y(boost::extents[x.size()+1][y.size()+1]);
  ValueType Mmax = 0;

  //initialize  
  for (uint i=0; i!=x.size()+1; ++i) {
    M[i][0]=0;
    X[i][0]=0;
    Y[i][0]=0;
  }

  for (uint j=0; j!=y.size()+1; ++j) {
    M[0][j]=0;
    X[0][j]=0;
    Y[0][j]=0;
  }

  //calculate
  for (uint i=1; i != x.size()+1; ++i) {
    for (uint j=1; j != y.size()+1; ++j) {
      M[i][j] = std::max(0.0, M[i-1][j-1]);
      M[i][j] = std::max(M[i][j], X[i-1][j-1]);
      M[i][j] = std::max(M[i][j], Y[i-1][j-1]);
      M[i][j] += score(x, y, i-1, j-1);
      Mmax = std::max(Mmax, M[i][j]);
      X[i][j] = std::max(M[i-1][j]+gap, X[i-1][j]+ext);
      Y[i][j] = std::max(std::max(M[i][j-1]+gap, X[i][j-1]+gap), Y[i][j-1]+ext);
    }
  }

  return Mmax;
}

template < class ValueType, class Data >
ValueType
BPLAKernel<ValueType,Data>::
operator()(const Data& x, const Data& y) const
{
  if (SW_)
    if (noBP_)
      return local_alignment_max(x, y, gap_, ext_, LAScore<ValueType,Data>(score_table_));
    else
      return local_alignment_max(x, y, gap_, ext_, BPLAScore<ValueType,Data>(score_table_, alpha_));
  else
    if (noBP_)
      return local_alignment_exp(x, y, beta_, gap_, ext_, LAScore<ValueType,Data>(score_table_));
    else
      return local_alignment_exp(x, y, beta_, gap_, ext_, BPLAScore<ValueType,Data>(score_table_, alpha_));
}


enum { M=0, IX=1, IY=2, LX=3, LY=4, RX=5, RY=6, N=7 };

template < class ValueType, class Data, class Table >
ValueType
BPLA_Forward(const Data& x, const Data& y,
             const boost::multi_array<ValueType,2>& score_table,
	     const std::vector<double>& param, Table& T)
{
  assert(T.size()==N);
  LAScore<ValueType, Data> la_score(score_table);
  ValueType alpha = param[0];
  ValueType beta = param[1];
  ValueType gap = param[2];
  ValueType ext = param[3];
  ValueType beta_gap = exp(beta*gap);
  ValueType beta_ext = exp(beta*ext);

  std::fill(T.data(), T.data()+T.num_elements(), 0.0);
  T[M][0][0]=1;
  T[LX][0][0]=1;
  T[LY][0][0]=1;
  
  for (uint i=1; i!=x.size()+1; ++i) {
    T[LX][i][0]+=T[LX][i-1][0];
  }

  for (uint j=1; j!=y.size()+1; ++j) {
    T[LY][0][j]+=T[LY][0][j-1];
  }

  for (uint i=1; i!= x.size()+1; ++i) {
    for (uint j=1; j!= y.size()+1; ++j) {
      ValueType s =
	alpha * (x.p_right[i-1]*y.p_right[j-1] + x.p_left[i-1] * y.p_left[j-1])
	+ x.p_unpair[i-1]*y.p_unpair[j-1] * la_score(x, y, i-1, j-1);
      ValueType beta_s = exp(beta * s);
      T[M][i][j] += beta_s*T[M][i-1][j-1];
      T[M][i][j] += beta_s*T[IX][i-1][j-1];
      T[M][i][j] += beta_s*T[IY][i-1][j-1];
      T[M][i][j] += beta_s*T[LX][i-1][j-1];
      T[M][i][j] += beta_s*T[LY][i-1][j-1];

      T[IX][i][j] += beta_gap*T[M][i-1][j];
      T[IX][i][j] += beta_ext*T[IX][i-1][j];

      T[IY][i][j] += beta_gap*T[M][i][j-1];
      T[IY][i][j] += beta_gap*T[IX][i][j-1];
      T[IY][i][j] += beta_ext*T[IY][i][j-1];

      T[LX][i][j] += T[LX][i-1][0];

      T[LY][i][j] += T[LX][i][j-1];
      T[LY][i][j] += T[LY][i][j-1];

      T[RX][i][j] += T[M][i-1][j];
      T[RX][i][j] += T[RX][i-1][j];

      T[RY][i][j] += T[M][i][j-1];
      T[RY][i][j] += T[RX][i][j-1];
      T[RY][i][j] += T[RY][i][j-1];
    }
  }

  return 1+T[M][x.size()][y.size()]
    +T[RX][x.size()][y.size()]+T[RY][x.size()][y.size()];
}

template < class ValueType, class Data, class Table >
ValueType
BPLA_Backward(const Data& x, const Data& y,
              const boost::multi_array<ValueType,2>& score_table,
	      const std::vector<double>& param, Table& T)
{
  assert(T.size()==N);
  LAScore<ValueType, Data> la_score(score_table);
  ValueType alpha = param[0];
  ValueType beta = param[1];
  ValueType gap = param[2];
  ValueType ext = param[3];
  ValueType beta_gap = exp(beta*gap);
  ValueType beta_ext = exp(beta*ext);

  std::fill(T.data(), T.data()+T.num_elements(), 0.0);
  T[M][x.size()][y.size()]=1;
  T[RX][x.size()][y.size()]=1;
  T[RY][x.size()][y.size()]=1;
  
  for (uint i=x.size(); i!=0; --i) {
    for (uint j=y.size(); j!=0; --j) {
      ValueType s =
	alpha * (x.p_right[i-1]*y.p_right[j-1] + x.p_left[i-1] * y.p_left[j-1])
	+ x.p_unpair[i-1]*y.p_unpair[j-1] * la_score(x, y, i-1, j-1);
      ValueType beta_s = exp(beta * s);
      T[M][i-1][j-1] += beta_s*T[M][i][j] ;
      T[IX][i-1][j-1] += beta_s*T[M][i][j];
      T[IY][i-1][j-1] += beta_s*T[M][i][j];
      T[LX][i-1][j-1] += beta_s*T[M][i][j];
      T[LY][i-1][j-1] += beta_s*T[M][i][j];

      T[M][i-1][j] += beta_gap*T[IX][i][j];
      T[IX][i-1][j] += beta_ext*T[IX][i][j];

      T[M][i][j-1] += beta_gap*T[IY][i][j];
      T[IX][i][j-1] += beta_gap*T[IY][i][j];
      T[IY][i][j-1] += beta_ext*T[IY][i][j];

      T[LX][i-1][0] += T[LX][i][j];

      T[LX][i][j-1] += T[LY][i][j];
      T[LY][i][j-1] += T[LY][i][j];

      T[M][i-1][j] += T[RX][i][j];
      T[RX][i-1][j] += T[RX][i][j];

      T[M][i][j-1] += T[RY][i][j];
      T[RX][i][j-1] += T[RY][i][j];
      T[RY][i][j-1] += T[RY][i][j];
    }
  }

  for (uint i=x.size(); i!=0; --i) {
    T[LX][i-1][0]+=T[LX][i][0];
  }

  for (uint j=y.size(); j!=0; --j) {
    T[LY][0][j-1]+=T[LY][0][j];
  }

  return 1+T[M][0][0]+T[LX][0][0]+T[LY][0][0];
}

template < class V >
inline
void
update_alpha_beta(double& d_a, double& d_b, V a, V b, V w_pair, V w_unpair, V v)
{
  d_a += b*w_pair*v;
  d_b += (a*w_pair+w_unpair)*v;
}

template < class V >
inline
void
update_beta_gap_ext(double& d_b, double& d_ge, V b, V ge, V v)
{
  d_b  += ge*v;
  d_ge +=  b*v;
}
  
template < class ValueType, class Data, class Table >
ValueType
BPLA_ForwardBackword(const Data& x, const Data& y,
                     const boost::multi_array<ValueType,2>& score_table,
		     const std::vector<double>& param,
		     const Table& F, const Table& B,
		     std::vector<double>& d)
{
  assert(F.size()==N);
  assert(B.size()==N);
  LAScore<ValueType, Data> la_score(score_table);
  ValueType alpha = param[0];
  ValueType beta = param[1];
  ValueType gap = param[2];
  ValueType ext = param[3];
  ValueType beta_gap = exp(beta*gap);
  ValueType beta_ext = exp(beta*ext);

  assert(d.size()==param.size());
  double& d_alpha = d[0];
  double& d_beta = d[1];
  double& d_gap = d[2];
  double& d_ext = d[3];

  for (uint i=1; i!= x.size()+1; ++i) {
    for (uint j=1; j!= y.size()+1; ++j) {
      ValueType w_pair=x.p_right[i-1]*y.p_right[j-1] + x.p_left[i-1] * y.p_left[j-1];
      ValueType w_unpair=x.p_unpair[i-1]*y.p_unpair[j-1] * la_score(x, y, i-1, j-1);
      ValueType beta_s = exp(beta * (alpha * w_pair + w_unpair));
      
      update_alpha_beta(d_alpha, d_beta, alpha, beta, w_pair, w_unpair,
			F[M][i-1][j-1]*beta_s*B[M][i][j]);
      update_alpha_beta(d_alpha, d_beta, alpha, beta, w_pair, w_unpair,
			F[IX][i-1][j-1]*beta_s*B[M][i][j]);
      update_alpha_beta(d_alpha, d_beta, alpha, beta, w_pair, w_unpair,
			F[IY][i-1][j-1]*beta_s*B[M][i][j]);
      update_alpha_beta(d_alpha, d_beta, alpha, beta, w_pair, w_unpair,
			F[LX][i-1][j-1]*beta_s*B[M][i][j]);
      update_alpha_beta(d_alpha, d_beta, alpha, beta, w_pair, w_unpair,
			F[LY][i-1][j-1]*beta_s*B[M][i][j]);

      update_beta_gap_ext(d_beta, d_gap, beta, gap, 
			  F[M][i-1][j]*beta_gap*B[IX][i][j]);
      update_beta_gap_ext(d_beta, d_ext, beta, ext, 
			  F[IX][i-1][j]*beta_ext*B[IX][i][j]);

      update_beta_gap_ext(d_beta, d_gap, beta, gap, 
			  F[M][i][j-1]*beta_gap*B[IY][i][j]);
      update_beta_gap_ext(d_beta, d_gap, beta, gap, 
			  F[IX][i][j-1]*beta_gap*B[IY][i][j]);
      update_beta_gap_ext(d_beta, d_ext, beta, ext, 
			  F[IY][i][j-1]*beta_ext*B[IY][i][j]);
    }
  }

  return 1+F[M][x.size()][y.size()]
    +F[RX][x.size()][y.size()]+F[RY][x.size()][y.size()];
}

// static
template < class ValueType, class Data >
ValueType
BPLAKernel<ValueType,Data>::
compute_gradients(const Data& x, const Data& y,
                  const boost::multi_array<value_type,2>& score_table,
		  const std::vector<double>& param,
		  std::vector<double>& d)
{
  std::fill(d.begin(), d.end(), 0.0);
  typedef boost::multi_array<ValueType,3> dp_type;
  dp_type F(boost::extents[N][x.size()+1][y.size()+1]);
  dp_type B(boost::extents[N][x.size()+1][y.size()+1]);
  BPLA_Forward<ValueType,Data,dp_type>(x, y, score_table, param, F);
  BPLA_Backward<ValueType,Data,dp_type>(x, y, score_table, param, B);
  return BPLA_ForwardBackword<ValueType,Data,dp_type>(x, y, score_table, param, F, B, d);
}

// instatiation
template class BPLAKernel<double,MData>;

template class BPLAKernel<double,AAMData>;
