// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <boost/multi_array.hpp>
#include "bpla_kernel.h"
#include "../common/pf_wrapper.h"

double RIBOSUM[5][5]={
  /* A    C      G      U       N */
  { 2.22, -1.86, -1.46, -1.39,  0}, //A
  {-1.86,  1.16, -2.48, -1.05,  0}, //C
  {-1.46, -2.48,  1.03, -1.74,  0}, //G
  {-1.39, -1.05, -1.74,  1.65,  0}, //U
  {    0,     0,     0,     0,  0}, //N
};

double RIBOSUM2[5][5]={
  /* A    C      G      U       N */
  { 5.846613, -1.860000, -1.460000, -1.390000, 0},
  {-1.860000,  4.786613, -2.480000, -1.050000, 0},
  {-1.460000, -2.480000,  4.656613, -1.740000, 0},
  {-1.390000, -1.050000, -1.740000,  5.276613, 0},
  { 0.000000,  0.000000,  0.000000,  0.000000, 0},
};

double CRF[5][5]={
  /* A    C      G      U       N */
  { 0.67, -1.11, -1.02, -1.17,  0.51}, //A
  {-1.11, -0.44, -1.96, -1.30,  0.19}, //C
  {-1.02, -1.96, -0.52, -2.03,  0.22}, //G
  {-1.17, -1.30, -2.03, -0.68,  0.25}, //U
  { 0.51,  0.19,  0.22,   0.25,    0}, //N
};

double SIMPLE[5][5]={
/* A C G U N */
  {1,0,0,0,0},//A
  {0,1,0,0,0},//C
  {0,0,1,0,0},//G
  {0,0,0,1,0},//U
  {0,0,0,0,0} //N
};


template < class ValueType, class Data >
static
ValueType
score(const Data& x, const Data& y, uint i, uint j, ValueType alpha) 
{
  ValueType v = 0.0;
  float n = 0;
  for (uint k=0; k!=N_RNA; ++k) {
    if (x.seq[i][k]==0) continue;
    for (uint l=0; l!=N_RNA; ++l) {
      if (y.seq[j][l]==0) continue;
      n += x.seq[i][k]*y.seq[j][l];
      v += RIBOSUM2[k][l]*x.seq[i][k]*y.seq[j][l];
    }
  }
  v = n==0 ? 0.0 : v/n;
  
  return alpha * (sqrt(x.p_right[i] * y.p_right[j] )+ sqrt(x.p_left[i] * y.p_left[j]))
    + sqrt(x.p_unpair[i] * y.p_unpair[j]) * v;
}


template < class ValueType, class Data >
ValueType
BPLAKernel<ValueType,Data>::
operator()(const Data& x, const Data& y) const
{
  typedef boost::multi_array<value_type,2> dp_type;
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
      M[i][j] = exp(beta_ * score(x, y, i-1, j-1, alpha_))
	* (1 + X[i-1][j-1] + Y[i-1][j-1] + M[i-1][j-1]);
      X[i][j] = beta_gap_ * M[i-1][j] 
	+ beta_ext_ * X[i-1][j] ;
      
      Y[i][j] = beta_gap_ * (M[i][j-1] + X[i][j-1]) 
	+ beta_ext_ * Y[i][j-1] ;

      X2[i][j]= M[i-1][j] + X2[i-1][j];
      Y2[i][j]= M[i][j-1] + X2[i][j-1] + Y2[i][j-1];
    }
  }

  return 1+X2[x.size()][y.size()]
    + Y2[x.size()][y.size()] + M[x.size()][y.size()];
}

template class BPLAKernel<double,MData>;
