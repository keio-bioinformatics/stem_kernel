// $Id: string_kernel.cpp 81 2006-09-15 11:41:21Z satoken $

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <boost/multi_array.hpp>
#include "string_kernel.h"
#include "../common/pf_wrapper.h"

double RIBOSUM[5][5]={
  /* A    C      G      U       N */
  {2.22,  -1.86, -1.46, -1.39,  0}, //A
  {-1.86, 1.16,  -2.48, -1.05,  0}, //C
  {-1.46, -2.48, 1.03,  -1.74,  0}, //G
  {-1.39, -1.05, -1.74, -1.65,  0}, //U
  {    0,     0,     0,     0,  0}, //N
};

template < class ValueType >
int
StringKernel<ValueType>::
char2rna(char c) const
{
  switch(c){
  case 'a':
    return 0;
    break;
  case 'c':
    return 1;
    break;
  case 'g':
    return 2;
    break;
  case 'u':
    return 3;
    break;
  default:
    return 4;
    break;
  }
}
/*
char rna2char(uint i)
{
  switch(i){
  case 0:
    return 'a';
    break;
  case 1:
    return 'c';
    break;
  case 2:
    return 'g';
    break;
  case 3:
    return 'u';
    break;
  defalut:
    return 'N';
    break;
  }
}
*/


template < class ValueType >
ValueType
StringKernel<ValueType>::
score(const Seq& x, const Seq& y, uint i, uint j) const
{  
  int base1 = char2rna(x.seq[i]);
  int base2 = char2rna(y.seq[j]);
  value_type p0x = 1 - (x.p_r[i] + x.p_l[i]);
  value_type p0y = 1 - (y.p_r[j] + y.p_l[j]);
  if (p0x < 0) p0x = 0;
  if (p0y < 0) p0y = 0;
  return alpha_ * (sqrt(x.p_r[i] * y.p_r[j]) + sqrt(x.p_l[i] * y.p_l[j]))
    + sqrt(p0x * p0y) * RIBOSUM[base1][base2];
}

template < class ValueType >
ValueType
StringKernel<ValueType>::
operator()(const Seq& x, const Seq& y) const
{
  // const value_type& g=gap_;
  // const value_type& g2 = g*g;
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
    //std::cout << rna2char(y.seq[i-1]) <<"\t\t"<< y.p_l[i-1]<<"\t\t"<< y.p_r[i-1]<< "\t\t"<< y.p_un[i-1]<< std::endl;
    for (uint j=1; j != y.size()+1; ++j) {
      M[i][j] = exp(beta_ * score(x, y, i-1, j-1))
	* (1 + X[i-1][j-1] + Y[i-1][j-1] + M[i-1][j-1]);
      X[i][j] = exp(beta_ * gap_) * M[i-1][j] 
	+ exp(beta_ * ext_) * X[i-1][j] ;
      
      Y[i][j] = exp(beta_ * gap_) * (M[i][j-1] + X[i][j-1]) 
	+ exp(beta_ * ext_) * Y[i][j-1] ;

      X2[i][j]= M[i-1][j] + X2[i-1][j];
      Y2[i][j]= M[i][j-1] + X2[i][j-1] + Y2[i][j-1];
    }
  }

  //std::cout << std::endl;

  /*
  for (uint i=1; i < x.size()+1; ++i) {
    for (uint j=1; j < y.size()+1; ++j) {
      std::cout << X[i][j] << " " ;
    }
    std::cout << std::endl;
  }
  */ 
  return 1+X2[x.size()][y.size()]
    + Y2[x.size()][y.size()] + M[x.size()][y.size()];
}

template < class ValueType >
ValueType
StringKernel<ValueType>::
operator()(const std::string& x, const std::string& y) const
{
  const value_type& g=gap_;
  value_type g2=g*g;
  typedef boost::multi_array<value_type,2> dp_type;
#if 1
 
  dp_type K0(boost::extents[x.size()+1][y.size()+1]);
  dp_type G0(boost::extents[x.size()+1][y.size()+1]);
  std::vector<value_type> K1(y.size()+1);
  std::vector<value_type> G1(y.size()+1);

  K0[0][0]=G0[0][0]=1.0;
  for (uint i=1; i!=x.size()+1; ++i) {
    K0[i][0]=1.0;
    G0[i][0]=G0[i-1][0]*g;
  }
  for (uint j=1; j!=y.size()+1; ++j) {
    K0[0][j]=1.0;
    G0[0][j]=G0[0][j-1]*g;
  }

  for (uint i=1; i!=x.size()+1; ++i) {
    K1[0]=G1[0]=0.0;
    for (uint j=1; j!=y.size()+1; ++j) {
      K1[j] = K1[j-1];
      G1[j] = G1[j-1]*g;
      if (x[i-1]==y[j-1]) {
	K1[j] += G0[i-1][j-1]*g2;
	K1[j] += G0[i-1][j-1]*g2;
      }
      K0[i][j] = K0[i-1][j] + K1[j];
      G0[i][j] = G0[i-1][j]*g + G1[j];
    }
  }

  return K0[x.size()][y.size()];
  
#else
  dp_type K0(boost::extents[x.size()+1][y.size()+1]);
  dp_type G0(boost::extents[x.size()+1][y.size()+1]);
  dp_type K1(boost::extents[x.size()+1][y.size()+1]);
  dp_type G1(boost::extents[x.size()+1][y.size()+1]);

  K0[0][0]=G0[0][0]=1.0;
  for (uint i=1; i!=x.size()+1; ++i) {
    K0[i][0]=1.0;
    G0[i][0]=G0[i-1][0]*g;
  }
  for (uint j=1; j!=y.size()+1; ++j) {
    K0[0][j]=1.0;
    G0[0][j]=G0[0][j-1]*g;
  }

  for (uint i=1; i!=x.size()+1; ++i) {
    K1[i][0]=G1[i][0]=0.0;
    for (uint j=1; j!=y.size()+1; ++j) {
      K1[i][j] = K1[i][j-1];
      G1[i][j] = G1[i][j-1]*g;
      if (x[i-1]==y[j-1]) {
	K1[i][j] += G0[i-1][j-1]*g2;
	G1[i][j] += G0[i-1][j-1]*g2;
      }
      K0[i][j] = K0[i-1][j] + K1[i][j];
      G0[i][j] = G0[i-1][j]*g + G1[i][j];
    }
  }

  return K0[x.size()][y.size()];
#endif
}

template class StringKernel<double>;
