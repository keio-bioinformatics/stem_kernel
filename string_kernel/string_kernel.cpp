// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <string>
#include <vector>
#include <boost/multi_array.hpp>
#include "string_kernel.h"

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
	G1[j] += G0[i-1][j-1]*g2;
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
