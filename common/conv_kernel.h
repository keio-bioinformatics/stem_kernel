// $Id:$

#ifndef __INC_CONV_KERNEL_H__
#define __INC_CONV_KERNEL_H__

#include <cmath>

// convolution kernels


// linear trasformation
template < class K1 >
class LTKernel 
{
public:
  typedef typename K1::value_type value_type;
  typedef typename K1::Data Data;
  
public:
  LTKernel(const K1& k1, value_type a, value_type b)
    : k1_(k1), a_(a), b_(b)
  {
  }

  value_type operator()(const Data& x, const Data& y) const
  {
    return a_*k1_(x,y)+b_;
  }

private:
  K1 k1_;
  value_type a_;
  value_type b_;
};

// add two kernels
template < class K1, class K2 >
class AddKernel
{
public:
  typedef typename K1::value_type value_type;
  typedef typename K1::Data Data;
  
public:
  AddKernel(const K1& k1, const K2& k2) : k1_(k1), k2_(k2)
  {
  }

  value_type operator()(const Data& x, const Data& y) const
  {
    return k1_(x,y)+k2_(x,y);
  }

private:
  K1 k1_;
  K2 k2_;
};

template < class K1 >
class ExpKernel
{
public:
  typedef typename K1::value_type value_type;
  typedef typename K1::Data Data;
  
public:
  ExpKernel(const K1& k1) : k1_(k1)
  {
  }

  value_type operator()(const Data& x, const Data& y) const
  {
    return exp(k1_(x,y));
  }

private:
  K1 k1_;
};

// note: LogKernel violates positive semi-definite
template < class K1 >
class LogKernel
{
public:
  typedef typename K1::value_type value_type;
  typedef typename K1::Data Data;
  
public:
  LogKernel(const K1& k1) : k1_(k1)
  {
  }

  value_type operator()(const Data& x, const Data& y) const
  {
    return log(k1_(x,y));
  }

private:
  K1 k1_;
};

#endif	// __INC_CONV_KERNEL_H__
