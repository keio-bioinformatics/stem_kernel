// $Id$

#ifndef __INC_STRING_KERNEL_H__
#define __INC_STRING_KERNEL_H__

#include <boost/multi_array.hpp>

template < class V, class D >
class StringKernel
{
public:
  typedef V value_type;
  typedef D Data;

public:
  StringKernel(const value_type& gap, const value_type& alpha);

  value_type operator()(const Data& xx, const Data& yy) const;

  template < class Seq >
  value_type operator()(const Seq& x, const Seq& y) const
  {
    return (*this)(Data(x), Data(y));
  }

private:
  value_type gap_;
  value_type alpha_;
  boost::multi_array<value_type,2> si_subst_;
};

#endif	// __INC_STRING_KERNEL_H__

// Local Variables:
// mode: C++
// End:
