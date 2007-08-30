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
  StringKernel(value_type gap, value_type alpha);
  StringKernel(value_type gap, value_type match, value_type mismatch);

  value_type operator()(const Data& xx, const Data& yy) const;

  template < class Seq >
  value_type operator()(const Seq& x, const Seq& y) const
  {
    return (*this)(Data(x), Data(y));
  }

private:
  value_type gap_;
  boost::multi_array<value_type,2> subst_;
};

#endif	// __INC_STRING_KERNEL_H__

// Local Variables:
// mode: C++
// End:
