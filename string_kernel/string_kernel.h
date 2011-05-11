// $Id$

#ifndef __INC_STRING_KERNEL_H__
#define __INC_STRING_KERNEL_H__

#include <string>

typedef unsigned int uint;

template < class ValueType >
class StringKernel
{
public:
  typedef ValueType value_type;

private:
  value_type gap_;

public:
  StringKernel(value_type gap=1) : gap_(gap) { }
  
  value_type operator()(const std::string& x, const std::string& y) const;
};

#endif
