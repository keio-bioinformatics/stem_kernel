// $Id$

#ifndef __INC_STRING_KERNEL_H__
#define __INC_STRING_KERNEL_H__

#include <string>
#include <vector>
#include <cmath>
#include "data.h"

template < class ValueType, class Data >
class BPLAKernel
{
public:
  typedef ValueType value_type;
  
public:
  BPLAKernel(value_type gap=1,value_type ext=1,
	     value_type alpha=1,value_type beta=1) 
    : gap_(gap), ext_(ext), alpha_(alpha), beta_(beta), 
      beta_gap_(exp(beta*gap)), beta_ext_(exp(beta*ext))
  {}

  value_type operator()(const Data& xx, const Data& yy) const;

private:
  value_type gap_;
  value_type ext_;
  value_type alpha_;
  value_type beta_;
  value_type beta_gap_;
  value_type beta_ext_;
};

#endif

// Local Variables:
// mode: C++
// End:
