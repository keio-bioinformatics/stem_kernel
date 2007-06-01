// $Id$

#ifndef __INC_SS_KERNEL_H__
#define __INC_SS_KERNEL_H__ 

#include "stem_kernel.h"
#include "string_kernel.h"

template < class ST, class D >
class StemStrKernel
{
public:
  typedef ST ScoreTable;
  typedef typename ScoreTable::value_type value_type;
  typedef D Data;

public:
  StemStrKernel(const ScoreTable& st, value_type gap, value_type alpha,
		uint len_band=0)
    : stem_(st, len_band), str_(gap,alpha)
  {
  }

  value_type operator()(const Data& xx, const Data& yy) const
  {
    return stem_(xx,yy)+str_(xx,yy);
  }

  template < class Seq >
  value_type operator()(const Seq& x, const Seq& y) const
  {
    return (*this)(Data(x), Data(y));
  }

private:
  StemKernel<ST,D> stem_;
  StringKernel<value_type,D> str_;
};

#endif	//  __INC_SS_KERNEL_H__

// Local Variables:
// mode: C++
// End:
