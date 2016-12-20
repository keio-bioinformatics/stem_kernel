// $Id$

#ifndef __INC_STRING_KERNEL_H__
#define __INC_STRING_KERNEL_H__

#include <string>
#include <vector>
#include <cmath>
#include <boost/multi_array.hpp>
#include "data.h"
#include "../optimizer/gradient.h"

template < class ValueType, class DataT >
class BPLAKernel
{
public:
  typedef ValueType value_type;
  typedef DataT Data;
  
public:
  BPLAKernel(const boost::multi_array<value_type,2>& score_table,
             bool noBP, bool SW=false, value_type gap=1,value_type ext=1,
	     value_type alpha=1,value_type beta=1) 
    : score_table_(score_table),
      noBP_(noBP), SW_(SW), gap_(gap), ext_(ext), alpha_(alpha), beta_(beta)
  {}

  value_type operator()(const Data& xx, const Data& yy) const;

  static
  value_type compute_gradients(const Data& xx, const Data& yy,
                               const boost::multi_array<value_type,2>& score_table,
			       const std::vector<double>& param,
			       std::vector<double>& d);

private:
  const boost::multi_array<value_type,2>& score_table_;
  bool noBP_;
  bool SW_;
  value_type gap_;
  value_type ext_;
  value_type alpha_;
  value_type beta_;
};

#endif

// Local Variables:
// mode: C++
// End:
