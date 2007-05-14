// $Id$

#ifndef __INC_STEM_KERNEL_H__
#define __INC_STEM_KERNEL_H__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <utility>
#include <boost/multi_array.hpp>
#include "dptable.h"

#ifdef HAVE_LIBRNA
#include "../common/pf_wrapper.h"
#endif

class NormalBasePair;
class WobbleBasePair;
#ifdef HAVE_LIBRNA
class BPMatrix;
#endif

template < class ValueType, class BPMat >
class StemKernel
{
public:
  typedef ValueType value_type;
  typedef DPTable<value_type,PoolAllocator> DPTable;

private:
  bool use_GU_;
  uint loop_;
  value_type gap_;
  value_type stack_;
  value_type subst_;
  uint band_;
  float bp_bound_;
  float ali_bound_;

public:
  StemKernel(bool use_GU, uint loop, value_type gap, value_type stack,
	     value_type subst, uint band, float ali_bound,
	     float bp_bound=1.0)
    : use_GU_(use_GU), loop_(loop), gap_(gap), stack_(stack), subst_(subst),
      band_(band), bp_bound_(bp_bound), ali_bound_(ali_bound)
  {
  }

  value_type operator()(const std::string& x, const std::string& y) const
  {
    return ali_bound_>0.0 || band_>0 ? partial_dp(x,y) : full_dp(x,y);
  }
  value_type with_subst_bp(const std::string& x, const std::string& y,
			   DPTable& dp) const;

private:
  value_type partial_dp(const std::string& x, const std::string& y) const;
  value_type full_dp(const std::string& x, const std::string& y) const;
  void alignment_constraints(const std::string& x, const std::string& y,
			     std::vector<uint>& c_low,
			     std::vector<uint>& c_high) const;
};

#endif

// Local Variables:
// mode: C++
// End:
