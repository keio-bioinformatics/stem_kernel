// $Id: pf_wrapper.h 48 2006-11-21 10:06:09Z satoken $

#ifndef __INC_PF_WRAPPER_H__
#define __INC_PF_WRAPPER_H__

#include <string>
#include <vector>
typedef unsigned int uint;

extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
};

class PFWrapper {
public:
  PFWrapper(const std::string& sequence,
	    bool noGU=false, bool noCloseGU=false, bool noLP=false)
    : seq_(sequence),
      sz_(seq_.size()+1), 
      pr_(sz_*(sz_+1)/2),
      iindx_(sz_),
      noGU_(noGU),
      noCloseGU_(noCloseGU),
      noLP_(noLP)
  {
  }

  ~PFWrapper()
  {
  }

  float fold(std::string& structure);

  float fold()
  {
    std::string s;
    return fold(s);
  }

  float operator()(int i, int j) const
  {
    return pr_[iindx_[i]-j];
  }

  uint size() const { return seq_.size(); }

  const std::string& sequence() const { return seq_; }

private:
  std::string seq_;
  uint sz_;
  std::vector<FLT_OR_DBL> pr_;
  std::vector<int> iindx_;
  bool noGU_;
  bool noCloseGU_;
  bool noLP_;
};

#endif	// __INC_PF_WRAPPER_H__
