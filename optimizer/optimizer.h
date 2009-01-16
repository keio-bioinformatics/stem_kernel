// $Id$

#ifndef __INC_OPTIMIZER_H__
#define __INC_OPTIMIZER_H__

#include <iosfwd>
#include <vector>

template < class CM, class GC >
class Optimizer : public CM
{
public:
  typedef typename CM::Data Data;

public:
  Optimizer() : CM() { }
  
  int optimize(const std::vector<int>& labels, const std::vector<Data>& data,
	       unsigned int ncv, std::vector<double>& param, double& C,
	       const std::vector<double>& lbd_param,
	       const std::vector<double>& ubd_param,
	       const std::vector<long>& nbd_param,
	       double eps, double factr, double pgtol, std::ostream* os);

  static void split(uint n, uint ncv, uint cv, 
		    std::vector<uint>& tr_i, std::vector<uint>& ts_i);
};

#include "optimizer.cpp"

#endif	//  __INC_OPTIMIZER_H__

// Local Variables:
// mode: C++
// End:
