// $Id$

#ifndef __INC_POLY_KERNEL_H__
#define __INC_POLY_KERNEL_H__

#include <iosfwd>
#include <vector>
#include <boost/shared_array.hpp>
#include "gradient.h"
#include "../libsvm/svm.h"

class PolyKernel
{
public:
  typedef boost::shared_array<struct svm_node> Data;

  PolyKernel() : degree_(2) { }
  
public:
  void set_degree(int degree) { degree_ = degree; }
  
  void calculate_matrix(const std::vector<Data>& vec,
			const std::vector<double>& param,
			Kmat& kmat, Gmat& gmat);

  static void read_data(std::istream& is,
			std::vector<int>& labels,
			std::vector<Data>& vec);

private:
  int degree_;
};

#endif	//  __INC_POLY_KERNEL_H__

// Local Variables:
// mode: C++
// End:
