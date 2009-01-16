// $Id$

#ifndef __INC_RBFKERNEL_H__
#define __INC_RBFKERNEL_H__

#include <iosfwd>
#include <vector>
#include <boost/shared_array.hpp>
#include "gradient.h"
#include "../libsvm/svm.h"

class RBFKernel
{
public:
  typedef boost::shared_array<struct svm_node> Data;
  
public:
  static void calculate_matrix(const std::vector<Data>& vec,
			       const std::vector<double>& param,
			       Kmat& kmat, Gmat& gmat);

  static void read_data(std::istream& is,
			std::vector<int>& labels,
			std::vector<Data>& vec);
};

#endif	//  __INC_RBFKERNEL_H__

// Local Variables:
// mode: C++
// End:
