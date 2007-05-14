// $Id$

#ifndef __INC_KERNEL_MATRIX_H__
#define __INC_KERNEL_MATRIX_H__

#include <iostream>
#include <vector>
#include <utility>
#include <boost/multi_array.hpp>

template <class ValueType>
class KernelMatrix {
public:
  typedef ValueType value_type;

private:
  uint row_;
  uint col_;
  boost::multi_array<value_type,2> matrix_;
  std::vector<value_type> self_;
  std::vector<std::string> label_;

public:
  KernelMatrix() : matrix_(), self_(), label_()
  {
  }
  
  KernelMatrix(uint row, uint col) :
    row_(row), col_(col),
    matrix_(boost::extents[row][col]), self_(row), label_(row)
  {
  }

  value_type& operator()(uint x, uint y)
  {
    return matrix_[x][y];
  }

  const value_type& operator()(uint x, uint y) const
  {
    return matrix_[x][y];
  }

  value_type& operator()(uint x)
  {
    return self_[x];
  }

  const value_type& operator()(uint x) const
  {
    return self_[x];
  }

  const std::vector<value_type>& self() const
  {
    return self_;
  }

  template < class Kernel, class ExampleSet >
  void calculate(const ExampleSet& train, const Kernel& kernel,
		 bool normalize=false, uint n_th=1);

  template < class Kernel, class ExampleSet >
  void calculate(const ExampleSet& test, const ExampleSet& train,
		 const Kernel& kernel,
		 bool norm_test=false, bool normalize=false, uint n_th=1);

  template < class Kernel, class ExampleSet >
  static void calculate(std::vector<value_type>& matrix,
			const typename ExampleSet::value_type& data,
			const ExampleSet& train,
			const std::vector<uint>& sv_index, const Kernel& kernel,
			uint n_th=1, value_type* data_self=NULL);

  template < class Kernel, class ExampleSet >
  static void calculate(std::vector<value_type>& matrix,
			const typename ExampleSet::value_type& data,
			const ExampleSet& train, const Kernel& kernel,
			uint n_th=1, value_type* data_self=NULL)
  {
    std::vector<uint> idx;
    calculate(matrix, data, train, idx, kernel, n_th, data_self);
  }
  
  template < class Kernel, class ExampleSet >
  static void diagonal(std::vector<value_type>& diag, const ExampleSet& train,
		       const std::vector<uint>& sv_index,
		       const Kernel& kernel, uint n_th=1);

  template < class Kernel, class ExampleSet >
  static void diagonal(std::vector<value_type>& diag, const ExampleSet& train,
			const Kernel& kernel, uint n_th=1)
  {
   std::vector<uint> idx;
    diagonal(diag, train, idx, kernel, n_th);
  }

  void print(std::ostream& out) const;
};

#include "kernel_matrix.cpp"

#endif // __INC_MATRIX_H__

// Local Variables:
// mode: C++
// End:
