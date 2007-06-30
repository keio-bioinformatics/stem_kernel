// $Id$

#ifndef __TMPL_KERNEL_MATRIX_CPP__
#define __TMPL_KERNEL_MATRIX_CPP__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cmath>
#include "kernel_matrix.h"
#ifdef HAVE_MPI
#include <mpi2c++/mpi++.h>
#else
#ifdef HAVE_BOOST_THREAD 
#include <boost/thread.hpp>
#endif
#endif

#ifndef HAVE_MPI
template < class Kernel, class Matrix, class ExampleSet >
class CalcTrainMatrix
{
  typedef typename Kernel::value_type value_type;
  uint n_th_;
  uint th_no_;
  const ExampleSet& train_;
  const Kernel& kernel_;
  Matrix& matrix_;

public:
  CalcTrainMatrix(uint n_th, uint th_no, const ExampleSet& train,
		  const Kernel& kernel, Matrix& matrix)
    : n_th_(n_th), th_no_(th_no), train_(train),
      kernel_(kernel), matrix_(matrix)
  {
  }

  void operator()()
  {
    uint cnt=0;
    uint n = train_.size();
    for (uint i=0; i!=n; ++i) {
      for (uint j=i; j!=n; ++j) {
	if (cnt++%n_th_==th_no_) {
	  matrix_(i,j)=kernel_(train_[i].second, train_[j].second);
	  if (i!=j) matrix_(j,i)=matrix_(i,j);
	}
      }
    }
  }
};

template < class Kernel, class ExampleSet >
class CalcDiagonal
{
  typedef typename Kernel::value_type value_type;
  uint n_th_;
  uint th_no_;
  const ExampleSet& train_;
  const std::vector<uint>& sv_index_;
  const Kernel& kernel_;
  std::vector<value_type>& diag_;

public:
  CalcDiagonal(uint n_th, uint th_no, const ExampleSet& train,
	       const std::vector<uint>& sv_index,
	       const Kernel& kernel, std::vector<value_type>& diag)
    : n_th_(n_th), th_no_(th_no), train_(train),
      sv_index_(sv_index), kernel_(kernel), diag_(diag)
  {
  }

  CalcDiagonal(uint n_th, uint th_no, const ExampleSet& train,
	       const Kernel& kernel, std::vector<value_type>& diag)
    : n_th_(n_th), th_no_(th_no), train_(train),
      sv_index_(), kernel_(kernel), diag_(diag)
  {
  }

  void operator()()
  {
    if (sv_index_.empty()) {
      for (uint i=0; i!=train_.size(); ++i) {
	if (i%n_th_==th_no_) {
	  diag_[i]=kernel_(train_[i].second, train_[i].second);
	}
      }
    } else {
      for (uint i=0; i!=sv_index_.size(); ++i) {
	if (i%n_th_==th_no_) {
	  uint x = sv_index_[i];
	  diag_[x] = kernel_(train_[x].second, train_[x].second);
	}
      }
    }
  }
};

template < class Kernel, class ExampleSet >
class CalcTestMatrix
{
  typedef typename Kernel::value_type value_type;
  typedef typename ExampleSet::value_type Example;

  uint n_th_;
  uint th_no_;
  const Example& data_;
  const ExampleSet& train_;
  const std::vector<uint>& sv_index_;
  const Kernel& kernel_;
  std::vector<value_type>& matrix_;
  value_type* data_self_;

public:
  CalcTestMatrix(uint n_th, uint th_no,
		 const Example& data, const ExampleSet& train,
		 const std::vector<uint>& sv_index,
		 const Kernel& kernel, std::vector<value_type>& matrix,
		 value_type* data_self=NULL)
    : n_th_(n_th), th_no_(th_no), data_(data), train_(train),
      sv_index_(sv_index), kernel_(kernel), matrix_(matrix),
      data_self_(data_self)
  {
  }

  CalcTestMatrix(uint n_th, uint th_no,
		 const Example& data, const ExampleSet& train,
		 const Kernel& kernel, std::vector<value_type>& matrix,
		 value_type* data_self=NULL)
    : n_th_(n_th), th_no_(th_no), data_(data), train_(train),
      sv_index_(), kernel_(kernel), matrix_(matrix),
      data_self_(data_self)
  {
  }

  void operator()()
  {
    uint cnt=0;
    if (sv_index_.empty()) {
      for (uint i=0; i!=train_.size(); ++i) {
	if (cnt++%n_th_==th_no_) {
	  matrix_[i] = kernel_(train_[i].second, data_.second);
	}
      }
    } else {
      for (uint i=0; i!=sv_index_.size(); ++i) {
	if (cnt++%n_th_==th_no_) {
	  uint x = sv_index_[i];
	  matrix_[x] = kernel_(train_[x].second, data_.second);
	}
      }
    }
    if (data_self_) {
      if (cnt++%n_th_==th_no_)
	*data_self_=kernel_(data_.second, data_.second);
      cnt++;
    }
  }
};

#else

template < class Kernel, class ExampleSet >
class CalcTrainMatrix
{
  typedef typename Kernel::value_type value_type;
  uint n_th_;
  uint th_no_;
  const ExampleSet& train_;
  const Kernel& kernel_;
  uint num_;
  std::vector<double> vals_;

public:
  CalcTrainMatrix(uint n_th, uint th_no, const ExampleSet& train,
		  const Kernel& kernel)
    : n_th_(n_th), th_no_(th_no), train_(train),
      kernel_(kernel), num_(0), vals_()
  {
    uint n=train.size();
    vals_.resize((n+1)*n/2/n_th_+1);
  }

  void operator()()
  {
    uint cnt=0;
    uint n = train_.size();
    for (uint i=0; i!=n; ++i) {
      for (uint j=i; j!=n; ++j) {
	if (cnt++%n_th_==th_no_) {
#if 0
	  std::cout << "calculating matrix["
		    << i << "][" << j << "]" << std::endl;
#endif
	  vals_[num_++]=kernel_(train_[i].second, train_[j].second);
	}
      }
    }
  }

  void send(uint dest, int tag)
  {
    MPI::COMM_WORLD.Ssend(&num_, 1, MPI::UNSIGNED, dest, tag);
    MPI::COMM_WORLD.Ssend(&vals_[0], num_, MPI::DOUBLE, dest, tag);
  }

  template < class Matrix >
  void recv_helper(uint src, const std::vector<double>& v, Matrix& matrix)
  {
    uint cnt=0;
    uint x=0;
    uint n = train_.size();
    for (uint i=0; i!=n; ++i) {
      for (uint j=i; j!=n; ++j) {
	if (cnt++%n_th_==src) {
	  matrix(i,j) = v[x++];
	  if (i!=j) matrix(j,i)=matrix(i,j);
	}
      }
    }
  }
  

  template < class Matrix >
  void recv(uint src, int tag, Matrix& matrix)
  {
    // receive
    if (src!=0) {
      uint num;
      MPI::COMM_WORLD.Recv(&num, 1, MPI::UNSIGNED, src, tag);
      std::vector<double> v(num);
      MPI::COMM_WORLD.Recv(&v[0], num, MPI::DOUBLE, src, tag);
      recv_helper(src, v, matrix);
    } else {
      recv_helper(src, vals_, matrix);
    }
  }  
};

template < class Kernel, class ExampleSet >
class CalcDiagonal
{
  typedef typename Kernel::value_type value_type;
  uint n_th_;
  uint th_no_;
  const ExampleSet& train_;
  const std::vector<uint> sv_index_;
  const Kernel& kernel_;
  uint num_;
  std::vector<double> vals_;

public:
  CalcDiagonal(uint n_th, uint th_no, const ExampleSet& train,
	       const std::vector<uint>& sv_index,
	       const Kernel& kernel)
    : n_th_(n_th), th_no_(th_no), train_(train), sv_index_(sv_index),
      kernel_(kernel), num_(0), vals_()
  {
    uint n=sv_index_.size();
    if (n==0) n=train_.size();
    vals_.resize(n/n_th_+1);
  }

  CalcDiagonal(uint n_th, uint th_no, const ExampleSet& train,
	       const Kernel& kernel)
    : n_th_(n_th), th_no_(th_no), train_(train), sv_index_(),
      kernel_(kernel), num_(0), vals_()
  {
    uint n=train_.size();
    vals_.resize(n/n_th_+1);
  }

  void operator()()
  {
    if (sv_index_.empty()) {
      for (uint i=0; i!=train_.size(); ++i) {
	if (i%n_th_==th_no_) {
	  vals_[num_++]=kernel_(train_[i].second, train_[i].second);
	}
      }
    } else {
      for (uint i=0; i!=sv_index_.size(); ++i) {
	if (i%n_th_==th_no_) {
	  uint x = sv_index_[i];
	  vals_[num_++]=kernel_(train_[x].second, train_[x].second);
	}
      }
    }
  }

  void send(uint dest, int tag)
  {
    MPI::COMM_WORLD.Ssend(&num_, 1, MPI::UNSIGNED, dest, tag);
    MPI::COMM_WORLD.Ssend(&vals_[0], num_, MPI::DOUBLE, dest, tag);
  }

  void recv_helper(uint src, std::vector<value_type>& v,
		   std::vector<value_type>& diag)
  {
    uint x=0;
    if (sv_index_.empty()) {
      for (uint i=0; i!=train_.size(); ++i) {
	if (i%n_th_==src) {
	  diag[i] = v[x++];
	}
      }
    } else {
      for (uint i=0; i!=sv_index_.size(); ++i) {
	if (i%n_th_==src) {
	  diag[sv_index_[i]] = v[x++];
	}
      }
    }
  }

  void recv(uint src, int tag, std::vector<value_type>& diag)
  {
    // receive
    if (src!=0) {
      uint num;
      MPI::COMM_WORLD.Recv(&num, 1, MPI::UNSIGNED, src, tag);
      std::vector<double> v(num);
      MPI::COMM_WORLD.Recv(&v[0], num, MPI::DOUBLE, src, tag);
      recv_helper(src, v, diag);
    } else {
      recv_helper(src, vals_, diag);
    }
  }  
};

template < class Kernel, class ExampleSet >
class CalcTestMatrix
{
  typedef typename Kernel::value_type value_type;
  typedef typename ExampleSet::value_type Example;

  uint n_th_;
  uint th_no_;
  const Example& data_;
  const ExampleSet& train_;
  const std::vector<uint>& sv_index_;
  const Kernel& kernel_;
  bool norm_test_;
  uint num_;
  std::vector<double> vals_;

public:
  CalcTestMatrix(uint n_th, uint th_no,
		 const Example& data, const ExampleSet& train,
		 const std::vector<uint>& sv_index,
		 const Kernel& kernel, bool norm_test=false)
    : n_th_(n_th), th_no_(th_no), data_(data), train_(train),
      sv_index_(sv_index), kernel_(kernel), norm_test_(norm_test),
      num_(0), vals_()
  {
    uint n = sv_index_.size();
    if (n==0) n=train_.size();
    if (norm_test_) n += 1;
    vals_.resize(n/n_th_+1);
  }

  CalcTestMatrix(uint n_th, uint th_no,
		 const Example& data, const ExampleSet& train,
		 const Kernel& kernel, bool norm_test=false)
    : n_th_(n_th), th_no_(th_no), data_(data), train_(train),
      sv_index_(), kernel_(kernel), norm_test_(norm_test),
      num_(0), vals_()
  {
    uint n = train_.size();
    if (norm_test_) n += 1;
    vals_.resize(n/n_th_+1);
  }

  void operator()()
  {
    uint cnt=0;
    if (sv_index_.empty()) {
      for (uint i=0; i!=train_.size(); ++i) {
	if (cnt++%n_th_==th_no_) {
	  vals_[num_++]=kernel_(train_[i].second, data_.second);
	}
      }
    } else {
      for (uint i=0; i!=sv_index_.size(); ++i) {
	if (cnt++%n_th_==th_no_) {
	  vals_[num_++]=kernel_(train_[sv_index_[i]].second, data_.second);
	}
      }
    }
    if (norm_test_) {
      if (cnt++%n_th_==th_no_) {
	vals_[num_++]=kernel_(data_.second, data_.second);
      }
    }
  }

  void send(uint dest, int tag)
  {
    MPI::COMM_WORLD.Ssend(&num_, 1, MPI::UNSIGNED, dest, tag);
    MPI::COMM_WORLD.Ssend(&vals_[0], num_, MPI::DOUBLE, dest, tag);
  }

  void recv_helper(uint src, const std::vector<double>& v,
		   std::vector<value_type>& matrix, value_type& data_self)
  {
    uint x=0;
    uint cnt=0;
    if (sv_index_.empty()) {
      for (uint i=0; i!=train_.size(); ++i) {
	if (cnt++%n_th_==src) matrix[i]=v[x++];
      }
    } else {
      for (uint i=0; i!=sv_index_.size(); ++i) {
	if (cnt++%n_th_==src) matrix[sv_index_[i]]=v[x++];
      }
    }
    if (norm_test_) {
      if (cnt++%n_th_==src) data_self=v[x++];
    }
  }
  
  void recv(uint src, int tag,
	    std::vector<value_type>& matrix, value_type& data_self)
  {
    if (src!=0) {
      uint num;
      MPI::COMM_WORLD.Recv(&num, 1, MPI::UNSIGNED, src, tag);
      std::vector<double> v(num);
      MPI::COMM_WORLD.Recv(&v[0], num, MPI::DOUBLE, src, tag);
      recv_helper(src, v, matrix, data_self);
    } else {
      recv_helper(src, vals_, matrix, data_self);
    }
  }

  void recv(uint src, int tag, std::vector<value_type>& matrix)
  {
    value_type x;
    recv(src, tag, matrix, x);
  }
  
};
#endif

template < class ValueType>
template < class Kernel, class ExampleSet >
void
KernelMatrix<ValueType>::
calculate(const ExampleSet& train, 
	  const Kernel& kernel, bool normalize, uint n_th)
{
  uint n=train.size();
  //resize(n,n);
#ifdef HAVE_MPI
  typedef CalcTrainMatrix<Kernel, ExampleSet > CalcMatrix;
  uint rank = MPI::COMM_WORLD.Get_rank();
  uint num_procs = MPI::COMM_WORLD.Get_size();
  CalcMatrix calc(num_procs, rank, train, kernel);
  calc();
  if (rank==0) {
    resize(n, n);
    for (uint i=0; i!=n; ++i) label_[i]=train[i].first;
    //label_.resize(n);
    //matrix_.resize(boost::extents[n][n]);
    for (uint t=0; t!=num_procs; ++t) {
      calc.recv(t, 0, *this);
    }
    if (normalize) {
      for (uint i=0; i!=n-1; ++i) {
	for (uint j=i+1; j!=n; ++j) {
	  matrix_[i][j]/=sqrt(matrix_[i][i]*matrix_[j][j]);
	  matrix_[j][i]=matrix_[i][j];
	}
      }
      for (uint i=0; i!=n; ++i)
	matrix_[i][i]=1;
    }
  } else {
    calc.send(0, 0);
  }
#else  // HAVE_MPI

  typedef CalcTrainMatrix<Kernel, KernelMatrix< value_type >, ExampleSet > CalcMatrix;
  resize(n);
  //matrix_.resize(boost::extents[n][n]);
  //label_.resize(n);
  for (uint i=0; i!=n; ++i) label_[i]=train[i].first;

#ifdef HAVE_BOOST_THREAD
  std::vector<boost::thread*> th(n_th);
  for (uint t=0; t!=n_th; ++t) {
    th[t]=new boost::thread(CalcMatrix(n_th, t, train, kernel, *this));
  }
  for (uint t=0; t!=n_th; ++t) {
    th[t]->join();
    delete th[t];
  }
#else  // HAVE_BOOST_THREAD
  CalcMatrix calc(1, 0, train, kernel, *this);
  calc();
#endif	// HAVE_BOOST_THREAD
  if (normalize) {
    for (uint i=0; i!=n-1; ++i) {
      for (uint j=i+1; j!=n; ++j) {
	matrix_[i][j]/=sqrt(matrix_[i][i]*matrix_[j][j]);
	matrix_[j][i]=matrix_[i][j];
      }
    }
    for (uint i=0; i!=n; ++i)
      matrix_[i][i]=1;
  }
#endif  // HAVE_MPI
}

//static
template < class ValueType>
template < class Kernel, class ExampleSet >
void
KernelMatrix<ValueType>::
diagonal(std::vector<value_type>& diag, const ExampleSet& train,
	 const std::vector<uint>& sv_index, const Kernel& kernel, uint n_th)
{
#ifdef HAVE_MPI
  typedef CalcDiagonal< Kernel, ExampleSet > CalcDiag;
  uint rank = MPI::COMM_WORLD.Get_rank();
  uint num_procs = MPI::COMM_WORLD.Get_size();
  CalcDiag calc_diag(num_procs, rank, train, sv_index, kernel);
  calc_diag();
  if (rank==0) {
    //diag.resize(train.size());
    for (uint t=0; t!=num_procs; ++t) {
      calc_diag.recv(t, 0, diag);
    }
  } else { 
    calc_diag.send(0, 0);
  }
#else
  typedef CalcDiagonal< Kernel, ExampleSet > CalcDiag;
  //diag.resize(train.size());
#ifdef HAVE_BOOST_THREAD
  std::vector<boost::thread*> th(n_th);
  for (uint t=0; t!=n_th; ++t) {
    th[t] = new boost::thread(CalcDiag(n_th, t, train, sv_index, kernel, diag));
  }
  for (uint t=0; t!=n_th; ++t) {
    th[t]->join();
    delete th[t];
  }
#else  // HAVE_BOOST_THREAD
  CalcDiag calc_diag(1, 0, train, sv_index, kernel, diag);
  calc_diag();
#endif
#endif  
}

template < class ValueType>
template < class Kernel, class ExampleSet >
void
KernelMatrix<ValueType>::
calculate(std::vector<value_type>& matrix,
	  const typename ExampleSet::value_type& data, const ExampleSet& train,
	  const std::vector<uint>& sv_index,
	  const Kernel& kernel, uint n_th, value_type* data_self)
{
#ifdef HAVE_MPI
  typedef CalcTestMatrix< Kernel, ExampleSet > CalcMatrix;
  uint rank = MPI::COMM_WORLD.Get_rank();
  uint num_procs = MPI::COMM_WORLD.Get_size();
  CalcMatrix calc(num_procs, rank, data, train, sv_index,
		  kernel, data_self!=NULL);
  calc();
  if (rank==0) {
    //matrix.resize(train.size());
    for (uint t=0; t!=num_procs; ++t) {
      if (data_self)
	calc.recv(t, 1, matrix, *data_self);
      else
	calc.recv(t, 1, matrix);
    }
  } else {
    calc.send(0, 1);
  }
#else
  typedef CalcTestMatrix< Kernel, ExampleSet > CalcMatrix;
  //matrix.resize(train.size());
#ifdef HAVE_BOOST_THREAD
  std::vector<boost::thread*> th(n_th);
  for (uint t=0; t!=n_th; ++t) {
    th[t] = new boost::thread(CalcMatrix(n_th, t, data, train, sv_index,
					 kernel, matrix, data_self));
  }
  for (uint t=0; t!=n_th; ++t) {
    th[t]->join();
    delete th[t];
  }
#else
  CalcMatrix calc(1, 0, data, train, sv_index, kernel, matrix, data_self);
  calc();
#endif
#endif
}

template < class ValueType>
template < class Kernel, class ExampleSet >
void
KernelMatrix<ValueType>::
calculate(const ExampleSet& test, const ExampleSet& train,
	  const Kernel& kernel, bool norm_test, bool normalize, uint n_th)
{
#ifdef HAVE_MPI
  if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
    resize(test.size(), train.size());
    //label_.resize(test.size());
    //matrix_.resize(boost::extents[test.size()][train.size()]);
    if (norm_test||normalize) self_.resize(train.size());
#ifdef HAVE_MPI
  }
#endif

  for (uint i=0; i!=test.size(); ++i) {
    std::vector<value_type> v(train.size());
    if (norm_test||normalize)
      calculate(v, test[i], train, kernel, n_th, &self_[i]);
    else
      calculate(v, test[i], train, kernel, n_th);

#ifdef HAVE_MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
    label_[i] = test[i].first;
    for (uint j=0; j!=train.size(); ++j)
      matrix_[i][j] = v[j];
#ifdef HAVE_MPI
    }
#endif
  }

  if (normalize) {
    std::vector<value_type> diag(train.size());
    diagonal(diag, train, kernel, n_th);
#ifdef HAVE_MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
      for (uint i=0; i!=test.size(); ++i) {
	for (uint j=0; j!=train.size(); ++j) {
	  matrix_[i][j] /= sqrt(self_[i]*diag[j]);
	}
      }
#ifdef HAVE_MPI
    }
#endif
  }
}

template < class ValueType>
void
KernelMatrix<ValueType>::
print(std::ostream& out) const
{
  uint n = matrix_.size();
  for (uint i=0; i!=n; ++i) {
    out << label_[i] << " 0:" << (i+1) << " ";
    uint m = matrix_[i].size();
    for (uint j=0; j!=m; ++j) {
      out << (j+1) << ":" << matrix_[i][j] << " ";
    }
    out << std::endl;
  }
}

#endif
