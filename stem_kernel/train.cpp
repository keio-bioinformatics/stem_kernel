// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#ifdef HAVE_BOOST_THREAD 
#include <boost/thread.hpp>
#endif
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "stem_kernel.h"
#include "example.h"
#include "lbfgs.h"

namespace po = boost::program_options;
namespace ub = boost::numeric::ublas;

// vector and matrix types 
// with the right storage order for FORTAN (column_major):
typedef ub::vector<double, std::vector<double> > DoubleVector;
typedef ub::vector<int,    std::vector<int>    > IntVector;
typedef ub::matrix<double, ub::column_major, std::vector<double > >  DoubleMatrix;

// here we declare the external FORTRAN routine from LAPACK
extern "C" {
  void dgesv_(int* n, int* nrhs, double* a, int* lda,
	      int* ipiv, double* b, int* ldb, int* info);
  void dsyev_ (char* jobz, char* uplo, int* n, double* a, int* lda,
	       double* w, double* work, int* lwork, int* info);
};

template < class Kernel >
class Func
{
  typedef typename Kernel::value_type value_type;
  const std::vector<value_type>& x_;
  std::vector<value_type> params_;
  DoubleMatrix& k_;
  std::vector<DoubleMatrix*>& grad_;
  const ExampleSet& train_;
  const Kernel& kernel_;
  uint n_th_;
  uint th_no_;

public:
  Func(const std::vector<value_type>& x,
       DoubleMatrix& k,
       std::vector<DoubleMatrix*>& grad,
       const ExampleSet& train, const Kernel& kernel,
       uint n_th, uint th_no)
    : x_(x), params_(x_.size()), k_(k), grad_(grad), train_(train),
      kernel_(kernel), n_th_(n_th), th_no_(th_no)
  {
    std::transform(x_.begin(), x_.end(), params_.begin(), ::exp);
  }

  void operator()()
  {
    uint cnt=0;
    uint n=train_.size();
    for (uint i=0; i!=n; ++i) {
      for (uint j=i; j!=n; ++j) {
	if (cnt++%n_th_==th_no_) {
	  value_type f;
	  std::vector<value_type> gg(grad_.size());
	  std::fill(gg.begin(), gg.end(), 0);
	  f = kernel_.inside_outside(train_[i].second, train_[j].second,
				     params_, gg);
	  k_(i,j) = k_(j,i) = f;
	  for (uint l=0; l!=gg.size(); ++l) {
	    (*grad_[l])(i,j) = (*grad_[l])(j,i) = gg[l];
	  }
	}
      }
    }
  }
};


static
DoubleMatrix
log(const DoubleMatrix& m)
{
  assert(m.size1()==m.size2());
  int n=m.size1();
  DoubleMatrix v(m);
  DoubleVector e(n);
  int lwork=4*n;
  DoubleVector w(lwork);
  int info;
  dsyev_("V", "L", &n, &v(0,0), &n, &e(0), &w(0), &lwork, &info);
  DoubleMatrix diag(n,n);
  for (int i=0; i!=n; ++i) diag(i,i)=::log(e(i));
  DoubleMatrix x(ub::prod(diag, ub::trans(v)));
  return ub::prod(v, x);
}

static
double
trace(const DoubleMatrix& m)
{
  assert(m.size1()==m.size2());
  double ret=0.0;
  for (uint i=0; i!=m.size1(); ++i) ret += m(i,i);
  return ret;
}

template < class Kernel >
void
object_func(const std::vector<typename Kernel::value_type>& x,
	    typename Kernel::value_type& f, 
	    std::vector<typename Kernel::value_type>& g,
	    const ExampleSet& train, const Kernel& kernel,
	    bool normalize, uint n_th)
{
  uint n=train.size();
  DoubleMatrix k(n,n);
  std::vector<DoubleMatrix*> grad(g.size());
  for (uint i=0; i!=grad.size(); ++i)
    grad[i] = new DoubleMatrix(n,n);
  
#ifdef HAVE_BOOST_THREAD
  std::vector<boost::thread*> th(n_th);
  for (uint t=0; t!=n_th; ++t) {
    th[t] = new boost::thread(Func<Kernel>(x, k, grad, train, kernel, n_th, t));
  }
  for (uint t=0; t!=n_th; ++t) {
    th[t]->join();
    delete th[t];
  }  
#else
  Func<Kernel> func(x, k, grad, train, kernel, 1, 0);
  func();
#endif

  if (normalize) {
    DoubleMatrix k2(n,n);
    std::vector<DoubleMatrix*> grad2(g.size());
    for (uint i=0; i!=grad2.size(); ++i)
      grad2[i] = new DoubleMatrix(n,n);

    for (uint i=0; i!=n; ++i) {
      k2(i,i) = 1.0;
      for (uint l=0; l!=g.size(); ++l) (*grad2[l])(i,i) = 0.0;
      for (uint j=i+1; j!=n; ++j) {
	k2(i,j) = k2(j,i) = k(i,j)/sqrt(k(i,i)*k(j,j));
	for (uint l=0; l!=grad2.size(); ++l) {
	  (*grad2[l])(i,j) = (*grad2[l])(j,i) =
	    (*grad[l])(i,j) / sqrt(k(i,i)*k(j,j))
	    - (*grad[l])(i,i) * k2(i,j) / k(i,i) / 2
	    - (*grad[l])(j,j) * k2(i,j) / k(j,j) / 2;
	}
      }
    }

    k = k2;
    for (uint i=0; i!=grad.size(); ++i) {
      (*grad[i]) = (*grad2[i]);
      delete grad2[i];
    }
  }

  f = 0;
  std::fill(g.begin(), g.end(), 0);
  DoubleMatrix log_k(log(k));
  DoubleMatrix ent(ub::prod(k, log_k));
  f = trace(ent);
  DoubleMatrix i_k(ub::identity_matrix<double>(n)+log_k);
  for (uint i=0; i!=grad.size(); ++i) {
    DoubleMatrix dg(ub::prod(*grad[i], i_k));
    g[i] = trace(dg);
    delete grad[i];
  }
}

template < class Kernel >
void
run(const std::vector<typename Kernel::value_type>& x_init,
    const ExampleSet& train, const Kernel& kernel, bool normalize,
    std::ostream& out, uint max_it, float eps, uint n_th)
{
  typedef typename Kernel::value_type value_type;
  value_type f=0, f_prev=0;
  std::vector<value_type> x(x_init.size());
  std::vector<value_type> g(x_init.size());
  std::copy(x_init.begin(), x_init.end(), x.begin());
  std::fill(g.begin(), g.end(), 0);

  out << "0: " << f << " NA, ";
  std::copy(x.begin(), x.end(),
	    std::ostream_iterator<value_type>(out, " "));
  out << ", ";
  std::copy(g.begin(), g.end(),
	    std::ostream_iterator<value_type>(out, " "));
  out << std::endl;

  LBFGS lbfgs;
  lbfgs.initialize(x.size(), 5);
  
  {
    object_func(x, f, g, train, kernel, normalize, n_th);
    lbfgs.update(&x[0], &f, &g[0]);
  }
  f_prev=f;

  out << "1: " << f << " NA, ";
  std::copy(x.begin(), x.end(),
	    std::ostream_iterator<value_type>(out, " "));
  out << ", ";
  std::copy(g.begin(), g.end(),
	    std::ostream_iterator<value_type>(out, " "));
  out << std::endl;

  for (uint i=1; i<max_it; ++i) {
    object_func(x, f, g, train, kernel, normalize, n_th);
    lbfgs.update(&x[0], &f, &g[0]);

    double conv_rate=std::abs((f-f_prev)/f);
    out << i+1 << ": " << f << " " << conv_rate << ", ";
    std::copy(x.begin(), x.end(),
	      std::ostream_iterator<value_type>(out, " "));
  out << ", ";
    std::copy(g.begin(), g.end(),
	      std::ostream_iterator<value_type>(out, " "));
    out << std::endl;

    if (conv_rate<eps) break;
    f_prev=f;
  }
}

int
main(int argc, char* argv[])
{
  typedef double value_type;
  float gap;
  float stack;
  float subst;
  uint loop;
  uint n_th;
  uint max_it;
  float eps;
  std::string training_file;
  std::string output_file;
  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("gap,g", po::value<float>(&gap)->default_value(0.8),
     "set the gap weight")
    ("stack,s", po::value<float>(&stack)->default_value(1.0),
     "set the stacking weight")
    ("loop,l", po::value<uint>(&loop)->default_value(0),
     "set minimum loop length")
    ("normalize,n", "normalize the kernel matrix")
    ("enable-wobble-pair,w", "allow wobble base pairs")
    ("substitution,v", po::value<float>(&subst)->default_value(0.5),
     "set substitution weight for base pairs")
    ("threads,t", po::value<uint>(&n_th)->default_value(1),
     "set the number of threads")
    ("max-iter,i", po::value<uint>(&max_it)->default_value(50),
     "set the maximum number of iteration for training process")
    ("conv,e", po::value<float>(&eps)->default_value(1e-5),
     "set the convergence rate to stop iteration for training process");
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("output-file", po::value<std::string>(&output_file),
     "set the output filename")
    ("training-file", po::value<std::string>(&training_file),
     "set the filename which includes training data");
  po::options_description cmd_opts("Options");
  cmd_opts.add(desc).add(hidden);
  po::variables_map vm;
  po::positional_options_description p;
  p.add("output-file", 1);
  p.add("training-file", 1);
  po::store(po::command_line_parser(argc, argv).
	    options(cmd_opts).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help") ||
      !vm.count("training-file") || !vm.count("output-file")) {
    std::cout << "Kernel Gram Matrix Calculator for Stem Kernel" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] output-file training-data\n\n"
	      << desc << std::endl;
    return 1;
  }

  // load examples
  ExampleSet train;
  uint n=0, m=0;
  {
    std::ifstream in(training_file.c_str());
    if (!in.is_open()) {
      perror(training_file.c_str());
      return 1;
    }
    load_examples(in, train);
    m = train.size();
    n = m;
  }

  std::vector<value_type> x_init(3);
  x_init[0]=::log(gap);
  x_init[1]=::log(stack);
  x_init[2]=::log(subst);
  std::ofstream out(output_file.c_str());
  if (!vm.count("enable-wobble-pair")) {
    StemKernel<value_type,NormalBasePair> kernel(loop, gap, stack, subst);
    run(x_init, train, kernel, vm.count("normalize"), out, max_it, eps, n_th);
  } else {
    StemKernel<value_type,WobbleBasePair> kernel(loop, gap, stack, subst);
    run(x_init, train, kernel, vm.count("normalize"), out, max_it, eps, n_th);
  }

  return 0;
}
