// $Id$

#ifndef __TMPL_OPTIMIZER_CPP__
#define __TMPL_OPTIMIZER_CPP__

#include <iostream>
#include <boost/multi_array.hpp>
#include "lbfgsb.h"
#include "optimizer.h"

template < class CM, class GC >
int
Optimizer<CM,GC>::
optimize(const std::vector<int>& labels, const std::vector<Data>& data,
	 uint ncv, std::vector<double>& param, double& C,
	 const std::vector<double>& lbd_param,
	 const std::vector<double>& ubd_param,
	 const std::vector<long>& nbd_param,
	 double eps, double factr, double pgtol, std::ostream* os)
{
  assert(param.size()==lbd_param.size());
  assert(param.size()==ubd_param.size());
  assert(param.size()==nbd_param.size());

  LBFGSB lbfgsb(factr, pgtol /*1e10, 1e-5*/);
  std::vector<double> lbd(param.size()+1);
  std::vector<double> ubd(param.size()+1);
  std::vector<long> nbd(param.size()+1);
  lbd[0]=1e-5; nbd[0]=LBFGSB::LOWER_BOUND;
  for (uint i=1; i!=nbd.size(); ++i) {
    lbd[i] = lbd_param[i-1];
    ubd[i] = ubd_param[i-1];
    nbd[i] = nbd_param[i-1];
  }
  lbfgsb.initialize(param.size()+1, 5, &lbd[0], &ubd[0], &nbd[0]);

  boost::multi_array<double,2> kmat(boost::extents[data.size()][data.size()]);
  boost::multi_array<double,3> gmat(boost::extents[param.size()][data.size()][data.size()]);

  uint count=0;
  int iflag=1;
  do {
    if (os) *os << "=== Step " << ++count << " ===" << std::endl;
    
    // calculate the kernel matrix and its gradients
    if (os) *os << "- calculate a kernel matrix" << std::flush;
    CM::calculate_matrix(data, param, kmat, gmat);
    if (os) *os << " done." << std::endl;

    if (os) *os << "- gradient calculation " << std::endl;
    std::vector<double> x(param.size()+1);
    x[0] = C; 
    for (uint i=1; i!=x.size(); ++i) x[i] = param[i-1];
    if (os) {
      *os << " x=[ ";
      std::copy(x.begin(), x.end(), std::ostream_iterator<double>(*os, " "));
      *os << "]" << std::endl << "  " << std::flush;
    }

    double f=0.0;
    std::vector<double> g(param.size()+1, 0.0);
    for (uint i=0; i!=ncv; ++i) {
      // make cross-validation data set
      std::vector<uint> tr_i, ts_i;
      split(labels.size(), ncv, i, tr_i, ts_i);

      // compute gradients
      std::vector<double> fg(param.size());
      double cg;
      GC gc(labels, tr_i, ts_i);
      double f0 = gc.compute(kmat, gmat, C, fg, cg, eps);
      f -= f0;
      g[0] -= cg;
      for (uint j=1; j!=g.size(); ++j)
	g[j] -= fg[j-1];
      if (os) *os << '.' << std::flush;
    }
    if (os) {
      *os << " done." << std::endl
	  << " f= " << -f << std::endl
	  << " g=[ ";
      std::transform(g.begin(), g.end(),
		     std::ostream_iterator<double>(*os, " "),
		     std::negate<double>());
      *os << "]" << std::endl << std::endl;
    }
    
    // update parameters
    iflag=lbfgsb.update(&x[0], f, &g[0]);
    C = x[0];
    for (uint i=1; i!=x.size(); ++i) param[i-1] = x[i];

  } while (iflag>0);

  return iflag;
}

//static
template < class CM, class GC >
void
Optimizer<CM,GC>::
split(uint n, uint ncv, uint cv, std::vector<uint>& tr_i, std::vector<uint>& ts_i)
{
  uint n_tr=0;
  for (uint j=0; j!=n; ++j) {
    if (j % ncv != cv) n_tr++;
  }
  tr_i.resize(n_tr);
  ts_i.resize(n-n_tr);
  for (uint j=0, k=0, l=0; j!=n; ++j) {
    if (j % ncv != cv)
      tr_i[k++] = j;
    else
      ts_i[l++] = j;
  }
}

#endif	//  __TMPL_OPTIMIZER_CPP__
