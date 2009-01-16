// $Id$

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <deque>
#include <vector>
#include "../libsvm/svm.h"
#include "poly_kernel.h"
#include "gradient.h"

static inline
double
powi(double base, int times)
{
  double tmp = base, ret = 1.0;

  for(int t=times; t>0; t/=2) {
    if(t%2==1) ret*=tmp;
    tmp = tmp * tmp;
  }
  return ret;
}

// inner vector
static inline
double
dot(const struct svm_node* x, const struct svm_node* y)
{
  double sum = 0;
  while(x->index != -1 && y->index != -1) {
    if(x->index == y->index){
      sum += x->value * y->value;
      ++x;
      ++y;
    } else {
      if(x->index > y->index)
	++y;
      else
	++x;
    }			
  }
  return sum;
}

// the RBF kernel
static inline
double
kernel(const struct svm_node* x, const struct svm_node* y, int degree,
       double gamma, double coef0)
{
  return powi(gamma*dot(x,y)+coef0,degree);
}

// a gradient of the RBF kernel with respect to gamma
static inline
void
kernel_g(const struct svm_node* x, const struct svm_node* y, int degree,
	 double gamma, double coef0, double& gamma_g, double& coef0_g)
{
  double d = dot(x,y);
  double v = powi(gamma*d+coef0,degree-1)*degree;
  gamma_g = v*d;
  coef0_g = v;
}

// static
void
PolyKernel::
calculate_matrix(const std::vector<Data>& vec,
		 const std::vector<double>& param,
		 Kmat& kmat, Gmat& gmat)
{
  double k, g0, g1;
  for (uint i=0; i!=vec.size(); ++i) {
    for (uint j=i; j!=vec.size(); ++j) {
      k = kernel(&vec[i][0], &vec[j][0], degree_, param[0], param[1]);
      kernel_g(&vec[i][0], &vec[j][0], degree_, param[0], param[1], g0, g1);
      kmat[i][j] = kmat[j][i] = k;
      gmat[0][i][j] = gmat[0][j][i] = g0;
      gmat[1][i][j] = gmat[1][j][i] = g1;
    }
  }
}

// static
void
PolyKernel::
read_data(std::istream& is, std::vector<int>& lb, std::vector<Data>& vc)
{
  std::deque<int> labels;
  std::deque<Data> vec;
  std::string l;
  while (std::getline(is, l)) {
    std::istringstream ss(l);
    std::string v;
    ss >> v;
    labels.push_back(atoi(v.c_str()));
    
    uint n=0;
    while (ss >> v) n++;

    Data na=Data(new svm_node[n+1]);
    vec.push_back(na);

    ss.clear(); ss.seekg(0);
    ss >> v;
    n=0;
    while (ss >> v) {
      sscanf(v.c_str(), "%d:%lf",
	     &na[n].index, &na[n].value);
      n++;
    }
    na[n].index=-1;
  }

  lb.resize(labels.size());
  std::copy(labels.begin(), labels.end(), lb.begin());
  vc.resize(vec.size());
  std::copy(vec.begin(), vec.end(), vc.begin());
}


