// $Id$

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <deque>
#include <vector>
#include "../libsvm/svm.h"
#include "rbf_kernel.h"
#include "gradient.h"

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
kernel(const struct svm_node* x, const struct svm_node* y, double gamma)
{
  return exp(-gamma * (dot(x,x)+dot(y,y)-2*dot(x,y)));
}

// a gradient of the RBF kernel with respect to gamma
static inline
double
kernel_g(const struct svm_node* x, const struct svm_node* y, double gamma)
{
  double w = dot(x,x)+dot(y,y)-2*dot(x,y);
  return -w * exp(-gamma * w);
}

// static
void
RBFKernel::
calculate_matrix(const std::vector<Data>& vec,
		 const std::vector<double>& param,
		 Kmat& kmat, Gmat& gmat)
{
  double k, g;
  for (uint i=0; i!=vec.size(); ++i) {
    for (uint j=i; j!=vec.size(); ++j) {
      k = kernel(&vec[i][0], &vec[j][0], param[0]);
      g = kernel_g(&vec[i][0], &vec[j][0], param[0]);
      kmat[i][j] = kmat[j][i] = k;
      gmat[0][i][j] = gmat[0][j][i] = g;
    }
  }
}

// static
void
RBFKernel::
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


