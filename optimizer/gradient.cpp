// $Id$

#include <vector>
#include <list>
#include <boost/multi_array.hpp>
#include <boost/range.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "../libsvm/solver.h"
#include "../libsvm/qmatrix.h"
#include "gradient.h"

namespace bnu = boost::numeric::ublas;

#define SIGMOID_CONST 10.0

static
void
svm_train(const svm_problem *prob, double C, double eps, double *alpha, double& b);

static
double
svm_predict(const svm_problem *prob, const svm_node *x,
	    const double* alpha, double b);

template < class K >
static
struct svm_node*
make_node(const K& k, uint i);

static
void
destroy_node(struct svm_node* x);

template < class K >
static
struct svm_problem*
make_problem(const std::vector<uint>& tr_i, const std::vector<int>& y, const K& k);

static
void
destroy_problem(struct svm_problem* prob);

static inline
double
sigmoid(double x, double s);

template < class V, class W >
static
void
calculate_avg_var(const V& v, W& avg, W& var);

static
void
find_support_vectors(const std::vector<uint>& tr_i,
		     double C,
		     const std::vector<double>& alpha,
		     bnu::vector<double>& alpha_u,
		     std::vector<uint>& u_idx,
		     std::vector<uint>& c_idx,
		     std::vector<uint>& z_idx);

template < class K >
static
void
solve_d(const K& kmat, const std::vector<int>& yvec,
	const std::vector<uint>& ts_i,
	const std::vector<uint>& u_idx, 
	const std::vector<double>& delta,
	bnu::vector<double>& d_u);

template < class K >
double
calculate_gradient_c(const K& k,
		     const std::vector<int>& y,
		     const std::vector<uint>& u_idx,
		     const std::vector<uint>& c_idx,
		     const std::vector<uint>& ts_i,
		     const bnu::vector<double>& d_u,
		     const std::vector<double>& delta);

template < class K, class G >
double
calculate_gradient_p(const K& k,
		     const G& g,
		     const std::vector<int>& y,
		     uint p, double C,
		     const std::vector<uint>& u_idx,
		     const std::vector<uint>& c_idx,
		     const std::vector<uint>& tr_i,
		     const std::vector<uint>& ts_i,
		     const bnu::vector<double>& d_u,
		     const std::vector<double>& delta,
		     const std::vector<double>& alpha,
		     const bnu::vector<double>& alpha_u,
		     double b);

template < class T >
void
conjugate_gradient(const bnu::matrix<T>& A, const bnu::vector<T>& b,
		   bnu::vector<T>& x, double tol=1e-10);

double
GradientComputation::
compute(const Kmat& k, const Gmat& g,
	double C, std::vector<double>& fg, double& cg, double eps/*=1e-3*/) const
{
  assert(y_.size()==k.size());
  assert(y_.size()==g[0].size());

  // 1. execute svm_train for tr using current paramters,
  //    and calculate coefficients alpha
  struct svm_problem* prob = make_problem(tr_i_, y_, k);
  std::vector<double> alpha(prob->l, 0.0);
  double b=0.0;
  svm_train(prob, C, eps, &alpha[0], b);
  
  // 2. execute svm_predict for ts using current parameters, alpha and C,
  //    and calculate decision values
  std::vector<double> dec_values(ts_i_.size());
  for (uint i=0; i!=ts_i_.size(); ++i) {
    struct svm_node* x = make_node(k, ts_i_[i]);
    dec_values[i] = svm_predict(prob, x, &alpha[0], b);
    destroy_node(x);
  }

  // 3. calculate f and its gradients delta with respect to decision values
  double f;
  std::vector<double> delta(dec_values.size());
  f = compute_delta(dec_values, delta);

  // 4. make the d from the kernel matrix and delta
  std::vector<uint> u_idx, c_idx, z_idx;
  bnu::vector<double> alpha_u;
  find_support_vectors(tr_i_, C, alpha, alpha_u, u_idx, c_idx, z_idx);
  //std::cerr << "(" << u_idx.size() << "," << c_idx.size() << "," << z_idx.size() << ")";
  bnu::vector<double> d_u;
  solve_d(k, y_, ts_i_, u_idx, delta, d_u);

  // 5. for each hyper parameter,
  //    calculate gradients of f from d, the kernel matrix, delta, and so forth.
  // 5.1. calculate a gradient of f with respect to C
  cg = calculate_gradient_c(k, y_, u_idx, c_idx, ts_i_, d_u, delta);
  // 5.2. calculate a gradient with respect to each hyperparameter
  for (uint p=0; p!=fg.size(); ++p) {
    fg[p] = calculate_gradient_p(k, g, y_, p, C, u_idx, c_idx, tr_i_,
				 ts_i_, d_u, delta, alpha, alpha_u, b);
  }

  destroy_problem(prob);

  return f;
}


double
GradientComputationAUC::
compute_delta(const std::vector<double>& dec_values,
	      std::vector<double>& delta) const
{
  assert(ts_i_.size()==dec_values.size());
  //assert(ts_i_.size()==delta.size());

  std::list<double> d;
  for (uint i=0; i!=ts_i_.size(); ++i) {
    if (y_[ts_i_[i]]>=0) {
      for (uint j=0; j!=ts_i_.size(); ++j) {
	if (y_[ts_i_[j]]<0) {
	  d.push_back(dec_values[i]-dec_values[j]);
	}
      }
    }
  }

  double avg, var;
  calculate_avg_var(d, avg, var);
  if (var<1e-10) var=1e-10;
  double rho = sqrt(var);
  double s  = SIGMOID_CONST/rho;
  double s2 = -SIGMOID_CONST/(d.size()*rho*var);

  double auc=0.0;
  delta.resize(ts_i_.size());
  std::fill(delta.begin(), delta.end(), 0.0);
  std::list<double>::const_iterator v=d.begin(); 
  for (uint i=0; i!=ts_i_.size(); ++i) {
    if (y_[ts_i_[i]]>=0) {
      for (uint j=0; j!=ts_i_.size(); ++j) {
	if (y_[ts_i_[j]]<0) {
	  double sig = sigmoid(*v, s);
	  auc += sig/d.size();
	  double w = sig*(1.0-sig)*(s+(*v)*s2*(*v-avg))/d.size();
	  delta[i] += w;
	  delta[j] -= w;
	  ++v;
	}
      }
    }
  }
  assert(v==d.end());
  
  return auc;
}

static
void
svm_train(const svm_problem *prob, double C, double eps, double *alpha, double& b)
{
  struct svm_parameter param;
  param.svm_type = C_SVC;
  param.kernel_type = PRECOMPUTED;
  param.degree = 3;
  param.gamma = 0;
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 40;
  param.C = C;
  param.probability = 1;
  param.eps = eps;
  param.p = 0.1;
  param.shrinking = 1;
  param.nr_weight = 0;
  
  int l = prob->l;
  double *minus_ones = new double[l];
  schar *y = new signed char[l];
  Solver::SolutionInfo si;
  int i;

  for(i=0;i<l;i++) {
    alpha[i] = 0;
    minus_ones[i] = -1;
    y[i] = prob->y[i]>0.0 ? +1 : -1;
  }

  Solver s;
  s.Solve(l, SVC_Q(*prob,param,y), minus_ones, y,
	  alpha, C, C, param.eps, &si, param.shrinking);
  b = si.rho;

  delete[] minus_ones;
  delete[] y;
}

static
double
svm_predict(const svm_problem *prob, const svm_node *x,
	    const double* alpha, double b)
{
  struct svm_parameter param;
  param.svm_type = C_SVC;
  param.kernel_type = PRECOMPUTED;
  param.degree = 3;
  param.gamma = 0;
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 40;
  param.C = 1;
  param.probability = 1;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.nr_weight = 0;

  int l = prob->l;
  double sum=0.0;
  for (int i=0; i<l; ++i) {
    if (alpha[i]!=0.0) {
      int y = prob->y[i]>0.0 ? +1 : -1;
      sum += alpha[i] * y * Kernel::k_function(x, prob->x[i], param);
    }
  }
  sum -= b;
  return sum;
}

template < class K >
static
struct svm_node*
make_node(const K& k, uint i)
{
  uint l=k[i].size();
  struct svm_node* x = new svm_node[l+2];
  uint c=0;
  x[c].index = 0;
  x[c].value = i+1;
  c++;
  for (uint j=0; j!=l; ++j) {
    x[c].index = j+1;
    x[c].value = k[i][j];
    c++;
  }
  x[c].index = -1;
  return x;
}

static
void
destroy_node(struct svm_node* x)
{
  delete[] x;
}

template < class K >
static
struct svm_problem*
make_problem(const std::vector<uint>& tr_i,
	     const std::vector<int>& y,
	     const K& k)
{
  struct svm_problem* prob = new struct svm_problem;
  prob->l = tr_i.size();
  prob->y = new double[prob->l];
  prob->x = new svm_node*[prob->l];
  for (uint i=0; i!=tr_i.size(); ++i) {
    prob->y[i] = y[tr_i[i]];
    prob->x[i] = make_node(k, tr_i[i]);
  }
  return prob;
}

static
void
destroy_problem(struct svm_problem* prob)
{
  for (int i=0; i!=prob->l; ++i)
    destroy_node(prob->x[i]);
  delete[] prob->x;
  delete[] prob->y;
  delete prob;
}

static inline
double
sigmoid(double x, double s)
{
  return 1.0/(1.0+exp(-s*x));
}

template < class T >
static inline
int 
sign(T v)
{
  return v>=0.0 ? +1 : -1;
}

template < class V, class W >
static
void
calculate_avg_var(const V& v, W& avg, W& var)
{
  double sum=0.0, sum2=0.0;
  uint n=v.size();
  typename boost::range_const_iterator<V>::type x;
  for (x=boost::begin(v); x!=boost::end(v); ++x) {
    sum += *x;
    sum2+= *x * *x;
  }
  avg = sum/n;
  var = sum2/n - avg*avg;
  if (var<0.0) var=0.0;		// for floating point error
}

static
void
find_support_vectors(const std::vector<uint>& tr_i,
		     double C,
		     const std::vector<double>& alpha,
		     bnu::vector<double>& alpha_u,
		     std::vector<uint>& u_idx,
		     std::vector<uint>& c_idx,
		     std::vector<uint>& z_idx)
{
  assert(tr_i.size()==alpha.size());

  uint nsv=0;
  uint csv=0;
  for (uint i=0; i!=tr_i.size(); ++i) {
    if (0.0<alpha[i] && alpha[i]<C)
      nsv++;
    else if (alpha[i]>=C)
      csv++;
  }
  alpha_u.resize(nsv);
  u_idx.resize(nsv);
  c_idx.resize(csv);
  z_idx.resize(tr_i.size()-nsv-csv);
  for (uint i=0, u=0, c=0, z=0; i!=tr_i.size(); ++i) {
    if (0.0<alpha[i] && alpha[i]<C) {
      u_idx[u] = tr_i[i];
      alpha_u[u] = alpha[i];
      u++;
    } else if (alpha[i]>=C) {
      c_idx[c++] = tr_i[i];
    } else {
      z_idx[z++] = tr_i[i];
    }
  }
}

template < class K >
static
void
solve_d(const K& kmat, const std::vector<int>& yvec,
	const std::vector<uint>& ts_i,
	const std::vector<uint>& u_idx, 
	const std::vector<double>& delta,
	bnu::vector<double>& d_u)
{
  typedef bnu::matrix<double> dmatrix;
  typedef bnu::vector<double> dvector;

  uint nsv = u_idx.size();
  if (nsv==0) return;

  // prepare P and r
  dmatrix p_u(nsv+1, nsv+1);
  dvector r_u(nsv+1);
  for (uint i=0; i!=nsv; ++i) {
    int yi = sign(yvec[u_idx[i]]);
    for (uint j=0; j!=nsv; ++j) {
      int yj = sign(yvec[u_idx[j]]);
      p_u(i,j) = yi*yj*kmat[u_idx[i]][u_idx[j]];
    }
    p_u(i,nsv)=p_u(nsv,i)=-yi;
    r_u(i)=0.0;
    for (uint j=0; j!=ts_i.size(); ++j) {
      r_u(i) += delta[j] * yi * kmat[u_idx[i]][ts_i[j]];
    }
  }

  p_u(nsv,nsv)=0.0;
  r_u(nsv)=0.0;
  for (uint j=0; j!=ts_i.size(); ++j) {
    r_u(nsv) -= delta[j];
  }

  // solve P d = r
#if 0
  try {
    bnu::permutation_matrix<> pm(p_u.size1());
    bnu::lu_factorize(p_u,pm);
    bnu::lu_substitute(p_u,pm,r_u);
    d_u = r_u;
  } catch (...) {
    std::cerr << "nsv=" << nsv << std::endl;
    throw;
  }
#else
  d_u = bnu::zero_vector<double>(nsv+1);
  conjugate_gradient(p_u, r_u, d_u, 1e-10);
#endif

#if 0
  dmatrix p_c(nsv+1, c_idx.size());
  dvector r_c(c_idx.size());
  for (uint i=0; i!=c_idx.size(); ++i) {
    int yi = sign(yvec[c_idx[i]]);

    for (uint j=0; j!=u_idx.size(); ++j) {
      int yj = sign(yvec[u_idx[j]]);
      p_c(j,i) = yi*yj*kmat[u_idx[j]][c_idx[i]];
    }
    p_c(nsv,i) = yi;

    r_c(i)=0.0;
    for (uint j=0; j!=ts_i.size(); ++j) {
      r_c(i) += delta[j] * yi * kmat[c_idx[i]][ts_i[j]];
    }
  }
  dvector d_c = r_c - bnu::prod(p_c, r_u);

  // copy the result and transform the indices
  std::vector<std::pair<uint,uint> > d_idx(alpha.size());
  dvector dx(alpha.size());
  uint j=0;
  for (uint i=0; i!=u_idx.size(); ++i) {
    dx(j)=r_u(i);
    d_idx[j]=std::make_pair(u_idx[i], j);
    j++;
  }
  for (uint i=0; i!=c_idx.size(); ++i) {
    dx(j)=d_c(i);
    d_idx[j]=std::make_pair(c_idx[i], j);
    j++;
  }
  for (uint i=0; i!=z_idx.size(); ++i) {
    dx(j)=0.0;
    d_idx[j]=std::make_pair(z_idx[i], j);
    j++;
  }
  std::sort(d_idx.begin(), d_idx.end());

  std::vector<std::pair<uint,uint> > tr_idx(tr_i.size());
  for (uint i=0; i!=tr_i.size(); ++i)
    tr_idx[i]=std::make_pair(tr_i[i], i);
  std::sort(tr_idx.begin(), tr_idx.end());
  
  d.resize(alpha.size());
  for (uint i=0; i!=d.size(); ++i) {
    assert(d_idx[i].first==tr_idx[i].first);
    d(tr_idx[i].second) = dx(d_idx[i].second);
  }
#endif
}

template < class K >
double
calculate_gradient_c(const K& k,
		     const std::vector<int>& y,
		     const std::vector<uint>& u_idx,
		     const std::vector<uint>& c_idx,
		     const std::vector<uint>& ts_i,
		     const bnu::vector<double>& d_u,
		     const std::vector<double>& delta)
{
  uint nsv=u_idx.size();
  double cg = 0.0;
  if (nsv>0) {
    bnu::vector<double> q_u_dot(nsv+1);
    for (uint i=0; i!=u_idx.size(); ++i) {
      int yi = sign(y[u_idx[i]]);
      q_u_dot(i) = 0.0;
      for (uint j=0; j!=c_idx.size(); ++j) {
	int yj = sign(y[c_idx[j]]);
	q_u_dot(i) -= yi*yj*k[u_idx[i]][c_idx[j]];
      }
    }
    q_u_dot(nsv) = 0.0;
    for (uint j=0; j!=c_idx.size(); ++j) {
      q_u_dot(nsv) += sign(y[c_idx[j]]);
    }
    cg += bnu::inner_prod(d_u, q_u_dot);
  }

  for (uint i=0; i!=delta.size(); ++i) {
    for (uint j=0; j!=c_idx.size(); ++j) {
      int yj = sign(y[c_idx[j]]);
      cg += delta[i] * yj * k[ts_i[i]][c_idx[j]];
    }
  }
  return cg;
}

template < class K, class G >
double
calculate_gradient_p(const K& k,
		     const G& g,
		     const std::vector<int>& y,
		     uint p, double C,
		     const std::vector<uint>& u_idx,
		     const std::vector<uint>& c_idx,
		     const std::vector<uint>& tr_i,
		     const std::vector<uint>& ts_i,
		     const bnu::vector<double>& d_u,
		     const std::vector<double>& delta,
		     const std::vector<double>& alpha,
		     const bnu::vector<double>& alpha_u,
		     double b)
{
  uint nsv=u_idx.size();
  double ret=0.0;

  if (nsv>0) {
    bnu::vector<double> q_u_dot(nsv+1);
    for (uint i=0; i!=u_idx.size(); ++i) {
      int yi = sign(y[u_idx[i]]);
      q_u_dot(i) = 0.0;
      for (uint j=0; j!=c_idx.size(); ++j) {
	int yj = sign(y[c_idx[j]]);
	q_u_dot(i) -= yi*yj*g[p][u_idx[i]][c_idx[j]]*C;
      }
    }
    q_u_dot(nsv) = 0.0;

    bnu::matrix<double> p_u_dot(nsv+1, nsv+1);
    for (uint i=0; i!=u_idx.size(); ++i) {
      int yi = sign(y[u_idx[i]]);
      for (uint j=0; j!=u_idx.size(); ++j) {
	int yj = sign(y[u_idx[j]]);
	p_u_dot(i,j) = yi*yj*g[p][u_idx[i]][u_idx[j]];
      }
      p_u_dot(i,nsv)=p_u_dot(nsv,i)=0.0;
    }
    p_u_dot(nsv,nsv)=0.0;

    bnu::vector<double> beta_u(nsv+1);
    for (uint i=0; i!=u_idx.size(); ++i) {
      beta_u(i)=alpha_u(i);
    }
    beta_u(nsv)=b;

    ret += bnu::inner_prod(d_u, q_u_dot - bnu::prod(p_u_dot, beta_u));
  }

  bnu::vector<double> dpsi_dot(tr_i.size()+1);
  for (uint j=0; j!=tr_i.size(); ++j) {
    dpsi_dot(j)=0.0;
    int yj = sign(y[tr_i[j]]);
    for (uint i=0; i!=delta.size(); ++i) {
      dpsi_dot(j) += delta[i] * yj * g[p][ts_i[i]][tr_i[j]];
    }
  }
  dpsi_dot(tr_i.size())=0.0;

  bnu::vector<double> beta(alpha.size()+1);
  for (uint i=0; i!=tr_i.size(); ++i) {
    beta(i)=alpha[i];
  }
  beta(alpha.size())=b;

  ret += bnu::inner_prod(dpsi_dot, beta);

  return ret;
    
}

// solve Ax=b or minimize ||Ax-b|| for a symmetric A
template < class T >
void
conjugate_gradient(const bnu::matrix<T>& A, const bnu::vector<T>& b,
		   bnu::vector<T>& x, double tol /*=1e-10*/)
{
  using namespace boost::numeric::ublas;
  vector<T> r = b - prod(A,x);
  if (inner_prod(r,r)<tol) return;
  vector<T> w = -r;
  vector<T> z = prod(A,w);
  T a = inner_prod(r,w)/inner_prod(w,z);
  x += a*w;
  for (uint i=0; i!=A.size1(); ++i) {
    r -= a*z;
    if (inner_prod(r,r)<tol) break;
    T B = inner_prod(r,z)/inner_prod(w,z);
    w = -r + B*w;
    z = prod(A,w);
    a = inner_prod(r,w)/inner_prod(w,z);
    x += a*w;
  }
}
