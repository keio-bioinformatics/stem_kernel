// $Id$

#ifndef __INC_SOLVER_H__
#define __INC_SOLVER_H__

#include "qmatrix.h"

typedef float Qfloat;
typedef signed char schar;

// An SMO algorithm in Fan et al., JMLR 6(2005), p. 1889--1918
// Solves:
//
//	min 0.5(\alpha^T Q \alpha) + p^T \alpha
//
//		y^T \alpha = \delta
//		y_i = +1 or -1
//		0 <= alpha_i <= Cp for y_i = 1
//		0 <= alpha_i <= Cn for y_i = -1
//
// Given:
//
//	Q, p, y, Cp, Cn, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping tolerance
//
// solution will be put in \alpha, objective value will be put in obj
//
class Solver {
public:
  Solver() {};
  virtual ~Solver() {};

  struct SolutionInfo {
    double obj;
    double rho;
    double upper_bound_p;
    double upper_bound_n;
    double r;	// for Solver_NU
  };

  void Solve(int l, const QMatrix& Q, const double *p_, const schar *y_,
	     double *alpha_, double Cp, double Cn, double eps,
	     SolutionInfo* si, int shrinking);

protected:
  int active_size;
  schar *y;
  double *G;		// gradient of objective function
  enum { LOWER_BOUND, UPPER_BOUND, FREE };
  char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
  double *alpha;
  const QMatrix *Q;
  const Qfloat *QD;
  double eps;
  double Cp,Cn;
  double *p;
  int *active_set;
  double *G_bar;		// gradient, if we treat free variables as 0
  int l;
  bool unshrinked;	// XXX

  double get_C(int i)
  {
    return (y[i] > 0)? Cp : Cn;
  }
  void update_alpha_status(int i)
  {
    if(alpha[i] >= get_C(i))
      alpha_status[i] = UPPER_BOUND;
    else if(alpha[i] <= 0)
      alpha_status[i] = LOWER_BOUND;
    else alpha_status[i] = FREE;
  }
  bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
  bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
  bool is_free(int i) { return alpha_status[i] == FREE; }
  void swap_index(int i, int j);
  void reconstruct_gradient();
  virtual int select_working_set(int &i, int &j);
  virtual double calculate_rho();
  virtual void do_shrinking();

private:
  bool be_shrunken(int i, double Gmax1, double Gmax2);	
};

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//
class Solver_NU : public Solver
{
public:
  Solver_NU() {}
  void Solve(int l, const QMatrix& Q, const double *p, const schar *y,
	     double *alpha, double Cp, double Cn, double eps,
	     SolutionInfo* si, int shrinking)
  {
    this->si = si;
    Solver::Solve(l,Q,p,y,alpha,Cp,Cn,eps,si,shrinking);
  }

private:
  SolutionInfo *si;
  int select_working_set(int &i, int &j);
  double calculate_rho();
  bool be_shrunken(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4);
  void do_shrinking();
};

#endif	// __INC_SOLVER_H__
