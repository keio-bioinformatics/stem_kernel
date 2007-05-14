/* $Id$ */

#ifndef __INC_LBFGS_H__
#define __INC_LBFGS_H__

#include <cfloat>
#include <cstring>

extern "C" {
  int lbfgs_(long *n, long *m, double *x, double *f, double *g,
	     long *diagco, double *diag, long *iprint, double *eps,
	     double *xtol, double *w, long *iflag);
  int setulb_(long *n, long *m, double *x, double *l, double *u, long *nbd,
	      double *f, double *g, double *factr, double *pgtol,
	      double *wa, long *iwa, char *task, long *iprint,
	      char *csave, long *lsave, long *isave, double *dsave,
	      short task_len, short csave_len);
};

class LBFGS
{
private:
  long n_;
  long m_;
  double *diag_;
  double eps_;
  double xtol_;
  double *w_;
  long iflag_;
  
public:
  LBFGS() : n_(0), m_(0), diag_(0),
	    eps_(DBL_EPSILON), xtol_(DBL_EPSILON),
	    w_(0), iflag_(0)
  {}
  
  ~LBFGS()
  {
    delete[] diag_;
    delete[] w_;
  }

  void initialize(int n, int m)
  {
    if (diag_ && n!=n_) delete[] diag_;
    diag_ = new double[n];
    if (w_ && (n!=n_ || m!=m_)) delete[] w_;
    w_=new double[n*(2*m+1)+2*m];
    n_ = n;
    m_ = m;
    iflag_ = 0;
  }
  
  int update(double *x, double *f, double *g)
  {
    long iprint[] = {-1,0};
    long diagco = 0;
    lbfgs_(&n_, &m_, x, f, g, &diagco, diag_, iprint,
	   &eps_, &xtol_, w_, &iflag_);
    return iflag_;
  }
};

class LBFGSB
{
private:
  long n_;
  long m_;
  double *l_;
  double *u_;
  long *nbd_;
  double *wa_;
  long *iwa_;
  char task_[60];
  char csave_[60];
  long lsave_[4];
  long isave_[44];
  double dsave_[29];
  
public:
  LBFGSB()
    : n_(0), m_(0), l_(0), u_(0), nbd_(0), wa_(0), iwa_(0)
  {
  }
  
  ~LBFGSB()
  {
    delete[] l_;
    delete[] u_;
    delete[] nbd_;
    delete[] wa_;
    delete[] iwa_;
  }
  
  void initialize(long n, long m, double *l, double *u, long *nbd)
  {
    if (wa_ && (n!=n_||m!=m_))  delete[] wa_;
    wa_ = new double[2*m*n+4*n+11*m*m+8*m];
    if (iwa_ && n!=n_) delete[] iwa_;
    iwa_ = new long[3*n];
    if (l_ && n!=n_) delete[] l_;
    l_ = new double[n];
    if (u_ && n!=n_) delete[] u_;
    u_ = new double[n];
    if (nbd_ && n!=n_) delete[] nbd_;
    nbd_ = new long[n];
    for (int i=0; i<n; ++i) {
      l_[i] = l[i];
      u_[i] = u[i];
      nbd_[i] = nbd[i];
    }
    n_ = n;
    m_ = m;
    strcpy(task_, "START");
  }
    
  int update(double *x, double *f, double *g)
  {
    long iprint = 1;
    double factr = 1.0;
    double pgtol = 0.0;
    while (1) {
      setulb_(&n_, &m_, x, l_, u_, nbd_, f, g, &factr, &pgtol, wa_, iwa_,
	      task_, &iprint, csave_, lsave_, isave_, dsave_, 60, 60);
      if (strncmp(task_, "FG", 2)==0) return 1;
      else if (strncmp(task_, "NEW_X", 5)==0) continue;
      else if (strncmp(task_, "CONV", 5)==0) return 0;
      else return -1;
    }
  }
};

#endif /* __INC_LBFGS_H__ */

// Local Variables:
// mode:C++
// End:
