// $Id$

#ifndef __INC_QMATRIX_H__
#define __INC_QMATRIX_H__

#include "../libsvm/svm.h"

typedef float Qfloat;
typedef signed char schar;
class Cache;

//
// Kernel evaluation
//
// the static method k_function is for doing single kernel evaluation
// the constructor of Kernel prepares to calculate the l*l kernel matrix
// the member function get_Q is for getting one column from the Q Matrix
//
class QMatrix
{
public:
  virtual Qfloat *get_Q(int column, int len) const = 0;
  virtual Qfloat *get_QD() const = 0;
  virtual void swap_index(int i, int j) const = 0;
  virtual ~QMatrix() {}
};

class Kernel : public QMatrix
{
public:
  Kernel(int l, svm_node * const * x, const svm_parameter& param);
  virtual ~Kernel();

  static double k_function(const svm_node *x, const svm_node *y,
			   const svm_parameter& param);
  virtual Qfloat *get_Q(int column, int len) const = 0;
  virtual Qfloat *get_QD() const = 0;
  virtual void swap_index(int i, int j) const;	// no so const...

protected:
  double (Kernel::*kernel_function)(int i, int j) const;

private:
  const svm_node **x;
  double *x_square;

  // svm_parameter
  const int kernel_type;
  const int degree;
  const double gamma;
  const double coef0;

  static double dot(const svm_node *px, const svm_node *py);
  double kernel_linear(int i, int j) const;
  double kernel_poly(int i, int j) const;
  double kernel_rbf(int i, int j) const;
  double kernel_sigmoid(int i, int j) const;
  double kernel_precomputed(int i, int j) const;
};

//
// Q matrices for various formulations
//
class SVC_Q: public Kernel
{ 
public:
  SVC_Q(const svm_problem& prob, const svm_parameter& param, const schar *y);
  Qfloat *get_Q(int i, int len) const;
  Qfloat *get_QD() const { return QD; }
  void swap_index(int i, int j) const;
  ~SVC_Q();

private:
  schar *y;
  Cache *cache;
  Qfloat *QD;
};

class ONE_CLASS_Q: public Kernel
{
public:
  ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param);
  Qfloat *get_Q(int i, int len) const;
  Qfloat *get_QD() const { return QD; }
  void swap_index(int i, int j) const;
  ~ONE_CLASS_Q();
  
private:
  Cache *cache;
  Qfloat *QD;
};

class SVR_Q: public Kernel
{ 
public:
  SVR_Q(const svm_problem& prob, const svm_parameter& param);
  void swap_index(int i, int j) const;
  Qfloat *get_Q(int i, int len) const;
  Qfloat *get_QD() const { return QD; }
  ~SVR_Q();

private:
  int l;
  Cache *cache;
  schar *sign;
  int *index;
  mutable int next_buffer;
  Qfloat *buffer[2];
  Qfloat *QD;
};

#endif	//  __INC_QMATRIX_H__
