// $Id$

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "qmatrix.h"

using std::swap;
using std::min;
using std::max;

template <class S, class T> inline void clone(T*& dst, S* src, int n)
{
  dst = new T[n];
  memcpy((void *)dst,(void *)src,sizeof(T)*n);
}
inline double powi(double base, int times)
{
  double tmp = base, ret = 1.0;

  for(int t=times; t>0; t/=2)
    {
      if(t%2==1) ret*=tmp;
      tmp = tmp * tmp;
    }
  return ret;
}


//
// Kernel Cache
//
// l is the number of total data items
// size is the cache size limit in bytes
//
class Cache
{
public:
  Cache(int l,long int size);
  ~Cache();

  // request data [0,len)
  // return some position p where [p,len) need to be filled
  // (p >= len if nothing needs to be filled)
  int get_data(const int index, Qfloat **data, int len);
  void swap_index(int i, int j);	// future_option

private:
  int l;
  long int size;
  struct head_t
  {
    head_t *prev, *next;	// a cicular list
    Qfloat *data;
    int len;		// data[0,len) is cached in this entry
  };

  head_t *head;
  head_t lru_head;
  void lru_delete(head_t *h);
  void lru_insert(head_t *h);
};

Cache::
Cache(int l_,long int size_):l(l_),size(size_)
{
  head = (head_t *)calloc(l,sizeof(head_t));	// initialized to 0
  size /= sizeof(Qfloat);
  size -= l * sizeof(head_t) / sizeof(Qfloat);
  size = max(size, (long int) 2*l);	// cache must be large enough for two columns
  lru_head.next = lru_head.prev = &lru_head;
}

Cache::
~Cache()
{
  for(head_t *h = lru_head.next; h != &lru_head; h=h->next)
    free(h->data);
  free(head);
}

void
Cache::
lru_delete(head_t *h)
{
  // delete from current location
  h->prev->next = h->next;
  h->next->prev = h->prev;
}

void
Cache::
lru_insert(head_t *h)
{
  // insert to last position
  h->next = &lru_head;
  h->prev = lru_head.prev;
  h->prev->next = h;
  h->next->prev = h;
}

int
Cache::
get_data(const int index, Qfloat **data, int len)
{
  head_t *h = &head[index];
  if(h->len) lru_delete(h);
  int more = len - h->len;

  if(more > 0)
    {
      // free old space
      while(size < more)
	{
	  head_t *old = lru_head.next;
	  lru_delete(old);
	  free(old->data);
	  size += old->len;
	  old->data = 0;
	  old->len = 0;
	}

      // allocate new space
      h->data = (Qfloat *)realloc(h->data,sizeof(Qfloat)*len);
      size -= more;
      swap(h->len,len);
    }

  lru_insert(h);
  *data = h->data;
  return len;
}

void
Cache::
swap_index(int i, int j)
{
  if(i==j) return;

  if(head[i].len) lru_delete(&head[i]);
  if(head[j].len) lru_delete(&head[j]);
  swap(head[i].data,head[j].data);
  swap(head[i].len,head[j].len);
  if(head[i].len) lru_insert(&head[i]);
  if(head[j].len) lru_insert(&head[j]);

  if(i>j) swap(i,j);
  for(head_t *h = lru_head.next; h!=&lru_head; h=h->next)
    {
      if(h->len > i)
	{
	  if(h->len > j)
	    swap(h->data[i],h->data[j]);
	  else
	    {
	      // give up
	      lru_delete(h);
	      free(h->data);
	      size += h->len;
	      h->data = 0;
	      h->len = 0;
	    }
	}
    }
}

Kernel::
Kernel(int l, svm_node * const * x_, const svm_parameter& param)
 : kernel_type(param.kernel_type), degree(param.degree),
   gamma(param.gamma), coef0(param.coef0)
{
  switch(kernel_type)
    {
    case LINEAR:
      kernel_function = &Kernel::kernel_linear;
      break;
    case POLY:
      kernel_function = &Kernel::kernel_poly;
      break;
    case RBF:
      kernel_function = &Kernel::kernel_rbf;
      break;
    case SIGMOID:
      kernel_function = &Kernel::kernel_sigmoid;
      break;
    case PRECOMPUTED:
      kernel_function = &Kernel::kernel_precomputed;
      break;
    }

  clone(x,x_,l);

  if(kernel_type == RBF)
    {
      x_square = new double[l];
      for(int i=0;i<l;i++)
	x_square[i] = dot(x[i],x[i]);
    }
  else
    x_square = 0;
}

Kernel::
~Kernel()
{
  delete[] x;
  delete[] x_square;
}

void
Kernel::
swap_index(int i, int j) const	// no so const...
{
  swap(x[i],x[j]);
  if(x_square) swap(x_square[i],x_square[j]);
}

double
Kernel::
dot(const svm_node *px, const svm_node *py)
{
  double sum = 0;
  while(px->index != -1 && py->index != -1)
    {
      if(px->index == py->index)
	{
	  sum += px->value * py->value;
	  ++px;
	  ++py;
	}
      else
	{
	  if(px->index > py->index)
	    ++py;
	  else
	    ++px;
	}			
    }
  return sum;
}

double
Kernel::
k_function(const svm_node *x, const svm_node *y, const svm_parameter& param)
{
  switch(param.kernel_type)
    {
    case LINEAR:
      return dot(x,y);
    case POLY:
      return powi(param.gamma*dot(x,y)+param.coef0,param.degree);
    case RBF:
      {
	double sum = 0;
	while(x->index != -1 && y->index !=-1)
	  {
	    if(x->index == y->index)
	      {
		double d = x->value - y->value;
		sum += d*d;
		++x;
		++y;
	      }
	    else
	      {
		if(x->index > y->index)
		  {	
		    sum += y->value * y->value;
		    ++y;
		  }
		else
		  {
		    sum += x->value * x->value;
		    ++x;
		  }
	      }
	  }

	while(x->index != -1)
	  {
	    sum += x->value * x->value;
	    ++x;
	  }

	while(y->index != -1)
	  {
	    sum += y->value * y->value;
	    ++y;
	  }
			
	return exp(-param.gamma*sum);
      }
    case SIGMOID:
      return tanh(param.gamma*dot(x,y)+param.coef0);
    case PRECOMPUTED:  //x: test (validation), y: SV
      return x[(int)(y->value)].value;
    default:
      return 0;	/* Unreachable */
    }
}

double
Kernel::
kernel_linear(int i, int j) const
{
  return dot(x[i],x[j]);
}

double
Kernel::
kernel_poly(int i, int j) const
{
  return powi(gamma*dot(x[i],x[j])+coef0,degree);
}

double
Kernel::
kernel_rbf(int i, int j) const
{
  return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
}

double
Kernel::
kernel_sigmoid(int i, int j) const
{
  return tanh(gamma*dot(x[i],x[j])+coef0);
}

double
Kernel::
kernel_precomputed(int i, int j) const
{
  return x[i][(int)(x[j][0].value)].value;
}


SVC_Q::
SVC_Q(const svm_problem& prob, const svm_parameter& param, const schar *y_)
  : Kernel(prob.l, prob.x, param)
{
  clone(y,y_,prob.l);
  cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
  QD = new Qfloat[prob.l];
  for(int i=0;i<prob.l;i++)
    QD[i]= (Qfloat)(this->*kernel_function)(i,i);
}

Qfloat*
SVC_Q::
get_Q(int i, int len) const
{
  Qfloat *data;
  int start;
  if((start = cache->get_data(i,&data,len)) < len)
    {
      for(int j=start;j<len;j++)
	data[j] = (Qfloat)(y[i]*y[j]*(this->*kernel_function)(i,j));
    }
  return data;
}

void
SVC_Q::
swap_index(int i, int j) const
{
  cache->swap_index(i,j);
  Kernel::swap_index(i,j);
  swap(y[i],y[j]);
  swap(QD[i],QD[j]);
}

SVC_Q::
~SVC_Q()
{
  delete[] y;
  delete cache;
  delete[] QD;
}

ONE_CLASS_Q::
ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param)
  : Kernel(prob.l, prob.x, param)
{
  cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
  QD = new Qfloat[prob.l];
  for(int i=0;i<prob.l;i++)
    QD[i]= (Qfloat)(this->*kernel_function)(i,i);
}
	
Qfloat*
ONE_CLASS_Q::
get_Q(int i, int len) const
{
  Qfloat *data;
  int start;
  if((start = cache->get_data(i,&data,len)) < len)
    {
      for(int j=start;j<len;j++)
	data[j] = (Qfloat)(this->*kernel_function)(i,j);
    }
  return data;
}

void
ONE_CLASS_Q::
swap_index(int i, int j) const
{
  cache->swap_index(i,j);
  Kernel::swap_index(i,j);
  swap(QD[i],QD[j]);
}

ONE_CLASS_Q::
~ONE_CLASS_Q()
{
  delete cache;
  delete[] QD;
}


SVR_Q::
SVR_Q(const svm_problem& prob, const svm_parameter& param)
  : Kernel(prob.l, prob.x, param)
{
  l = prob.l;
  cache = new Cache(l,(long int)(param.cache_size*(1<<20)));
  QD = new Qfloat[2*l];
  sign = new schar[2*l];
  index = new int[2*l];
  for(int k=0;k<l;k++)
    {
      sign[k] = 1;
      sign[k+l] = -1;
      index[k] = k;
      index[k+l] = k;
      QD[k]= (Qfloat)(this->*kernel_function)(k,k);
      QD[k+l]=QD[k];
    }
  buffer[0] = new Qfloat[2*l];
  buffer[1] = new Qfloat[2*l];
  next_buffer = 0;
}

void
SVR_Q::
swap_index(int i, int j) const
{
  swap(sign[i],sign[j]);
  swap(index[i],index[j]);
  swap(QD[i],QD[j]);
}


Qfloat*
SVR_Q::
get_Q(int i, int len) const
{
  Qfloat *data;
  int real_i = index[i];
  if(cache->get_data(real_i,&data,l) < l)
    {
      for(int j=0;j<l;j++)
	data[j] = (Qfloat)(this->*kernel_function)(real_i,j);
    }

  // reorder and copy
  Qfloat *buf = buffer[next_buffer];
  next_buffer = 1 - next_buffer;
  schar si = sign[i];
  for(int j=0;j<len;j++)
    buf[j] = si * sign[j] * data[index[j]];
  return buf;
}

SVR_Q::
~SVR_Q()
{
  delete cache;
  delete[] sign;
  delete[] index;
  delete[] buffer[0];
  delete[] buffer[1];
  delete[] QD;
}

