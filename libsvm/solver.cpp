// $Id$

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <cmath>
#include <cctype>
#include <cfloat>
#include "solver.h"

using std::swap;
using std::min;
using std::max;

template <class S, class T> inline void clone(T*& dst, S* src, int n)
{
  dst = new T[n];
  memcpy((void *)dst,(void *)src,sizeof(T)*n);
}

#define INF HUGE_VAL
#define TAU 1e-12
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

#if 0
static void info(const char *fmt,...)
{
  va_list ap;
  va_start(ap,fmt);
  vprintf(fmt,ap);
  va_end(ap);
}
static void info_flush()
{
  fflush(stdout);
}
#else
static void info(const char *fmt,...) {}
static void info_flush() {}
#endif

// Solver
void
Solver::
swap_index(int i, int j)
{
  Q->swap_index(i,j);
  swap(y[i],y[j]);
  swap(G[i],G[j]);
  swap(alpha_status[i],alpha_status[j]);
  swap(alpha[i],alpha[j]);
  swap(p[i],p[j]);
  swap(active_set[i],active_set[j]);
  swap(G_bar[i],G_bar[j]);
}

void
Solver::
reconstruct_gradient()
{
  // reconstruct inactive elements of G from G_bar and free variables

  if(active_size == l) return;

  int i;
  for(i=active_size;i<l;i++)
    G[i] = G_bar[i] + p[i];
	
  for(i=0;i<active_size;i++)
    if(is_free(i))
      {
	const Qfloat *Q_i = Q->get_Q(i,l);
	double alpha_i = alpha[i];
	for(int j=active_size;j<l;j++)
	  G[j] += alpha_i * Q_i[j];
      }
}

void
Solver::
Solve(int l, const QMatrix& Q, const double *p_, const schar *y_,
      double *alpha_, double Cp, double Cn, double eps,
      SolutionInfo* si, int shrinking)
{
  this->l = l;
  this->Q = &Q;
  QD=Q.get_QD();
  clone(p, p_,l);
  clone(y, y_,l);
  clone(alpha,alpha_,l);
  this->Cp = Cp;
  this->Cn = Cn;
  this->eps = eps;
  unshrinked = false;

  // initialize alpha_status
  {
    alpha_status = new char[l];
    for(int i=0;i<l;i++)
      update_alpha_status(i);
  }

  // initialize active set (for shrinking)
  {
    active_set = new int[l];
    for(int i=0;i<l;i++)
      active_set[i] = i;
    active_size = l;
  }

  // initialize gradient
  {
    G = new double[l];
    G_bar = new double[l];
    int i;
    for(i=0;i<l;i++)
      {
	G[i] = p[i];
	G_bar[i] = 0;
      }
    for(i=0;i<l;i++)
      if(!is_lower_bound(i))
	{
	  const Qfloat *Q_i = Q.get_Q(i,l);
	  double alpha_i = alpha[i];
	  int j;
	  for(j=0;j<l;j++)
	    G[j] += alpha_i*Q_i[j];
	  if(is_upper_bound(i))
	    for(j=0;j<l;j++)
	      G_bar[j] += get_C(i) * Q_i[j];
	}
  }

  // optimization step

  int iter = 0;
  int counter = min(l,1000)+1;

  while(1)
    {
      // show progress and do shrinking

      if(--counter == 0)
	{
	  counter = min(l,1000);
	  if(shrinking) do_shrinking();
	  info("."); info_flush();
	}

      int i,j;
      if(select_working_set(i,j)!=0)
	{
	  // reconstruct the whole gradient
	  reconstruct_gradient();
	  // reset active set size and check
	  active_size = l;
	  info("*"); info_flush();
	  if(select_working_set(i,j)!=0)
	    break;
	  else
	    counter = 1;	// do shrinking next iteration
	}
		
      ++iter;

      // update alpha[i] and alpha[j], handle bounds carefully
		
      const Qfloat *Q_i = Q.get_Q(i,active_size);
      const Qfloat *Q_j = Q.get_Q(j,active_size);

      double C_i = get_C(i);
      double C_j = get_C(j);

      double old_alpha_i = alpha[i];
      double old_alpha_j = alpha[j];

      if(y[i]!=y[j])
	{
	  double quad_coef = Q_i[i]+Q_j[j]+2*Q_i[j];
	  if (quad_coef <= 0)
	    quad_coef = TAU;
	  double delta = (-G[i]-G[j])/quad_coef;
	  double diff = alpha[i] - alpha[j];
	  alpha[i] += delta;
	  alpha[j] += delta;
			
	  if(diff > 0)
	    {
	      if(alpha[j] < 0)
		{
		  alpha[j] = 0;
		  alpha[i] = diff;
		}
	    }
	  else
	    {
	      if(alpha[i] < 0)
		{
		  alpha[i] = 0;
		  alpha[j] = -diff;
		}
	    }
	  if(diff > C_i - C_j)
	    {
	      if(alpha[i] > C_i)
		{
		  alpha[i] = C_i;
		  alpha[j] = C_i - diff;
		}
	    }
	  else
	    {
	      if(alpha[j] > C_j)
		{
		  alpha[j] = C_j;
		  alpha[i] = C_j + diff;
		}
	    }
	}
      else
	{
	  double quad_coef = Q_i[i]+Q_j[j]-2*Q_i[j];
	  if (quad_coef <= 0)
	    quad_coef = TAU;
	  double delta = (G[i]-G[j])/quad_coef;
	  double sum = alpha[i] + alpha[j];
	  alpha[i] -= delta;
	  alpha[j] += delta;

	  if(sum > C_i)
	    {
	      if(alpha[i] > C_i)
		{
		  alpha[i] = C_i;
		  alpha[j] = sum - C_i;
		}
	    }
	  else
	    {
	      if(alpha[j] < 0)
		{
		  alpha[j] = 0;
		  alpha[i] = sum;
		}
	    }
	  if(sum > C_j)
	    {
	      if(alpha[j] > C_j)
		{
		  alpha[j] = C_j;
		  alpha[i] = sum - C_j;
		}
	    }
	  else
	    {
	      if(alpha[i] < 0)
		{
		  alpha[i] = 0;
		  alpha[j] = sum;
		}
	    }
	}

      // update G

      double delta_alpha_i = alpha[i] - old_alpha_i;
      double delta_alpha_j = alpha[j] - old_alpha_j;
		
      for(int k=0;k<active_size;k++)
	{
	  G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
	}

      // update alpha_status and G_bar

      {
	bool ui = is_upper_bound(i);
	bool uj = is_upper_bound(j);
	update_alpha_status(i);
	update_alpha_status(j);
	int k;
	if(ui != is_upper_bound(i))
	  {
	    Q_i = Q.get_Q(i,l);
	    if(ui)
	      for(k=0;k<l;k++)
		G_bar[k] -= C_i * Q_i[k];
	    else
	      for(k=0;k<l;k++)
		G_bar[k] += C_i * Q_i[k];
	  }

	if(uj != is_upper_bound(j))
	  {
	    Q_j = Q.get_Q(j,l);
	    if(uj)
	      for(k=0;k<l;k++)
		G_bar[k] -= C_j * Q_j[k];
	    else
	      for(k=0;k<l;k++)
		G_bar[k] += C_j * Q_j[k];
	  }
      }
    }

  // calculate rho

  si->rho = calculate_rho();

  // calculate objective value
  {
    double v = 0;
    int i;
    for(i=0;i<l;i++)
      v += alpha[i] * (G[i] + p[i]);

    si->obj = v/2;
  }

  // put back the solution
  {
    for(int i=0;i<l;i++)
      alpha_[active_set[i]] = alpha[i];
  }

  // juggle everything back
  /*{
    for(int i=0;i<l;i++)
    while(active_set[i] != i)
    swap_index(i,active_set[i]);
    // or Q.swap_index(i,active_set[i]);
    }*/

  si->upper_bound_p = Cp;
  si->upper_bound_n = Cn;

  info("\noptimization finished, #iter = %d\n",iter);

  delete[] p;
  delete[] y;
  delete[] alpha;
  delete[] alpha_status;
  delete[] active_set;
  delete[] G;
  delete[] G_bar;
}

// return 1 if already optimal, return 0 otherwise
int
Solver::
select_working_set(int &out_i, int &out_j)
{
  // return i,j such that
  // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
  // j: minimizes the decrease of obj value
  //    (if quadratic coefficeint <= 0, replace it with tau)
  //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
	
  double Gmax = -INF;
  double Gmax2 = -INF;
  int Gmax_idx = -1;
  int Gmin_idx = -1;
  double obj_diff_min = INF;

  for(int t=0;t<active_size;t++)
    if(y[t]==+1)	
      {
	if(!is_upper_bound(t))
	  if(-G[t] >= Gmax)
	    {
	      Gmax = -G[t];
	      Gmax_idx = t;
	    }
      }
    else
      {
	if(!is_lower_bound(t))
	  if(G[t] >= Gmax)
	    {
	      Gmax = G[t];
	      Gmax_idx = t;
	    }
      }

  int i = Gmax_idx;
  const Qfloat *Q_i = NULL;
  if(i != -1) // NULL Q_i not accessed: Gmax=-INF if i=-1
    Q_i = Q->get_Q(i,active_size);

  for(int j=0;j<active_size;j++)
    {
      if(y[j]==+1)
	{
	  if (!is_lower_bound(j))
	    {
	      double grad_diff=Gmax+G[j];
	      if (G[j] >= Gmax2)
		Gmax2 = G[j];
	      if (grad_diff > 0)
		{
		  double obj_diff; 
		  double quad_coef=Q_i[i]+QD[j]-2*y[i]*Q_i[j];
		  if (quad_coef > 0)
		    obj_diff = -(grad_diff*grad_diff)/quad_coef;
		  else
		    obj_diff = -(grad_diff*grad_diff)/TAU;

		  if (obj_diff <= obj_diff_min)
		    {
		      Gmin_idx=j;
		      obj_diff_min = obj_diff;
		    }
		}
	    }
	}
      else
	{
	  if (!is_upper_bound(j))
	    {
	      double grad_diff= Gmax-G[j];
	      if (-G[j] >= Gmax2)
		Gmax2 = -G[j];
	      if (grad_diff > 0)
		{
		  double obj_diff; 
		  double quad_coef=Q_i[i]+QD[j]+2*y[i]*Q_i[j];
		  if (quad_coef > 0)
		    obj_diff = -(grad_diff*grad_diff)/quad_coef;
		  else
		    obj_diff = -(grad_diff*grad_diff)/TAU;

		  if (obj_diff <= obj_diff_min)
		    {
		      Gmin_idx=j;
		      obj_diff_min = obj_diff;
		    }
		}
	    }
	}
    }

  if(Gmax+Gmax2 < eps)
    return 1;

  out_i = Gmax_idx;
  out_j = Gmin_idx;
  return 0;
}

bool
Solver::
be_shrunken(int i, double Gmax1, double Gmax2)
{
  if(is_upper_bound(i))
    {
      if(y[i]==+1)
	return(-G[i] > Gmax1);
      else
	return(-G[i] > Gmax2);
    }
  else if(is_lower_bound(i))
    {
      if(y[i]==+1)
	return(G[i] > Gmax2);
      else	
	return(G[i] > Gmax1);
    }
  else
    return(false);
}

void
Solver::
do_shrinking()
{
  int i;
  double Gmax1 = -INF;		// max { -y_i * grad(f)_i | i in I_up(\alpha) }
  double Gmax2 = -INF;		// max { y_i * grad(f)_i | i in I_low(\alpha) }

  // find maximal violating pair first
  for(i=0;i<active_size;i++)
    {
      if(y[i]==+1)	
	{
	  if(!is_upper_bound(i))	
	    {
	      if(-G[i] >= Gmax1)
		Gmax1 = -G[i];
	    }
	  if(!is_lower_bound(i))	
	    {
	      if(G[i] >= Gmax2)
		Gmax2 = G[i];
	    }
	}
      else	
	{
	  if(!is_upper_bound(i))	
	    {
	      if(-G[i] >= Gmax2)
		Gmax2 = -G[i];
	    }
	  if(!is_lower_bound(i))	
	    {
	      if(G[i] >= Gmax1)
		Gmax1 = G[i];
	    }
	}
    }

  // shrink

  for(i=0;i<active_size;i++)
    if (be_shrunken(i, Gmax1, Gmax2))
      {
	active_size--;
	while (active_size > i)
	  {
	    if (!be_shrunken(active_size, Gmax1, Gmax2))
	      {
		swap_index(i,active_size);
		break;
	      }
	    active_size--;
	  }
      }

  // unshrink, check all variables again before final iterations

  if(unshrinked || Gmax1 + Gmax2 > eps*10) return;
	
  unshrinked = true;
  reconstruct_gradient();

  for(i=l-1;i>=active_size;i--)
    if (!be_shrunken(i, Gmax1, Gmax2))
      {
	while (active_size < i)
	  {
	    if (be_shrunken(active_size, Gmax1, Gmax2))
	      {
		swap_index(i,active_size);
		break;
	      }
	    active_size++;
	  }
	active_size++;
      }
}

double
Solver::
calculate_rho()
{
  double r;
  int nr_free = 0;
  double ub = INF, lb = -INF, sum_free = 0;
  for(int i=0;i<active_size;i++)
    {
      double yG = y[i]*G[i];

      if(is_upper_bound(i))
	{
	  if(y[i]==-1)
	    ub = min(ub,yG);
	  else
	    lb = max(lb,yG);
	}
      else if(is_lower_bound(i))
	{
	  if(y[i]==+1)
	    ub = min(ub,yG);
	  else
	    lb = max(lb,yG);
	}
      else
	{
	  ++nr_free;
	  sum_free += yG;
	}
    }

  if(nr_free>0)
    r = sum_free/nr_free;
  else
    r = (ub+lb)/2;

  return r;
}

// return 1 if already optimal, return 0 otherwise
int Solver_NU::select_working_set(int &out_i, int &out_j)
{
  // return i,j such that y_i = y_j and
  // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
  // j: minimizes the decrease of obj value
  //    (if quadratic coefficeint <= 0, replace it with tau)
  //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)

  double Gmaxp = -INF;
  double Gmaxp2 = -INF;
  int Gmaxp_idx = -1;

  double Gmaxn = -INF;
  double Gmaxn2 = -INF;
  int Gmaxn_idx = -1;

  int Gmin_idx = -1;
  double obj_diff_min = INF;

  for(int t=0;t<active_size;t++)
    if(y[t]==+1)
      {
	if(!is_upper_bound(t))
	  if(-G[t] >= Gmaxp)
	    {
	      Gmaxp = -G[t];
	      Gmaxp_idx = t;
	    }
      }
    else
      {
	if(!is_lower_bound(t))
	  if(G[t] >= Gmaxn)
	    {
	      Gmaxn = G[t];
	      Gmaxn_idx = t;
	    }
      }

  int ip = Gmaxp_idx;
  int in = Gmaxn_idx;
  const Qfloat *Q_ip = NULL;
  const Qfloat *Q_in = NULL;
  if(ip != -1) // NULL Q_ip not accessed: Gmaxp=-INF if ip=-1
    Q_ip = Q->get_Q(ip,active_size);
  if(in != -1)
    Q_in = Q->get_Q(in,active_size);

  for(int j=0;j<active_size;j++)
    {
      if(y[j]==+1)
	{
	  if (!is_lower_bound(j))	
	    {
	      double grad_diff=Gmaxp+G[j];
	      if (G[j] >= Gmaxp2)
		Gmaxp2 = G[j];
	      if (grad_diff > 0)
		{
		  double obj_diff; 
		  double quad_coef = Q_ip[ip]+QD[j]-2*Q_ip[j];
		  if (quad_coef > 0)
		    obj_diff = -(grad_diff*grad_diff)/quad_coef;
		  else
		    obj_diff = -(grad_diff*grad_diff)/TAU;

		  if (obj_diff <= obj_diff_min)
		    {
		      Gmin_idx=j;
		      obj_diff_min = obj_diff;
		    }
		}
	    }
	}
      else
	{
	  if (!is_upper_bound(j))
	    {
	      double grad_diff=Gmaxn-G[j];
	      if (-G[j] >= Gmaxn2)
		Gmaxn2 = -G[j];
	      if (grad_diff > 0)
		{
		  double obj_diff; 
		  double quad_coef = Q_in[in]+QD[j]-2*Q_in[j];
		  if (quad_coef > 0)
		    obj_diff = -(grad_diff*grad_diff)/quad_coef;
		  else
		    obj_diff = -(grad_diff*grad_diff)/TAU;

		  if (obj_diff <= obj_diff_min)
		    {
		      Gmin_idx=j;
		      obj_diff_min = obj_diff;
		    }
		}
	    }
	}
    }

  if(max(Gmaxp+Gmaxp2,Gmaxn+Gmaxn2) < eps)
    return 1;

  if (y[Gmin_idx] == +1)
    out_i = Gmaxp_idx;
  else
    out_i = Gmaxn_idx;
  out_j = Gmin_idx;

  return 0;
}

bool
Solver_NU::
be_shrunken(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4)
{
  if(is_upper_bound(i))
    {
      if(y[i]==+1)
	return(-G[i] > Gmax1);
      else	
	return(-G[i] > Gmax4);
    }
  else if(is_lower_bound(i))
    {
      if(y[i]==+1)
	return(G[i] > Gmax2);
      else	
	return(G[i] > Gmax3);
    }
  else
    return(false);
}

void Solver_NU::do_shrinking()
{
  double Gmax1 = -INF;	// max { -y_i * grad(f)_i | y_i = +1, i in I_up(\alpha) }
  double Gmax2 = -INF;	// max { y_i * grad(f)_i | y_i = +1, i in I_low(\alpha) }
  double Gmax3 = -INF;	// max { -y_i * grad(f)_i | y_i = -1, i in I_up(\alpha) }
  double Gmax4 = -INF;	// max { y_i * grad(f)_i | y_i = -1, i in I_low(\alpha) }

  // find maximal violating pair first
  int i;
  for(i=0;i<active_size;i++)
    {
      if(!is_upper_bound(i))
	{
	  if(y[i]==+1)
	    {
	      if(-G[i] > Gmax1) Gmax1 = -G[i];
	    }
	  else	if(-G[i] > Gmax4) Gmax4 = -G[i];
	}
      if(!is_lower_bound(i))
	{
	  if(y[i]==+1)
	    {	
	      if(G[i] > Gmax2) Gmax2 = G[i];
	    }
	  else	if(G[i] > Gmax3) Gmax3 = G[i];
	}
    }

  // shrinking

  for(i=0;i<active_size;i++)
    if (be_shrunken(i, Gmax1, Gmax2, Gmax3, Gmax4))
      {
	active_size--;
	while (active_size > i)
	  {
	    if (!be_shrunken(active_size, Gmax1, Gmax2, Gmax3, Gmax4))
	      {
		swap_index(i,active_size);
		break;
	      }
	    active_size--;
	  }
      }

  // unshrink, check all variables again before final iterations

  if(unshrinked || max(Gmax1+Gmax2,Gmax3+Gmax4) > eps*10) return;
	
  unshrinked = true;
  reconstruct_gradient();

  for(i=l-1;i>=active_size;i--)
    if (!be_shrunken(i, Gmax1, Gmax2, Gmax3, Gmax4))
      {
	while (active_size < i)
	  {
	    if (be_shrunken(active_size, Gmax1, Gmax2, Gmax3, Gmax4))
	      {
		swap_index(i,active_size);
		break;
	      }
	    active_size++;
	  }
	active_size++;
      }
}

double
Solver_NU::
calculate_rho()
{
  int nr_free1 = 0,nr_free2 = 0;
  double ub1 = INF, ub2 = INF;
  double lb1 = -INF, lb2 = -INF;
  double sum_free1 = 0, sum_free2 = 0;

  for(int i=0;i<active_size;i++)
    {
      if(y[i]==+1)
	{
	  if(is_upper_bound(i))
	    lb1 = max(lb1,G[i]);
	  else if(is_lower_bound(i))
	    ub1 = min(ub1,G[i]);
	  else
	    {
	      ++nr_free1;
	      sum_free1 += G[i];
	    }
	}
      else
	{
	  if(is_upper_bound(i))
	    lb2 = max(lb2,G[i]);
	  else if(is_lower_bound(i))
	    ub2 = min(ub2,G[i]);
	  else
	    {
	      ++nr_free2;
	      sum_free2 += G[i];
	    }
	}
    }

  double r1,r2;
  if(nr_free1 > 0)
    r1 = sum_free1/nr_free1;
  else
    r1 = (ub1+lb1)/2;
	
  if(nr_free2 > 0)
    r2 = sum_free2/nr_free2;
  else
    r2 = (ub2+lb2)/2;
	
  si->r = (r1+r2)/2;
  return (r1-r2)/2;
}
