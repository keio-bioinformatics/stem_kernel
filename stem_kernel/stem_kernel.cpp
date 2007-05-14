// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "stem_kernel.h"

#include <iostream>
#include <list>
#include <boost/multi_array.hpp>
#include "phmm.h"

template < class ValueType, class BPMat >
void
StemKernel<ValueType, BPMat>::
alignment_constraints(const std::string& x, const std::string& y,
		      std::vector<uint>& c_low,
		      std::vector<uint>& c_high) const
{
  c_low.resize(x.size()+1);
  c_high.resize(x.size()+1);
  std::fill(c_low.begin(), c_low.end(), 0);
  std::fill(c_high.begin(), c_high.end(), y.size());

  if (ali_bound_>0.0) {
    PairHMM<> phmm;
    PairHMM<>::FBTable tb;
    PairHMM<>::DPPath path;
    phmm.forward_backward(x, y, tb);
    phmm.map_path(tb, path);
#if 0
    for (PairHMM<>::DPPath::const_iterator p=path.begin(); p!=path.end(); ++p) {
      if (p->s==PairHMM<>::M && tb[p->s][p->x][p->y]>=ali_bound_) {
	std::cout << "(" << p->s << ","
		  << p->x << "," << p->y << ") "
		  << tb[p->s][p->x][p->y] << std::endl;
      }
    }
#endif
    uint low_x=0, low_y=0;
    PairHMM<>::DPPath::const_iterator p;
    for (p=path.begin(); p!=path.end(); ++p) {
      if (p->s==PairHMM<>::M && tb[p->s][p->x][p->y]>=ali_bound_) {
	for (uint i=low_x; i!=p->x; ++i) {
	  c_low[i] = low_y;
	  c_high[i] = p->y;
	}
	c_low[p->x]=p->y;
	c_high[p->x]=p->y;
	low_x = p->x+1;
	low_y = p->y;
      }
    }
    for (uint i=low_x; i!=x.size()+1; ++i) {
      c_low[i] = low_y;
      c_high[i] = y.size();
    }
    
    if (band_>0) {
      for (uint i=0; i!=x.size()+1; ++i) {
	if (c_high[i]-c_low[i]<band_*2) {
	  uint j = (c_high[i]+c_low[i])/2;
	  c_low[i] = j<band_ ? 0 : j-band_;
	  c_high[i] = j+band_>y.size() ? y.size() : j+band_;
	}
      }
    }
  } else if (band_>0) {
    for (uint i=0; i!=x.size()+1; ++i) {
      uint j=static_cast<uint>(static_cast<double>(i)/x.size()*y.size()+0.5);
      c_low[i] = j<band_ ? 0 : j-band_;
      c_high[i] = j+band_>y.size() ? y.size() : j+band_;
    }
  }

#if 0
  for (uint i=0; i!=x.size()+1; ++i) {
    std::cout << i << ": " << c_low[i] << " " << c_high[i] << std::endl;
  }
#endif
}

enum { K0=0, K1=1, K2=2, K3=3, G0=4, G1=5, G2=6, G3=7 };

template <class DPTable, class ValueType>
static inline
void
dp_init(DPTable& dp, uint i, uint j, uint k, uint l, ValueType g)
{
  dp(K0,i,j,k,l) = dp(K0,i,j-1,k,l);
  dp(G0,i,j,k,l) = dp(G0,i,j-1,k,l)*g;
  dp(K1,i,j,k,l) = dp(K1,i+1,j,k,l);
  dp(G1,i,j,k,l) = dp(G1,i+1,j,k,l)*g;
  dp(K2,i,j,k,l) = dp(K2,i,j,k,l-1);
  dp(G2,i,j,k,l) = dp(G2,i,j,k,l-1)*g;
  dp(K3,i,j,k,l) = dp(K3,i,j,k+1,l);
  dp(G3,i,j,k,l) = dp(G3,i,j,k+1,l)*g;
}

template <class DPTable, class ValueType>
static inline
void
dp_update(DPTable& dp, uint i, uint j, uint k, uint l, ValueType g)
{
  dp(K2,i,j,k,l) += dp(K3,i,j,k,l);
  dp(G2,i,j,k,l) += dp(G3,i,j,k,l);
  dp(K1,i,j,k,l) += dp(K2,i,j,k,l);
  dp(G1,i,j,k,l) += dp(G2,i,j,k,l);
  dp(K0,i,j,k,l) += dp(K1,i,j,k,l);
  dp(G0,i,j,k,l) += dp(G1,i,j,k,l);
}

template < class ValueType, class BPMat >
ValueType
StemKernel<ValueType, BPMat>::
partial_dp(const std::string& x, const std::string& y) const
{
  const value_type& g = gap_;
  std::vector<value_type> g_pow(y.size()+1);
  g_pow[0]=1.0;
  for (uint i=1; i!=g_pow.size(); ++i) g_pow[i]=g_pow[i-1]*g;

  DPTable dp(8, x.size()+1, y.size()+1);

  std::vector<uint> c_low, c_high;
  alignment_constraints(x, y, c_low, c_high);

  BPMat bp_x(x,loop_,use_GU_);
  BPMat bp_y(y,loop_,use_GU_);

  for (uint j=0; j!=x.size()+1; ++j) {
    dp.create(K0,j,j, 1.0); dp.create(G0,j,j);
    dp.create(K1,j,j, 0.0); dp.create(G1,j,j, 0.0);
    dp.create(K2,j,j, 0.0); dp.create(G2,j,j, 0.0);
    dp.create(K3,j,j, 0.0); dp.create(G3,j,j, 0.0);

    for (uint l=0; l!=y.size()+1; ++l) {
      dp(G0,j,j,l,l) = 1.0;
      if (l==0) continue;
      for (uint k=l-1; ; --k) {
	dp(G0,j,j,k,l)=dp(G0,j,j,k+1,l)*g;
	if (k==0) break;
      }
    }
    if (j==0) continue;

    for (uint i=j-1; ; --i) {
      float bp_ij = bp_x.prob(i,j-1);
      dp.create(K0,i,j, 0.0); dp.create(G0,i,j, 0.0);
      dp.create(K1,i,j, 0.0); dp.create(G1,i,j, 0.0);
      dp.create(K2,i,j, 0.0); dp.create(G2,i,j, 0.0);
      dp.create(K3,i,j, 0.0); dp.create(G3,i,j, 0.0);
#if 0
      for (uint l=0; l!=c_low[j]; ++l) {
	dp(K0,i,j,l,l) = 1.0;
	dp(G0,i,j,l,l) = dp(G0,i+1,j,l,l)*g;
	if (l==0) continue;
	
	for (uint k=l-1; ; --k) {
	  dp_init(dp,i,j,k,l,g); dp_update(dp,i,j,k,l,g);
	  if (k==0) break;
	}
      }
#endif
      for (uint l=c_low[j]; l!=c_high[j]+1; ++l) {
	dp(K0,i,j,l,l) = 1.0;
	dp(G0,i,j,l,l) = dp(G0,i+1,j,l,l)*g;
	if (l==0) continue;
#if 0
	for (uint k=l-1; k>c_high[i]; --k) {
	  dp_init(dp,i,j,k,l,g); dp_update(dp,i,j,k,l,g);
	  //if (k<=c_low[i]) continue;
        }
#endif
	for (uint k=std::min(l-1,c_high[i]); k>=c_low[i]; --k) {
	  if (l<=c_high[j-1]) {
	    dp(K0,i,j,k,l) = dp(K0,i,j-1,k,l);
	    dp(G0,i,j,k,l) = dp(G0,i,j-1,k,l)*g;
	  } else {
	    // approximation
	    dp(K0,i,j,k,l) = dp(K0,i,j-1,k,c_high[j-1]);
	    dp(G0,i,j,k,l) = dp(G0,i,j-1,k,c_high[j-1])*g*g;
	  }

	  if (k>=c_low[i+1]) {
	    dp(K1,i,j,k,l) = dp(K1,i+1,j,k,l);
	    dp(G1,i,j,k,l) = dp(G1,i+1,j,k,l)*g;
	  } else {
#if 1				// approximation
	    dp(K1,i,j,k,l) = dp(K1,i+1,j,c_low[i+1],l);
	    dp(G1,i,j,k,l) = dp(G1,i+1,j,c_low[i+1],l)*g*g;
#else
	    value_type k1=0.0, g1=0.0;
	    for (uint ii=j; ; --ii) {
	      value_type k2=dp(K2,ii,j,c_low[ii],l);
	      value_type g2=dp(G2,ii,j,c_low[ii],l)*g_pow[c_low[ii]-k];
	      for (uint ll=k; ll!=c_low[ii]; ++ll) {
		k2 += dp(K3,ii,j,ll,ll);
		g2 += dp(G3,ii,j,ll,ll)*g_pow[c_low[ii]-k];
	      }
	      k1 = k1   + k2;
	      g1 = g1*g + g2;
	      if (ii==i+1) break;
	    }
	    dp(K1,i,j,k,l) = k1;
	    dp(G1,i,j,k,l) = g1*g;
#endif
	  }

	  if (l-1>=c_low[j] || k==l-1) {
	    dp(K2,i,j,k,l) = dp(K2,i,j,k,l-1);
	    dp(G2,i,j,k,l) = dp(G2,i,j,k,l-1)*g;
	  } else {
	    dp(K2,i,j,k,l) = 0.0;
	    dp(G2,i,j,k,l) = 0.0;
	    for (uint ll=k; ll!=l; ++ll) {
	      dp(K2,i,j,k,l) += dp(K3,i,j,ll,ll);
	      dp(G2,i,j,k,l) += dp(G3,i,j,ll,ll)*g_pow[l-k];
	    }
	  }

	  if (k+1<=c_high[i]) {
	    dp(K3,i,j,k,l) = dp(K3,i,j,k+1,l);
	    dp(G3,i,j,k,l) = dp(G3,i,j,k+1,l)*g;
	  } else {
	    dp(K3,i,j,k,l) = dp(K3,i,j,l,l);
	    dp(G3,i,j,k,l) = dp(G3,i,j,l,l)*g_pow[l-k];
	  }
	  
	  if (bp_ij>bp_bound_) {
	    float bp_kl=bp_y.prob(k,l-1);
	    if (bp_kl>bp_bound_) {
	      if (x[i]==y[k] && x[j-1]==y[l-1]) {
		dp(K3,i,j,k,l) +=
		  dp(G0,i+1,j-1,k+1,l-1) * stack_ * bp_ij * bp_kl;
		dp(G3,i,j,k,l) +=
		  dp(G0,i+1,j-1,k+1,l-1);
	      } else {
		dp(K3,i,j,k,l) +=
		  dp(G0,i+1,j-1,k+1,l-1) * stack_ * subst_ * bp_ij * bp_kl;
	      }
	    }
	  }
	  dp_update(dp,i,j,k,l,g);
	  if (k==c_low[i]) break;
	}

	if (c_low[i]==0) continue;
#if 0
	for (uint k=c_low[i]-1; ; --k) {
	  dp_init(dp,i,j,k,l,g); dp_update(dp,i,j,k,l,g);
	  if (k==0) break;
	}
#endif
      }
#if 0
      for (uint l=c_high[j]+1; l!=y.size()+1; ++l) {
	dp(K0,i,j,l,l) = 1.0;
	dp(G0,i,j,l,l) = dp(G0,i+1,j,l,l)*g;
	if (l==0) continue;

	for (uint k=l-1; ; --k) {
	  dp_init(dp,i,j,k,l,g); dp_update(dp,i,j,k,l,g);
	  if (k==0) break;
        }
      }
#endif
      if (i==0) break;
    }

    for (uint i=0; i<j; ++i) {
      dp.destroy(K0,i,j-1); dp.destroy(G0,i,j-1);
      dp.destroy(K1,i,j-1); dp.destroy(G1,i,j-1);
      dp.destroy(K2,i,j-1); dp.destroy(G2,i,j-1);
      dp.destroy(K3,i,j-1); dp.destroy(G3,i,j-1);
    }
  }

  return dp(K0, 0, x.size(), 0, y.size());
}

template < class ValueType, class BPMat >
ValueType
StemKernel<ValueType, BPMat>::
full_dp(const std::string& x, const std::string& y) const
{
  const value_type& g = gap_;
  DPTable dp(8, x.size()+1, y.size()+1);

  BPMat bp_x(x,loop_,use_GU_);
  BPMat bp_y(y,loop_,use_GU_);

  for (uint j=0; j!=x.size()+1; ++j) {
    dp.create(K0,j,j, 1.0); dp.create(G0,j,j);
    dp.create(K1,j,j, 0.0); dp.create(G1,j,j, 0.0);
    dp.create(K2,j,j, 0.0); dp.create(G2,j,j, 0.0);
    dp.create(K3,j,j, 0.0); dp.create(G3,j,j, 0.0);

    for (uint l=0; l!=y.size()+1; ++l) {
      dp(G0,j,j,l,l) = 1.0;
      if (l==0) continue;
      for (uint k=l-1; ; --k) {
	dp(G0,j,j,k,l)=dp(G0,j,j,k+1,l)*g;
	if (k==0) break;
      }
    }
    if (j==0) continue;

    for (uint i=j-1; ; --i) {
      float bp_ij = bp_x.prob(i,j-1);
      dp.create(K0,i,j, 0.0); dp.create(G0,i,j, 0.0);
      dp.create(K1,i,j, 0.0); dp.create(G1,i,j, 0.0);
      dp.create(K2,i,j, 0.0); dp.create(G2,i,j, 0.0);
      dp.create(K3,i,j, 0.0); dp.create(G3,i,j, 0.0);
      for (uint l=0; l!=y.size()+1; ++l) {
	dp(K0,i,j,l,l) = 1.0;
	dp(G0,i,j,l,l) = dp(G0,i+1,j,l,l)*g;
	if (l==0) continue;
	for (uint k=l-1; ; --k) {
	  dp_init(dp,i,j,k,l,g);
 	  if (bp_ij>bp_bound_) {
	    float bp_kl=bp_y.prob(k,l-1);
	    if (bp_kl>bp_bound_) {
	      if (x[i]==y[k] && x[j-1]==y[l-1]) {
		dp(K3,i,j,k,l) +=
		  dp(G0,i+1,j-1,k+1,l-1) * stack_ * bp_ij * bp_kl;
		dp(G3,i,j,k,l) +=
		  dp(G0,i+1,j-1,k+1,l-1);
	      } else {
		dp(K3,i,j,k,l) +=
		  dp(G0,i+1,j-1,k+1,l-1) * stack_ * subst_ * bp_ij * bp_kl;
	      }
	    }
	  }
	  dp_update(dp,i,j,k,l,g);
	  if (k==0) break;
	}
      }
      if (i==0) break;
    }

    for (uint i=0; i<j; ++i) {
      dp.destroy(K0,i,j-1); dp.destroy(G0,i,j-1);
      dp.destroy(K1,i,j-1); dp.destroy(G1,i,j-1);
      dp.destroy(K2,i,j-1); dp.destroy(G2,i,j-1);
      dp.destroy(K3,i,j-1); dp.destroy(G3,i,j-1);
    }
  }

  return dp(K0, 0, x.size(), 0, y.size());
}

class NormalBasePair
{
public:
  NormalBasePair(const std::string& seq, uint loop, bool useGU)
    : seq_(seq), loop_(loop)
  {
    assert(useGU==false);
  }

  float prob(uint i, uint j)
  {
    const char& a=seq_[i];
    const char& b=seq_[j];
    return i+1+loop_<=j && 
      (a=='a' && b=='u' || a=='u' && b=='a' ||
       a=='g' && b=='c' || a=='c' && b=='g') ? 1.0 : 0.0;
  }

private:
  const std::string& seq_;
  const uint loop_;
};

class WobbleBasePair
{
public:
  WobbleBasePair(const std::string& seq, uint loop, bool useGU)
    : seq_(seq), loop_(loop)
  {
    assert(useGU==true);
  }

  float prob(uint i, uint j)
  {
    const char& a=seq_[i];
    const char& b=seq_[j];
    return i+1+loop_<=j &&
      (a=='a' && b=='u' || a=='u' && b=='a' ||
       a=='g' && b=='c' || a=='c' && b=='g' ||
       a=='g' && b=='u' || a=='u' && b=='g') ? 1.0 : 0.0;
  }

private:
  const std::string& seq_;
  const uint loop_;
};

#ifdef HAVE_LIBRNA
class BPMatrix
{
public:
  BPMatrix(const std::string& seq, uint loop, bool useGU)
    : bp_(seq, useGU), seq_(seq), loop_(loop), useGU_(useGU)
  {
    bp_.fold();
  }

  float prob(uint i, uint j)
  {
    return bp_(i+1, j+1);
  }

private:
  PFWrapper bp_;
  const std::string& seq_;
  const uint loop_;
  bool useGU_;
};
#endif

// instantiation
template class StemKernel<double,NormalBasePair>;
template class StemKernel<double,WobbleBasePair>;
#ifdef HAVE_LIBRNA
template class StemKernel<double,BPMatrix>;
#endif
