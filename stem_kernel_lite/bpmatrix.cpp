// $Id: bpmatrix.cpp 38 2006-10-24 10:39:00Z satoken $

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "bpmatrix.h"
#include "../common/rna.h"
#include <cmath>
#include <algorithm>

BPMatrix::
BPMatrix(uint sz)
#ifdef HAVE_BOOST_RANDOM
  : sz_(sz),
    table_(sz_+1),
    rand_(RNG(time(NULL)), Dist(0,1))
#else
  : sz_(sz),
    table_(sz_+1)
#endif
{
  table_.fill(0.0);
}

template < class V >
BPMatrix::
BPMatrix(uint sz, const V* pr, const int* iindx)
#ifdef HAVE_BOOST_RANDOM
  : sz_(sz),
    table_(sz_+1),
    rand_(RNG(time(NULL)), Dist(0,1))
#else
  : sz_(sz),
    table_(sz_+1)
#endif
{
  table_.fill(0.0);
  for (uint j=2; j!=sz_+1; ++j) {
    for (uint i=j-1; ; --i) {
      table_(i,j) = pr[iindx[i]-j];
      if (i==1) break;
    }
  }
}

template <class Seq>
BPMatrix::
BPMatrix(const std::list<Seq>& ali,
	 const std::list< boost::shared_ptr<BPMatrix> >& bps)
#ifdef HAVE_BOOST_RANDOM
  : sz_(ali.begin()->size()),
    table_(sz_+1),
    rand_(RNG(time(NULL)), Dist(0,1))
#else
  : sz_(ali.begin()->size()),
    table_(sz_+1)
#endif
{
  assert(ali.size()==bps.size());
  typedef typename Seq::value_type rna_type;
  const rna_type GAP = RNASymbol<rna_type>::GAP;
  uint n_seq=ali.size();
  table_.fill(0.0);
  typename std::list<Seq>::const_iterator a;
  std::list< boost::shared_ptr<BPMatrix> >::const_iterator b;
  for (a=ali.begin(), b=bps.begin(); a!=ali.end() && b!=bps.end(); ++a, ++b) {
    std::vector<uint> idxmap(a->size());
    for (uint i=0, j=0; i!=a->size(); ++i) {
      if ((*a)[i]!=GAP) idxmap[j++]=i;
    }
    for (uint j=1; j!=(*b)->size(); ++j) {
      for (uint i=j-1; /*i>=0*/; --i) {
	(*this)(idxmap[i]+1,idxmap[j]+1) += (**b)(i+1,j+1);
	if (i==0) break;
      }
    }
  }
  for (uint j=1; j!=sz_; ++j) {
    for (uint i=j-1; /*i>=0*/; --i) {
      (*this)(i+1,j+1) = (*this)(i+1,j+1) / n_seq;
      if (i==0) break;
    }
  }
}

void
BPMatrix::
traceback(std::string& str, float bp_w /*=1.0*/) const
{
  std::vector<double> ss(sz_);
  std::fill(ss.begin(), ss.end(), 1.0);
  for (uint j=1; j!=sz_; ++j) {
    for (uint i=0; i!=j; ++i) {
      double p_ij = (*this)(i+1,j+1);
      ss[i]-=p_ij;
      ss[j]-=p_ij;
    }
  }

  CYKTable<double> bp(sz_);
  CYKTable<double> bp_m(sz_);
  bp.fill(0.0);
  bp_m.fill(0.0);
  for (uint j=1; j!=sz_; ++j) {
    for (uint i=j-1; /*i>=0*/; --i) {
      double p_ij = (*this)(i+1,j+1);
      if (i+2<=j && p_ij>0.0)
	bp(i,j)=bp(i+1,j-1)+2*bp_w*p_ij;
      if (i+1<=j && bp(i,j)<bp(i+1,j)+ss[i])
	bp(i,j)=bp(i+1,j)+ss[i];
      if (i<=j-1 && bp(i,j)<bp(i,j-1)+ss[j])
	bp(i,j)=bp(i,j-1)+ss[j];
      for (uint k=i; k<j; ++k) {
	if (bp_m(i,j)<bp(i,k)+bp(k+1,j))
	  bp_m(i,j)=bp(i,k)+bp(k+1,j);
      }
      if (bp(i,j)<bp_m(i,j))
	bp(i,j)=bp_m(i,j);
      if (i==0) break;
    }
  }

  str.resize(sz_);
  std::fill(str.begin(), str.end(), '.');
  //std::cout << bp(0,sz_-1) << std::endl;
  traceback(bp, bp_m, ss, str, 0, sz_-1, bp_w);
}


void
BPMatrix::
traceback(const CYKTable<double>& bp,
	  const CYKTable<double>& bp_m,
	  const std::vector<double>& ss,
	  std::string& str,
	  uint i, uint j,
	  float bp_w /*=1.0*/) const
{
  if (i==j) return;
  double p=0.0, l=0.0, r=0.0, m=0.0;

  double p_ij = (*this)(i+1,j+1);
  if (i+2<=j && p_ij>0.0) p=exp(bp(i+1,j-1)+2*bp_w*p_ij);
  if (i+1<=j) l=exp(ss[i]+bp(i+1,j));
  if (i<=j-1) r=exp(ss[j]+bp(i,j-1));
  m=exp(bp_m(i,j));

  double rnd=rand01()*(p+l+r+m);
  if (/*i+2<=j && p_ij>0.0 &&*/ rnd<p) {
    str[i]='('; 
    str[j]=')';
    //std::cout << "p: " << i << "," << j << std::endl;
    traceback(bp, bp_m, ss, str, i+1, j-1);
    return;
  }
  if (/*i+1<=j &&*/ rnd<p+l) {
    str[i]='.'; 
    //std::cout << "l: " << i << "," << j << std::endl;
    traceback(bp, bp_m, ss, str, i+1, j);
    return;
  }
  if (/*i<=j-1 &&*/ rnd<p+l+r) {
    str[j]='.'; 
    //std::cout << "r: " << i << "," << j << std::endl;
    traceback(bp, bp_m, ss, str, i, j-1);
    return;
  }
  {
    std::vector<double> br(j-i);
    br[0]=exp(bp(i,i)+bp(i+1,j));
    for (uint k=i+1; k<j; ++k)
      br[k-i]=br[k-i-1]+exp(bp(i,k)+bp(k+1,j));
    double rnd_m=rand01()*br[j-i-1];
    for (uint k=i; k<j; ++k) {
      if (rnd_m<br[k-i]) {
	//std::cout << "b: " << i << "," << j << "," << k << std::endl;
	traceback(bp, bp_m, ss, str, i, k);
	traceback(bp, bp_m, ss, str, k+1, j);
	return;
      }
    }
  }
}

#ifndef HAVE_BOOST_RANDOM
double
BPMatrix::
rand01() const 
{
  return rand()/(static_cast<double>(RAND_MAX)+1);
}
#else
double
BPMatrix::
rand01() const 
{
  return rand_();
}
#endif

// template instantiation

template 
BPMatrix::
BPMatrix(const std::list<std::string>& ali,
	 const std::list< boost::shared_ptr<BPMatrix> >& bps);

template 
BPMatrix::
BPMatrix(const std::list<RNASequence>& ali,
	 const std::list< boost::shared_ptr<BPMatrix> >& bps);

extern "C" {
#include <ViennaRNA/fold_vars.h>
};

template
BPMatrix::
BPMatrix(uint sz, const FLT_OR_DBL* pr, const int* iindx);
