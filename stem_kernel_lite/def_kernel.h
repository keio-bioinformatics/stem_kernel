// $Id:$

#ifndef __INC_DEF_KERNEL_H__
#define __INC_DEF_KERNEL_H__

#include "score_table.h"
#include "stem_kernel.h"
#include "string_kernel.h"
#include "../common/conv_kernel.h"

template <class V, class D>
class SiStemKernel
{
public:
  typedef V value_type;
  typedef D Data;

public:
  SiStemKernel(value_type gap, value_type stack, value_type covar,
	       value_type loop_gap, uint len_band)
    : st_(gap, stack, covar, loop_gap),
      k_(st_, len_band)
  {
  }

  V operator()(const Data& x, const Data& y) const
  {
    return k_(x,y);
  }

private:
  SimpleScoreTable<D,V> st_;
  StemKernel<SimpleScoreTable<D,V>, D> k_;
};

template <class V, class D>
class SuStemKernel
{
public:
  typedef V value_type;
  typedef D Data;

public:
  SuStemKernel(value_type gap, value_type beta, value_type loop_gap,
	       uint len_band)
    : st_(gap, beta, loop_gap),
      k_(st_, len_band)
  {
  }

  V operator()(const Data& x, const Data& y) const
  {
    return k_(x,y);
  }

private:
  SubstScoreTable<D,V> st_;
  StemKernel<SubstScoreTable<D,V>, D> k_;
};

template <class V, class D>
class SuStemStrKernel
{
public:
  typedef V value_type;
  typedef D Data;

public:
  SuStemStrKernel(value_type alpha, value_type beta,
		  value_type gap, value_type loop_gap, uint len_band)
    : stem_(gap, beta, loop_gap, len_band),
      str_(loop_gap, alpha),
      k_(stem_,str_)
  {
  }

  value_type operator()(const D& x, const D& y) const
  {
    return k_(x,y);
  }

private:
  SuStemKernel<V,D> stem_;
  StringKernel<V,D> str_;
  AddKernel<SuStemKernel<V,D>, StringKernel<V,D> > k_;
};

template <class V, class D>
class LSuStemKernel
{
public:
  typedef V value_type;
  typedef D Data;
  
 public:
  LSuStemKernel(value_type gap, value_type beta, value_type loop_gap,
		uint len_band)
    : stem_(gap, beta, loop_gap, len_band),
      log_stem_(stem_),
      k_(log_stem_, beta, 0.0)
    {
    }

  value_type operator()(const D& x, const D& y) const
  {
    return k_(x,y);
  }

 private:
  SuStemKernel<V,D> stem_;
  LogKernel<SuStemKernel<V,D> > log_stem_;
  LTKernel<LogKernel<SuStemKernel<V,D> > > k_;
};

template <class V, class D>
class LSuStrKernel
{
public:
  typedef V value_type;
  typedef D Data;

 public:
  LSuStrKernel(value_type loop_gap, value_type alpha)
    : str_(loop_gap, alpha),
      log_str_(str_),
      k_(log_str_, alpha, 0.0)
    {
    }

  value_type operator()(const D& x, const D& y) const
  {
    return k_(x,y);
  }
 
 private:
  StringKernel<V,D> str_;
  LogKernel<StringKernel<V,D> > log_str_;
  LTKernel<LogKernel<StringKernel<V,D> > > k_;
};

template <class V, class D>
class LSuStemStrKernel
{
public:
  typedef V value_type;
  typedef D Data;

public:
  LSuStemStrKernel(value_type alpha, value_type beta,
		   value_type gap, value_type loop_gap, uint len_band)
    : stem_(gap, beta, loop_gap, len_band),
      str_(loop_gap, alpha),
      k_(stem_,str_)
  {
  }

  value_type operator()(const D& x, const D& y) const
  {
    return k_(x,y);
  }

private:
  LSuStemKernel<V,D> stem_;
  LSuStrKernel<V,D> str_;
  AddKernel<LSuStemKernel<V,D>, LSuStrKernel<V,D> > k_;
};

#endif
