// $Id$

#ifndef __INC_STEM_KERNEL_H__
#define __INC_STEM_KERNEL_H__

template < class ST, class D >
class StemKernel
{
public:
  typedef ST ScoreTable;
  typedef typename ScoreTable::value_type value_type;
  typedef D Data;

public:
  StemKernel(const ScoreTable& st, uint len_band=0)
    : st_(st), len_band_(len_band) { }

  value_type operator()(const Data& xx, const Data& yy) const;

  template < class Seq >
  value_type operator()(const Seq& x, const Seq& y) const
  {
    return (*this)(Data(x), Data(y));
  }

private:
  const ScoreTable& st_;
  uint len_band_;
};

#endif // __INC_STEM_KERNEL_H__

// Local Variables:
// mode: C++
// End:
