// $Id: bpmatrix.h 35 2006-10-17 10:53:34Z satoken $

#ifndef __INC_BPMATRIX_H__
#define __INC_BPMATRIX_H__

#include <vector>
#include <string>
#ifdef HAVE_BOOST_RANDOM
#include <boost/random.hpp>
#endif
#include "../common/cyktable.h"
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <list>

class BPMatrix {
public:
  BPMatrix(uint sz);

  template < class V >
  BPMatrix(uint sz, const V* pr, const int* iindx);

  template <class Seq>
  BPMatrix(const std::list<Seq>& ali,
	   const std::list< boost::shared_ptr<BPMatrix> >& bps);

  double operator()(uint i, uint j) const
  {
    return table_(i,j);
  }

  double& operator()(uint i, uint j)
  {
    return table_(i,j);
  }

  void set_value(uint i, uint j, double v)
  {
    table_(i,j)=v;
  }

  uint size() const
  {
    return sz_;
  }

  void traceback(std::string& str, float bp_w=1.0) const;

private:
  void traceback(const CYKTable<double>& bp,
		 const CYKTable<double>& bp_m,
		 const std::vector<double>& ss,
		 std::string& str,
		 uint i, uint j,
		 float bp_w=1.0) const;

  double rand01() const;

private:
  uint sz_;
  CYKTable<double> table_;
#ifdef HAVE_BOOST_RANDOM
  typedef boost::mt19937 RNG;
  typedef boost::uniform_real<> Dist;
  typedef boost::variate_generator<RNG,Dist> Gen;
  mutable Gen rand_;
#endif
};

#endif	// __INC_BPMATRIX_H__

// Local Variables:
// mode: C++
// End:
