// $Id$

#ifndef __INC_GRADIENT_H__
#define __INC_GRADIENT_H__

#include <vector>
#include <boost/multi_array.hpp>

typedef boost::multi_array<double,2> Kmat;
typedef boost::multi_array<double,3> Gmat;
typedef unsigned int uint;

class GradientComputation
{
public:
  GradientComputation(const std::vector<int>& y,
		      const std::vector<uint>& tr_i,
		      const std::vector<uint>& ts_i)
    : y_(y), tr_i_(tr_i), ts_i_(ts_i)
  {
  }

  virtual ~GradientComputation()
  {
  }
  
  double compute(const Kmat& k, const Gmat& g,
		 double C, std::vector<double>& fg, double& cg, double eps=1e-3) const;

protected:
  virtual double compute_delta(const std::vector<double>& dec_values,
			       std::vector<double>& delta) const = 0;

protected:
  const std::vector<int>& y_;
  std::vector<uint> tr_i_;
  std::vector<uint> ts_i_;
};

class GradientComputationAUC : public GradientComputation
{
public:
  GradientComputationAUC(const std::vector<int>& y,
			 const std::vector<uint>& tr_i,
			 const std::vector<uint>& ts_i)
    : GradientComputation(y, tr_i, ts_i)
  {
  }

  virtual ~GradientComputationAUC()
  {
  }
  
protected:
  double compute_delta(const std::vector<double>& dec_values,
		       std::vector<double>& delta) const;
  
};

#endif	// __INC_GRADIENT_H__

// Local Variables:
// mode: C++
// End:
