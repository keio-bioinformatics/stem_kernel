// $Id: bpmatrix.h 35 2006-10-17 10:53:34Z satoken $

#ifndef __INC_BPMATRIX_H__
#define __INC_BPMATRIX_H__

#include <string>
#include "../common/cyktable.h"
#include "../common/rna.h"
#include <iostream>
#include <list>
#include <boost/program_options.hpp>

class BPMatrix {
public:
  enum { FOLD, LFOLD, SFOLD, ALIFOLD }; // available methods
  struct Options
  {
    bool alifold;
    float th;
    uint n_samples;

    void add_options(boost::program_options::options_description& desc);
    uint method() const;
  };
  
public:
  BPMatrix(const std::string& s, const Options& opts);

  BPMatrix(const MASequence<std::string>& ma, const Options& opts);

  double operator()(uint i, uint j) const
  {
    return table_(i,j);
  }

  double& operator()(uint i, uint j)
  {
    return table_(i,j);
  }

  uint size() const
  {
    return sz_;
  }

private:
  uint sz_;
  CYKTable<double> table_;
};

#endif	// __INC_BPMATRIX_H__

// Local Variables:
// mode: C++
// End:
