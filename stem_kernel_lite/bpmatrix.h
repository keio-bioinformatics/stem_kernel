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
  typedef std::list< boost::shared_ptr<BPMatrix> >::const_iterator matrix_iterator;
  typedef std::list< std::vector<uint> >::const_iterator idx_map_iterator;
  
  
public:
  BPMatrix(const std::string& s, const Options& opts);

  BPMatrix(const std::list<std::string>& ma, const Options& opts);

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

  void add_matrix(const boost::shared_ptr<BPMatrix>& bp,
		  const std::vector<uint>& idx_map)
  {
    bps_.push_back(bp);
    idx_map_.push_back(idx_map);
  }

  matrix_iterator matrix_begin() const { return bps_.begin(); }
  matrix_iterator matrix_end() const { return bps_.end(); }
  idx_map_iterator idx_map_begin() const { return idx_map_.begin(); }
  idx_map_iterator idx_map_end() const { return idx_map_.end(); }

  uint n_matrices() const { return bps_.size(); }

private:
  uint sz_;
  CYKTable<double> table_;

  std::list< boost::shared_ptr<BPMatrix> > bps_;
  std::list< std::vector<uint> > idx_map_;
};

#endif	// __INC_BPMATRIX_H__

// Local Variables:
// mode: C++
// End:
