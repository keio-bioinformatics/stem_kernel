// $Id$
#ifndef __INC_DPTABLE_H__
#define __INC_DPTABLE_H__

#include <vector>
#include <list>

typedef unsigned int uint;

template < class V >
class DPTable
{
public:
  typedef V value_type;
  
  DPTable(uint sz_x, uint sz_y)
    : sz_y_(sz_y), tbl_(sz_x, NULL), pool_()
  {
  }

  ~DPTable()
  {
    for (uint i=0; i!=tbl_.size(); ++i) {
      if (tbl_[i]) delete tbl_[i];
    }
    typename std::list< std::vector<value_type>* >::iterator x;
    for (x=pool_.begin(); x!=pool_.end(); ++x)
      delete *x;
  }

  void allocate(uint i)
  {
    if (pool_.empty()) {
      tbl_[i] = new std::vector<value_type>(sz_y_/*, 0.0*/);
    } else {
      tbl_[i] = pool_.front();
      //std::fill(tbl_[i]->begin(), tbl_[i]->end(), 0.0);
      pool_.pop_front();
    }
  }

  void deallocate(uint i)
  {
    if (tbl_[i]) {
      pool_.push_back(tbl_[i]);
      tbl_[i] = NULL;
    }
  }

  value_type operator()(uint x, uint y) const
  {
    return (*tbl_[x])[y];
  }

  value_type& operator()(uint x, uint y)
  {
    return (*tbl_[x])[y];
  }

  const std::vector<value_type>& operator[](uint x) const
  {
    return *tbl_[x];
  }

  std::vector<value_type>& operator[](uint x)
  {
    return *tbl_[x];
  }

private:
  uint sz_y_;
  std::vector< std::vector<value_type>* > tbl_;
  std::list< std::vector<value_type>* > pool_;
};

#endif	//  __INC_DPTABLE_H__
