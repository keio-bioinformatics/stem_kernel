// $Id$

#ifndef __INC_DPTABLE_H__
#define __INC_DPTABLE_H__

#include <vector>
#include <cassert>
#include "../common/cyktable.h"
#include "allocator.h"

template < class Cell, template <class> class Allocator = NewAllocator >
class DPTableTmpl : Allocator< CYKTable< Cell > >
{
private:
  typedef CYKTable< Cell > TableY;
  typedef CYKTable< TableY* > TableX;
  typedef std::vector< TableX* > Table;

  uint size_x_;
  uint size_y_;
  Table table_;

public:
  DPTableTmpl(uint n_states, uint size_x, uint size_y)
    : size_x_(size_x), size_y_(size_y), table_(n_states)
  {
    for (uint s=0; s!=table_.size(); ++s) {
      table_[s] = new TableX(size_x_);
      table_[s]->fill(NULL);
    }
  }
  
  ~DPTableTmpl()
  {
    destroy_all();
    for (uint i=0; i!=table_.size(); ++i) {
      delete table_[i];
    }
  }

  void create(uint s, uint i, uint j)
  {
    if (table_[s]->get(i,j)==NULL) {
      TableY* t = Allocator<TableY>::malloc();
      t->resize(size_y_);
      table_[s]->put(i,j,t);
    }
  }

  void create(uint s, uint i, uint j,
	      uint k_from, uint k_to, uint l_from, uint l_to)
  {
    if (table_[s]->get(i,j)==NULL) {
      TableY* t = Allocator<TableY>::malloc();
      t->resize(k_from, k_to, l_from, l_to);
      table_[s]->put(i,j,t);
    }
  }

  void create(uint s, uint i, uint j, const Cell& val)
  {
    create(s, i,j);
    table_[s]->get(i,j)->fill(val);
  }

  void create(uint s, uint i, uint j,
	      uint k_from, uint k_to, uint l_from, uint l_to,
	      const Cell& val)
  {
    create(s, i,j, k_from, k_to, l_from, l_to);
    table_[s]->get(i,j)->fill(val);
  }

  void create_all()
  {
    for (uint s=0; s!=table_.size(); ++s) {
      for (uint j=0; j!=size_x_; ++j) {
	for (uint i=j; ; --i) {
	  create(s, i,j);
	  if (i==0) break;
	}
      }
    }
  }

  void create_all(const Cell& val)
  {
    for (uint s=0; s!=table_.size(); ++s) {
      for (uint j=0; j!=size_x_; ++j) {
	for (uint i=j; ; --i) {
	  create(s, i,j, val);
	  if (i==0) break;
	}
      }
    }
  }

  void destroy(uint s, uint i, uint j)
  {
    if (table_[s]->get(i,j)!=NULL) {
      Allocator<TableY>::free(table_[s]->get(i,j));
      table_[s]->put(i,j,NULL);
    }
  }

  void destroy_all()
  {
    for (uint s=0; s!=table_.size(); ++s) {
      for (uint j=0; j!=size_x_; ++j) {
	for (uint i=j; ; --i) {
	  destroy(s,i,j);
	  if (i==0) break;
	}
      }
    }
  }

  const Cell& operator()(uint s, uint i, uint j, uint k, uint l) const
  {
    assert(table_[s]);
    assert(table_[s]->get(i,j));
    return table_[s]->get(i,j)->get(k,l);
  }

  Cell& operator()(uint s, uint i, uint j, uint k, uint l)
  {
    assert(table_[s]);
    assert(table_[s]->get(i,j));
    return table_[s]->get(i,j)->get(k,l);
  }
};

#endif

// Local Variables:
// mode: C++
// End:
