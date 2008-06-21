/*
 * $Id$
 *
 * PHMMTS -- an implementation of Pair Hidden Markov Models on Tree Structures,
 *           which is based on "Pair Hidden Markov Models on Tree Structures",
 *           Yasubumi Sakakibara, ISMB 2003.
 *
 * Copyright (C) 2003-2005 Kengo Sato
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#ifndef __INC_CYKTABLE_H__
#define __INC_CYKTABLE_H__

#include <cassert>
#include <utility>

typedef unsigned int uint;
typedef std::pair<uint,uint> Pos;

/**
 * @brief A templete class of CYK tables
 *        used for Cocke-Younger-Kasami (CYK) algorithm
 */
template <class T>
class CYKTable
{
public:
  typedef T value_type;

private:
  T* table_;
  T** ptr_;
  uint i_from_;
  uint i_to_;
  uint j_from_;
  uint j_to_;
  uint allocated_size_;
  
public:
  CYKTable()
    : table_(NULL), ptr_(NULL), j_to_(0), allocated_size_(0)
  {
  }

  CYKTable(uint i_from, uint i_to, uint j_from, uint j_to)
    : table_(NULL), ptr_(NULL), j_to_(0), allocated_size_(0)
  {
    setup(i_from, i_to, j_from, j_to);
  }

  CYKTable(uint size)
    : table_(NULL), ptr_(NULL), j_to_(0), allocated_size_(0)
  {
    setup(0, size, 0, size);
  }
    
  ~CYKTable() { delete[] table_; delete[] ptr_; }

  uint estimate_size(uint i_from, uint i_to, uint j_from, uint j_to)
  {
    uint ret=0;
    for (uint j=j_from; j!=j_to; ++j) {
      if (j==i_from) {
	ret++;
	continue;
      }
#if 0
      for (uint i=std::min(j,i_to-1); ; --i) {
	ret++;
	if (i==i_from) break;
      }
#else
      ret += std::min(j,i_to-1)-i_from+1;
#endif
    }
    return ret;
  }

  void setup(uint i_from, uint i_to, uint j_from, uint j_to)
  {
    assert(i_from<=i_to);
    assert(j_from<=j_to);
    assert(i_from<=j_from);
    assert(i_to<=j_to);

    uint new_size = estimate_size(i_from, i_to, j_from, j_to);
    if (allocated_size_!=new_size) {
      delete[] table_;
      allocated_size_ = new_size;
      table_ = new T[allocated_size_];
    }
    if (j_to_!=j_to) {
      delete[] ptr_;
      ptr_ = new T*[j_to];
    }
    std::fill(ptr_, ptr_+j_to, static_cast<T*>(NULL));

    uint x=0;
    for (uint j=j_from; j!=j_to; ++j) {
      ptr_[j] = &table_[x - i_from];
      if (j==i_from) {
	x++;
	continue;
      }
      for (uint i=std::min(j,i_to-1); ; --i) {
	x++;
	if (i==i_from) break;
      }
    }
    assert(x==allocated_size_);

    i_from_ = i_from;
    i_to_ = i_to;
    j_from_ = j_from;
    j_to_ = j_to;
  }

  void resize(uint i_from, uint i_to, uint j_from, uint j_to)
  {
    setup(i_from, i_to, j_from, j_to);
  }

  void resize(uint size)
  {
    resize(0, size, 0, size);
  }

  void fill(const T& val)
  {
    std::fill(table_, table_+allocated_size_, val);
  }

  inline
  void put(uint i, uint j, const T& val)
  {
    assert(i>=i_from_);
    assert(i<i_to_);
    assert(j>=j_from_);
    assert(j<j_to_);
    assert(i<=j);
    assert(ptr_[j]!=NULL);
    ptr_[j][i]=val;
  }

  inline
  void put(const Pos& pos, const T& val)
  {
    put(pos.first, pos.second, val);
  }

  inline
  const T& get(uint i, uint j) const
  {
    assert(i>=i_from_);
    assert(i<=i_to_);
    assert(j>=j_from_);
    assert(j<=j_to_);
    assert(i<=j);
    assert(ptr_[j]!=NULL);
    return ptr_[j][i];
  }

  inline
  const T& get(const Pos& pos) const
  {
    return get(pos.first, pos.second);
  }

  inline
  T& get(uint i, uint j)
  {
    assert(i>=i_from_);
    assert(i<=i_to_);
    assert(j>=j_from_);
    assert(j<=j_to_);
    assert(i<=j);
    assert(ptr_[j]!=NULL);
    return ptr_[j][i];
  }

  inline
  T& get(const Pos& pos)
  {
    return get(pos.first, pos.second);
  }


  inline
  const T& operator()(uint i, uint j) const
  {
    return get(i,j);
  }

  inline
  const T& operator()(const Pos& pos) const
  {
    return get(pos);    
  }

  inline
  T& operator()(uint i, uint j)
  {
    return get(i,j);
  }

  inline
  T& operator()(const Pos& pos)
  {
    return get(pos);    
  }
};

#endif // __INC_CYKTABLE_H__

// Local Variables:
// mode: C++
// End:
