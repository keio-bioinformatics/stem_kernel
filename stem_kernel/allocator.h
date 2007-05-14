// $Id$

#ifndef __INC_ALLOCATOR_H__
#define __INC_ALLOCATOR_H__

#include <list>

template < class object_type >
class NewAllocator
{
public:
  static object_type* malloc()
  {
    return new object_type;
  }

  static void free(object_type* obj)
  {
    delete obj;
  }
};

template < class object_type >
class PoolAllocator
{
private:
  std::list<object_type*> free_list_;
  
public:
  PoolAllocator() : free_list_() { }

  ~PoolAllocator()
  {
    typename std::list<object_type*>::iterator x;
    for (x=free_list_.begin(); x!=free_list_.end(); ++x)
      delete *x;
  }

  object_type* malloc()
  {
    if (free_list_.empty()) {
      return new object_type;
    } else {
      object_type* obj = free_list_.front();
      free_list_.pop_front();
      return obj;
    }
  }

  void free(object_type* obj)
  {
    free_list_.push_back(obj);
  }
};

#endif
