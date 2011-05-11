// $Id$

#ifndef __INC_GLOB_WRAPPER_H__
#define __INC_GLOB_WRAPPER_H__

#include <glob.h>
#include <list>

typedef unsigned int uint;

class Glob
{
public:
  typedef std::list<std::string>::const_iterator const_iterator;

  Glob() : path_() { }
  
  Glob(const char* pattern, uint flag=0) : path_()
  {
    expand(pattern, flag);
  }

  uint expand(const char* pattern, uint flag=0)
  {
    glob_t buf;
    glob(pattern, flag, NULL, &buf);
    for (uint i=0; i!=buf.gl_pathc; ++i) {
      path_.push_back(std::string(buf.gl_pathv[i]));
    }
    globfree(&buf);
    return path_.size();
  }
  const_iterator begin() const { return path_.begin(); }
  const_iterator end() const { return path_.end(); }
  bool empty() const { return path_.empty(); }
  uint size() const { return path_.size(); }
  
private:
  std::list<std::string> path_;
};

#endif	//  __INC_GLOB_WRAPPER_H__
