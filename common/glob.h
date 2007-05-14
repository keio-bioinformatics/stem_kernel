// $Id:$

#ifndef __INC_GLOB_H__
#define __INC_GLOB_H__

#include <glob.h>
#include <list>
#include <string>

class Glob
{
public:
  typedef std::list<std::string>::const_iterator const_iterator;
  
  Glob(const char* pattern, uint flag=0) : path_()
  {
    glob_t buf;
    glob(pattern, flag, NULL, &buf);
    for (uint i=0; i!=buf.gl_pathc; ++i) {
      path_.push_back(std::string(buf.gl_pathv[i]));
    }
    globfree(&buf);
  }
  const_iterator begin() const { return path_.begin(); }
  const_iterator end() const { return path_.end(); }
  bool empty() const { return path_.empty(); }
  
private:
  std::list<std::string> path_;
};

#endif	// __INC_GLOB_H__

// Local Variables:
// mode: C++
// End:
