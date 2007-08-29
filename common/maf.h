// $Id:$

#ifndef __INC_MAF_H__
#define __INC_MAF_H__

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <boost/spirit/iterator/file_iterator.hpp>
#include "../common/rna.h"

template < class Seq >
bool
load_maf(MASequence<Seq>& ma, boost::spirit::file_iterator<>& fi);

template < class Seq >
bool
load_maf(std::list<Seq>& ma, boost::spirit::file_iterator<>& fi);

template < class Seq >
bool
load_maf(std::list< MASequence<Seq> >& ma, const char* filename);

class MAF
{
private:
  boost::shared_ptr<std::istream> in_;
  std::string cur_;
  
public:
  MAF() : cur_() { };

  MAF(const MAF& x) : in_(x.in_), cur_(x.cur_) { };

  MAF(const char* f)
    : in_(boost::shared_ptr<std::istream>(new std::ifstream(f))), cur_()
  {
    init();
  }

  bool open(const char* f)
  {
    in_ = boost::shared_ptr<std::ifstream>(new std::ifstream(f));
    init();
    return true;
  }

  template < class Seq >
  bool get(MASequence<Seq>& ma);

private:
  void init();
};


#endif //  __INC_MAF_H__

// Local Variables:
// mode: C++
// End:
