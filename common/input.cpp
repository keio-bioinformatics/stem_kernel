// $Id$

#include <istream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include "input.h"

const char  HEADER = '>';
const char  TERM   = '*';
const char* BLANK  = " \t\n\r";

InputSeq::
InputSeq(std::istream& in)
  : in_(in), seq_(), name_(), str_()
{
  std::string buf;
  std::getline(in_, buf, HEADER);
}

InputSeq::
~InputSeq()
{
}

bool
InputSeq::
get_next_seq()
{
  std::string buf;
  std::stringstream ss;
  bool in_seq=true;
  name_.clear();
  seq_.clear();
  str_.clear();

  std::getline(in_, buf, HEADER);
  if (!in_) return false;
  if (*buf.rbegin()==TERM)
    buf.erase(buf.length()-1, 1);

  ss << buf;
  while (ss.get()==' ');
  ss.unget();
  std::getline(ss, name_);
  std::string line;
  while (std::getline(ss, line)) {
    if (static_cast<int>(line.find_first_of("()"))>=0) {
      in_seq=false;
      std::stringstream l(line);
      std::string p, e;
      double w=1.0;
      l >> p >> e;
      if (!e.empty()) w=exp(-atof(e.c_str()));
      str_.push_back(std::make_pair(p,w));
    } else if (in_seq) {
      seq_ += line;
    } else {
      throw ParseError();
    }
  }

  double total=0.0;
  std::list<StrInfo>::iterator x;
  for (x=str_.begin(); x!=str_.end(); ++x) total += x->second;
  for (x=str_.begin(); x!=str_.end(); ++x) x->second /= total;
  
  return true;
}

const char*
InputSeq::ParseError::
what() const throw()
{
  return "parser error";
}
