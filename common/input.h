// $Id$

#ifndef __INC_INPUT_H__
#define __INC_INPUT_H__

#include <iosfwd>
#include <string>
#include <list>
#include <exception>

class InputSeq
{
public:
  typedef std::pair<std::string, double> StrInfo;

  class ParseError : public std::exception
  {
  public:
    ParseError() { }
    virtual ~ParseError() throw() { }
    virtual const char* what() const throw();
  };

public:
  InputSeq(std::istream& in);
  ~InputSeq();

  bool get_next_seq();
  const std::string& seq() const { return seq_; }
  const std::string& name() const { return name_; }
  const std::list<StrInfo>& structs() const { return str_; }

private:
  std::istream& in_;
  std::string seq_;
  std::string name_;
  std::list<StrInfo> str_;
};


#endif	//  __INC_INPUT_H__

// Local Variables:
// mode: C++
// End:
