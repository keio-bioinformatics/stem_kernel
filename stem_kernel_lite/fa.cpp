// $Id:$

#include "fa.h"
#include <iostream>
#include <sstream>
#include <string>
#include <boost/spirit/core.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/iterator/file_iterator.hpp>

using namespace boost::spirit;

template < class Seq >
struct fa_parser : public grammar< fa_parser<Seq> >
{
  fa_parser(std::string& n, std::string& s) : name(n), seq(s) { }

  std::string& name;
  std::string& seq;

  struct append_seq
  {
    std::string& seq;

    append_seq(std::string& s) : seq(s) { }

    template <class Ite>
    void operator()(Ite i1, Ite i2) const
    {
      std::string s(i1, i2);
      seq += s;
    }
  };


  template <class ScannerT>
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t fa;
    rule_t head;
    rule_t seq_l;
    rule_t seq;

    definition(const fa_parser& self)
    {
      fa = head >> seq;
      head = ch_p('>') >> (*(blank_p | graph_p))[assign_a(self.name)] >> eol_p;
      seq_l = *(graph_p - '>' - eol_p);
      seq = +(seq_l[append_seq(self.seq)] >> eol_p);
    }

    const rule_t& start() const { return fa; }
  };
};

template < class Seq >
bool
load_fa(Seq& s, file_iterator<>& fi)
{
  std::string name, seq;
  fa_parser<Seq> parser(name, seq);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) return false;
  fi = info.stop;
  char2rna(s, seq);
  return true;
}

template < class Seq >
bool
load_fa(std::list<Seq>& fa, const char* filename)
{
  uint n=0;
  file_iterator<> fi(filename);
  if (!fi) {
    std::ostringstream os;
    os << filename << ": no such file";
    throw os.str().c_str();
    //return false;
  }
  while (1) {
    Seq s;
    if (load_fa(s, fi)) {
      n++;
      fa.push_back(s);
    } else {
      break;
    }
  }
  return n!=0;
}

#if 0
template
bool
load_fa(std::list<RNASequence>& fa, const char* filename);
#else
template
bool
load_fa(std::list<std::string>& fa, const char* filename);
#endif
