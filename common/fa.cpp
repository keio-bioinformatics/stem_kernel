// $Id:$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "fa.h"
#include <iostream>
#include <sstream>
#include <string>

using namespace BOOST_SPIRIT_CLASSIC_NS;

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
load_fa(MASequence<Seq>& ma, file_iterator<>& fi)
{
  file_iterator<> st = fi;
  std::string name, seq;
  fa_parser<Seq> parser(name, seq);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) {
    fi = st;
    return false;
  }
  fi = info.stop;
  Seq s;
  char2rna(s, seq);
  ma.add_seq(s);
  return true;
}

template < class Seq >
bool
load_fa(std::list<Seq>& ma, file_iterator<>& fi)
{
  file_iterator<> st = fi;
  std::string name, seq;
  fa_parser<Seq> parser(name, seq);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) {
    fi = st;
    return false;
  }
  fi = info.stop;
  Seq s;
  char2rna(s, seq);
  ma.push_back(s);
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

template < class Seq >
bool
load_fa(std::list< MASequence<Seq> >& ma, const char* filename)
{
  std::list<std::string> fa;
  bool res=load_fa(fa, filename);
  if (!res) return res;
  std::list<std::string>::const_iterator x;
  MASequence<Seq> m;
  for (x=fa.begin(); x!=fa.end(); ++x) {
    Seq r;
    char2rna(r, *x);
    m.add_seq(r);
  }
  ma.push_back(m);
  return true;
}


#if 0
template
bool
load_fa(std::list<RNASequence>& fa, const char* filename);
#else
template
bool
load_fa(std::string& seq, file_iterator<>& fi);

template
bool
load_fa(MASequence<std::string>& ma, file_iterator<>& fi);

template
bool
load_fa(std::list<std::string>& ma, file_iterator<>& fi);

template
bool
load_fa(std::list<std::string>& fa, const char* filename);

template
bool
load_fa(std::list<MASequence<std::string> >& ma, const char* filename);
#endif
