// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "aln.h"
#include <iostream>
#include <sstream>
#include <string>
#include <deque>
#include <map>
#include <stdexcept>
#include <boost/spirit.hpp>

using namespace boost::spirit;

struct aln_parser : public grammar< aln_parser >
{
  class format_error : public std::logic_error
  {
  public:
    format_error(const std::string& msg) : std::logic_error(msg) {}
  };
  
  struct WA
  {
    WA() : cur_index(0), cur_name(), cur_seq(), names(), seqs() { }

    uint cur_index;
    std::string cur_name;
    std::string cur_seq;
    std::deque<std::string> names;
    std::deque<std::string> seqs;
  };

  WA& wa;
  aln_parser(WA& x) : wa(x) {}

  struct push_seq
  {
    WA& wa;
    push_seq(WA& x) : wa(x) { }

    template < class Ite >
    void operator()(Ite i1, Ite i2) const
    {
      assert(wa.names.size()==wa.seqs.size());
      if (wa.cur_index >= wa.names.size()) {
	wa.names.push_back(wa.cur_name);
	wa.seqs.push_back(wa.cur_seq);
      } else if (wa.names[wa.cur_index] == wa.cur_name) {
	wa.seqs[wa.cur_index] += wa.cur_seq;
      } else {
	throw format_error("format error: broken sequence name consistency");
      }
      wa.cur_index++;
    }
  };

  struct reset_index
  {
    WA& wa;
    reset_index(WA& x) : wa(x) { }

    template < class Ite >
    void operator()(Ite i1, Ite i2) const
    {
      uint l=wa.seqs[0].size();
      for (uint i=1; i!=wa.seqs.size(); ++i) {
	if (l!=wa.seqs[i].size())
	  throw format_error("format error: broken sequence length consistency");
      }
      wa.cur_index = 0;
    }
  };

  template <class ScannerT>
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t aln;
    rule_t header;
    rule_t head_word;
    rule_t body;
    rule_t body_part;
    rule_t empty;
    rule_t seq;
    rule_t status;
    
    definition(const aln_parser& self)
    {
      aln = header >> +empty >> body;
      head_word = str_p("CLUSTAL") | str_p("PROBCONS");
      header =  head_word >> +print_p >> eol_p;
      body_part = +seq[push_seq(self.wa)] >> !status >> +empty;
      body = +(body_part[reset_index(self.wa)]);
      empty = *blank_p >> eol_p;
      seq
	= (+graph_p - head_word)[assign_a(self.wa.cur_name)]
	>> +blank_p >> (+graph_p)[assign_a(self.wa.cur_seq)]
	>> *blank_p >> eol_p;
      status = *(chset<>("*:.") | blank_p) >> eol_p;
    }

    const rule_t& start() const { return aln; }
  };
};

template < class Seq >
bool
load_aln(MASequence<Seq>& ma, file_iterator<>& fi)
{
  aln_parser::WA wa;
  aln_parser parser(wa);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) return false;
  fi = info.stop;

  std::deque<std::string>::const_iterator x;
  for (x=wa.seqs.begin(); x!=wa.seqs.end(); ++x) {
    Seq r;
    char2rna(r, *x);
    ma.add_seq(r);
  }
  return true;
}

template < class Seq >
bool
load_aln(std::list<Seq>& ma, file_iterator<>& fi)
{
  aln_parser::WA wa;
  aln_parser parser(wa);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) return false;
  fi = info.stop;

  std::deque<std::string>::const_iterator x;
  for (x=wa.seqs.begin(); x!=wa.seqs.end(); ++x) {
    Seq r;
    char2rna(r, *x);
    ma.push_back(r);
  }
  return true;
}

template < class Seq >
bool
load_aln(std::list< MASequence<Seq> >& ma, const char* filename)
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
    MASequence<Seq> m;
    if (load_aln(m, fi)) {
      n++;
      ma.push_back(m);
    } else {
      break;
    }
  }
  return n!=0;
}

#if 0
template
bool
load_aln(std::list< MASequence<RNASequence> >& ma, const char* filename);
#endif

template
bool
load_aln(std::list< MASequence<std::string> >& ma, const char* filename);

template
bool
load_aln(MASequence<std::string>& ma, file_iterator<>& fi);

template
bool
load_aln(std::list<std::string>& ma, file_iterator<>& fi);
