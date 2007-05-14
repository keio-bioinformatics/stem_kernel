// $Id$

#include "aln.h"
#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <boost/spirit.hpp>

using namespace boost::spirit;

struct aln_parser : public grammar< aln_parser >
{
  struct WA
  {
    std::string cur_name;
    std::string cur_seq;
    std::list<std::string> names;
    std::map<std::string, std::string> name_seq_map;
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
      if (wa.name_seq_map.find(wa.cur_name)==wa.name_seq_map.end()) {
	wa.names.push_back(wa.cur_name);
	wa.name_seq_map.insert(std::make_pair(wa.cur_name,std::string()));
      }
      wa.name_seq_map[wa.cur_name] += wa.cur_seq;
    }
  };

  template <class ScannerT>
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t aln;
    rule_t header;
    rule_t empty;
    rule_t seq;
    rule_t status;
    
    definition(const aln_parser& self)
    {
      aln = header >> *empty >> +(+seq[push_seq(self.wa)] >> !status >> *empty);
      header = (str_p("CLUSTAL") | str_p("PROBCONS")) >> +print_p >> eol_p;
      empty = *blank_p >> eol_p;
      seq
	= (+graph_p)[assign_a(self.wa.cur_name)]
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

  std::list<std::string>::const_iterator x;
  for (x=wa.names.begin(); x!=wa.names.end(); ++x) {
    Seq r;
    char2rna(r, wa.name_seq_map[*x]);
    ma.add_seq(r);
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
