// $Id:$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "maf.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>

using namespace BOOST_SPIRIT_CLASSIC_NS;

struct maf_parser : public grammar< maf_parser >
{
  std::list<std::string>& seqs;

  maf_parser(std::list<std::string>& x) : seqs(x) {}

  template <class ScannerT>
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t maf;
    rule_t header;
    rule_t comment;
    rule_t empty;
    rule_t ali;
    rule_t seq, seq_s1, seq_s2, seq_i;

    definition(const maf_parser& self)
    {
      maf = !header >> *(comment | empty) >> ali >> +seq >> *empty;
      header = str_p("##maf") >> *print_p >> eol_p;
      comment = comment_p("#");
      empty = *blank_p >> eol_p;
      ali = ch_p('a') >> +blank_p >> *print_p >> eol_p;
      seq = ((seq_s1 >> seq_s2) | seq_i) >> eol_p;
      seq_s1 = ch_p('s') >> +blank_p >> +graph_p >> +blank_p
			 >> uint_p >> +blank_p >> uint_p >> +blank_p;
      seq_s2 = sign_p >> +blank_p >> uint_p >> +blank_p
		      >> (+graph_p)[push_back_a(self.seqs)] >> *blank_p;
      seq_i = (ch_p('i') | ch_p('e')) >> +blank_p >> +print_p >> eol_p;
    }

    const rule_t& start() const { return maf; }
  };
};

template < class Seq >
bool
load_maf(MASequence<Seq>& ma, file_iterator<>& fi)
{
  file_iterator<> s = fi;
  std::list<std::string> seqs;
  maf_parser parser(seqs);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) {
    fi = s;
    return false;
  }
  fi = info.stop;

  MASequence<Seq> m;
  std::list<std::string>::const_iterator x;
  for (x=seqs.begin(); x!=seqs.end(); ++x) {
    Seq r;
    char2rna(r, *x);
    ma.add_seq(r);
  }
  return true;
}

template < class Seq >
bool
load_maf(std::list<Seq>& ma, file_iterator<>& fi)
{
  file_iterator<> s = fi;
  std::list<std::string> seqs;
  maf_parser parser(seqs);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) {
    fi = s;
    return false;
  }
  fi = info.stop;

  MASequence<Seq> m;
  std::list<std::string>::const_iterator x;
  for (x=seqs.begin(); x!=seqs.end(); ++x) {
    Seq r;
    char2rna(r, *x);
    ma.push_back(r);
  }
  return true;
}

template < class Seq >
bool
load_maf(std::list< MASequence<Seq> >& ma, const char* filename)
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
    if (load_maf(m, fi)) {
      n++;
      ma.push_back(m);
    } else {
      break;
    }
  }
  return n!=0;
}

void
MAF::
init()
{
  while (std::getline(*in_, cur_)) {
    if (cur_[0]=='a') return;
  }
  cur_.clear();
}

template < class Seq >
bool
MAF::
get(MASequence<Seq>& ma)
{
  uint n=0;
  while (std::getline(*in_, cur_)) {
    if (cur_[0]=='a') break;
    if (cur_[0]=='s') {
      std::istringstream s(cur_);
      std::vector<std::string> v;
      std::copy(std::istream_iterator<std::string>(s),
		std::istream_iterator<std::string>(),
		std::back_inserter(v));
      if (v.size()==7) {
	Seq r;
	char2rna(r, v[6]);
	ma.add_seq(r);
	++n;
      }
    }
  }
  return n>0;
}

#if 0
template
bool
load_maf(std::list< MASequence<RNASequence> >& ma, const char* filename);
#else
template
bool
load_maf(std::list< MASequence<std::string> >& ma, const char* filename);
#endif

template
bool
MAF::
get(MASequence<std::string>& ma);

template
bool
load_maf(MASequence<std::string>& ma, file_iterator<>& fi);

template
bool
load_maf(std::list<std::string>& ma, file_iterator<>& fi);
