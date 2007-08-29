// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <list>
#include <stack>
#include <algorithm>
#include <cctype>
#include <cassert>
#include "rna.h"

static
char
char_table[N_IUPAC+1] = {
  RNASymbol<char>::A, RNASymbol<char>::C,
  RNASymbol<char>::G, RNASymbol<char>::U,
  RNASymbol<char>::T, RNASymbol<char>::GAP,
  RNASymbol<char>::R, RNASymbol<char>::Y,
  RNASymbol<char>::M, RNASymbol<char>::K,
  RNASymbol<char>::S, RNASymbol<char>::W,
  RNASymbol<char>::B, RNASymbol<char>::D,
  RNASymbol<char>::H, RNASymbol<char>::V,
  RNASymbol<char>::N,
};

static
rna_t
rna_table[N_IUPAC+1] = {
  RNASymbol<rna_t>::A, RNASymbol<rna_t>::C,
  RNASymbol<rna_t>::G, RNASymbol<rna_t>::U,
  RNASymbol<rna_t>::T, RNASymbol<rna_t>::GAP,
  RNASymbol<rna_t>::R, RNASymbol<rna_t>::Y,
  RNASymbol<rna_t>::M, RNASymbol<rna_t>::K,
  RNASymbol<rna_t>::S, RNASymbol<rna_t>::W,
  RNASymbol<rna_t>::B, RNASymbol<rna_t>::D,
  RNASymbol<rna_t>::H, RNASymbol<rna_t>::V,
  RNASymbol<rna_t>::N,
};

rna_t
char2rna(char rna)
{
#if 0
  switch (tolower(rna)) {
  case RNASymbol<char>::A:
    return static_cast<rna_t>(RNASymbol<rna_t>::A);
    break;
  case RNASymbol<char>::C:
    return static_cast<rna_t>(RNASymbol<rna_t>::C);
    break;
  case RNASymbol<char>::G:
    return static_cast<rna_t>(RNASymbol<rna_t>::G);
    break;
  case RNASymbol<char>::U:
  case RNASymbol<char>::T:
    return static_cast<rna_t>(RNASymbol<rna_t>::U);
    break;
  default:
  case RNASymbol<char>::GAP:
    return static_cast<rna_t>(RNASymbol<rna_t>::GAP);
    break;
  }
#else
  char r = tolower(rna);
  for (uint i=0; i!=N_IUPAC+1; ++i) {
    if (r==char_table[i]) return rna_table[i];
  }
  return static_cast<rna_t>(RNASymbol<rna_t>::GAP);
#endif
}

char
rna2char(const rna_t& rna)
{
#if 0
  switch (rna) {
  case RNASymbol<rna_t>::A:
    return RNASymbol<char>::A;
    break;
  case RNASymbol<rna_t>::C:
    return RNASymbol<char>::C;
    break;
  case RNASymbol<rna_t>::G:
    return RNASymbol<char>::G;
    break;
  case RNASymbol<rna_t>::U:
    return RNASymbol<char>::U;
    break;
  default:
  case RNASymbol<rna_t>::GAP:
    return RNASymbol<char>::GAP;
    break;
  }
#else
  for (uint i=0; i!=N_IUPAC+1; ++i) {
    if (rna==rna_table[i]) return char_table[i];
  }
  return RNASymbol<char>::GAP;
#endif
}

struct c2r {
  RNASequence::value_type operator()(const std::string::value_type& c) const
  {
    return char2rna(c);
  }
};

struct r2c {
  std::string::value_type operator()(const RNASequence::value_type& c) const
  {
    return rna2char(c);
  }
};

void
char2rna(RNASequence& out, const std::string& in)
{
  out.resize(in.size());
  std::transform(in.begin(), in.end(), out.begin(), c2r());
}

void
char2rna(std::string& out, const std::string& in)
{
  out.resize(in.size());
  std::copy(in.begin(), in.end(), out.begin());
}

void
rna2char(std::string& out, const RNASequence& in)
{
  out.resize(in.size());
  std::transform(in.begin(), in.end(), out.begin(), r2c());
}

void
rna2char(std::string& out, const std::string& in)
{
  out.resize(in.size());
  std::copy(in.begin(), in.end(), out.begin());
}

template < class V >
Column<V>::
Column(const std::deque<V>& val)
  : val_(val), cnt_(N_RNA, 0)
{
  typename std::deque<V>::const_iterator x;
  for (x=val_.begin(); x!=val_.end(); ++x) {
    if (*x != RNASymbol<V>::GAP)
      cnt_[index(*x)]++;
  }
}

template < class V >
void
Column<V>::
push_back(const V& v)
{
  val_.push_back(v);
  if (v != RNASymbol<V>::GAP)
    cnt_[index(v)]++;
}

template < class Seq >
Seq
MASequence<Seq>::
get_seq(uint i) const
{
  Seq r;
  r.resize(this->size());
  typename MASequence<Seq>::const_iterator x;
  typename Seq::iterator y;
  for (x=this->begin(), y=r.begin(); x!=this->end(); ++x, ++y)
    *y = (*x)[i];
  return r;
}

template
RNASequence
MASequence<RNASequence>::get_seq(uint i) const;

template
std::string
MASequence<std::string>::get_seq(uint i) const;


template < class Seq >
Column<typename Seq::value_type>
gap_symbol(const MASequence<Seq>& ma_seq)
{
  std::deque<typename Seq::value_type> v(ma_seq[0].n_seqs());
  typename Seq::value_type GAP = RNASymbol<typename Seq::value_type>::GAP;
  std::fill(v.begin(), v.end(), GAP);
  Column<typename Seq::value_type> c(v);
  return c;
}

template
Column<RNASequence::value_type>
gap_symbol(const MASequence<RNASequence>& ma_seq);

template
Column<std::string::value_type>
gap_symbol(const MASequence<std::string>& ma_seq);


template <class Seq>
Seq
erase_gap(const Seq& seq)
{
  const typename Seq::value_type GAP = RNASymbol<typename Seq::value_type>::GAP;
  Seq ret;
  typename Seq::const_iterator x;
  for (x=seq.begin(); x!=seq.end(); ++x) {
    if (*x!=GAP) ret.push_back(*x);
  }
  return ret;
}

template
RNASequence
erase_gap(const RNASequence& seq);

template
std::string
erase_gap(const std::string& seq);


static
void
make_paren_list(const std::string& str,
		std::list< std::pair<uint,uint> >& pa_list)
{
  std::stack<uint> st;
  for (uint i=0; i!=str.size(); ++i) {
    if (str[i]=='(') {
      st.push(i);
    } else if (str[i]==')') {
      pa_list.push_back(std::make_pair(st.top(), i));
      st.pop();
    }
  }
  assert(st.empty());
}

static
void
make_paren_map(const std::string& str, std::vector<uint>& pa_map)
{
  std::list< std::pair<uint,uint> > pa_list;
  make_paren_list(str, pa_list);
  pa_map.resize(str.size());
  std::fill(pa_map.begin(), pa_map.end(), static_cast<uint>(-1));
  std::list< std::pair<uint,uint> >::const_iterator x;
  for (x=pa_list.begin(); x!=pa_list.end(); ++x) {
    pa_map[x->first]  = x->second;
    pa_map[x->second] = x->first;
  }
}

template <class Seq>
void
erase_gap(const Seq& seq, const std::string& str,
	  Seq& r_seq, std::string& r_str)
{
  const typename Seq::value_type A = RNASymbol<typename Seq::value_type>::A;
  const typename Seq::value_type C = RNASymbol<typename Seq::value_type>::C;
  const typename Seq::value_type G = RNASymbol<typename Seq::value_type>::G;
  const typename Seq::value_type U = RNASymbol<typename Seq::value_type>::U;
  const typename Seq::value_type GAP = RNASymbol<typename Seq::value_type>::GAP;

  assert(seq.size()==str.size());
  std::vector<uint> pa_map;
  make_paren_map(str, pa_map);

  for (uint i=0; i!=str.size(); ++i) {
    if (pa_map[i]!=static_cast<uint>(-1)) {
      if (seq[i]==GAP) {
	pa_map[pa_map[i]]=static_cast<uint>(-1);
	pa_map[i]=static_cast<uint>(-1);
      } else if (seq[i]==A && seq[pa_map[i]]==U ||
		 seq[i]==U && seq[pa_map[i]]==A ||
		 seq[i]==G && seq[pa_map[i]]==C ||
		 seq[i]==C && seq[pa_map[i]]==G ||
		 seq[i]==G && seq[pa_map[i]]==U ||
		 seq[i]==U && seq[pa_map[i]]==G ) {
      } else {
	pa_map[pa_map[i]]=static_cast<uint>(-1);
	pa_map[i]=static_cast<uint>(-1);
      }
    }
  }

  for (uint i=0; i!=str.size(); ++i) {
    if (seq[i]!=GAP) {
      r_seq.push_back(seq[i]);
      if (pa_map[i]==static_cast<uint>(-1)) {
	r_str.push_back(' ');
      } else if (pa_map[i]>i) {
	r_str.push_back('(');
      } else {
	r_str.push_back(')');
      }
    }
  }
}

// instantiation

template
class Column<std::string::value_type>;

template
class Column<RNASequence::value_type>;

template
void
erase_gap(const RNASequence& seq, const std::string& str,
	  RNASequence& r_seq, std::string& r_str);

template
void
erase_gap(const std::string& seq, const std::string& str,
	  std::string& r_seq, std::string& r_str);

template < class Seq >
MASequence<Seq>&
MASequence<Seq>::
add_seq(const Seq& seq)
{
  if (this->empty()) this->resize(seq.size());
  assert(this->size()==seq.size());
  typename Seq::const_iterator x;
  typename MASequence<Seq>::iterator y;
  for (x=seq.begin(), y=this->begin(); x!=seq.end(); ++x, ++y) {
    y->push_back(*x);
  }
  return *this;
}

template
MASequence<RNASequence>&
MASequence<RNASequence>::add_seq(const RNASequence& seq);

template
MASequence<std::string>&
MASequence<std::string>::add_seq(const std::string& seq);
