// $Id:$

#include <map>
#include <cassert>
#include "dag.h"

using namespace DAG;

template < class T, class V >
static
void
make_cnt(std::list< std::pair< std::pair<rna_t,rna_t>, V > >& cnt,
	 const T& c_first, const T& c_last)
{
  std::pair<rna_t,rna_t> k = std::make_pair(index(c_first), index(c_last));
  cnt.push_back(std::make_pair(k,1));
}

template < class T, class V >
static
void
make_cnt(std::list< std::pair< std::pair<rna_t,rna_t>, V > >& cnt,
	 const Column<T>& c_first, const Column<T>& c_last)
{
  assert(c_first.n_seqs()==c_last.n_seqs());
  T GAP=RNASymbol<T>::GAP;
  std::map<std::pair<rna_t,rna_t>, uint> c;
  for (uint i=0; i!=c_first.n_seqs(); ++i) {
    if (c_first[i]!=GAP && c_last[i]!=GAP) {
      std::pair<rna_t,rna_t> k=std::make_pair(index(c_first[i]),
					      index(c_last[i]));
      std::map<std::pair<rna_t,rna_t>, uint>::iterator x = c.find(k);
      if (x==c.end()) {
	c.insert(std::make_pair(k,1));
      } else {
	x->second++;
      }
    }
  }
  std::map<std::pair<rna_t,rna_t>, uint>::const_iterator x;
  for (x=c.begin(); x!=c.end(); ++x) {
    cnt.push_back(*x);
  }
}
    
template < class E >
template < class C >
Node<E>::
Node(uint first, uint last, const C& c_first, const C& c_last,
     float weight, uint n_edges)
  : first_(first), last_(last), edges_(n_edges), weight_(weight), cnt_()
{
  make_cnt(cnt_, c_first, c_last);
}

template < class E >
template < class C >
Node<E>::
Node(const Pos& pos, const C& c_first, const C& c_last,
     float weight, uint n_edges)
  : first_(pos.first), last_(pos.second),
    edges_(n_edges), weight_(weight), cnt_()
{
  make_cnt(cnt_, c_first, c_last);
}

// instantiation

template
class Node<Edge>;

#if 0
template
Node<Edge>::
Node(const Pos& pos, const rna_t& c_first, const rna_t& c_last, float weight);

template
Node<Edge>::
Node(const Pos& pos, const Column<rna_t>& c_first, const Column<rna_t>& c_last,
     float weight);
#endif

template
Node<Edge>::
Node(const Pos& pos, const char& c_first, const char& c_last,
     float weight, uint n_edges);

template
Node<Edge>::
Node(const Pos& pos, const Column<char>& c_first, const Column<char>& c_last,
     float weight, uint n_edges);
