// $Id:$

#include <map>
#include <cassert>
#include "dag.h"
#include "bpmatrix.h"

using namespace DAG;

#if 0
template < class BPM >
static
void
make_cnt(std::list<bp_freq_t>& bp_freq, uint i, uint j,
	 const std::string& seq, const BPM& bpm)
{
  float f = bpm(i+1, j+1);
  if (f==0.0) return;
  for (uint a=0; a!=N_RNA; ++a) {
    float aw=iupac_weight[index(seq[i])][a];
    if (aw==0.0) continue;
    for (uint b=0; b!=N_RNA; ++b) {
      float bw=iupac_weight[index(seq[j])][b];
      if (bw==0.0) continue;
      bp_freq.push_back(std::make_pair(std::make_pair(a,b), f*aw*bw));
    }
  }
}

template < class BPM >
static
void
make_cnt(std::list<bp_freq_t>& bp_freq, uint i, uint j,
	 const MASequence<std::string>& seq, const BPM& bpm)
{
  char GAP=RNASymbol<char>::GAP;
  std::map<bp_t, float> c;
  typename BPM::matrix_iterator mtx = bpm.matrix_begin();
  typename BPM::idx_map_iterator idx = bpm.idx_map_begin();
  assert(seq[i].n_seqs() == bpm.n_matrices());

  for (uint x=0; x!=seq[i].n_seqs(); ++x, ++mtx, ++idx) {
    if (seq[i][x]!=GAP && seq[j][x]!=GAP) {
      float f = (**mtx)((*idx)[i]+1, (*idx)[j]+1);
      if (f==0.0) continue;
      for (uint a=0; a!=N_RNA; ++a) {
	float aw=iupac_weight[index(seq[i][x])][a];
	if (aw==0.0) continue;
	for (uint b=0; b!=N_RNA; ++b) {
	  float bw=iupac_weight[index(seq[j][x])][b];
	  if (bw==0.0) continue;
	  bp_t k = std::make_pair(a,b);
	  std::map<bp_t, float>::iterator y = c.find(k);
	  if (y==c.end()) {
	    c.insert(std::make_pair(k, f*aw*bw));
	  } else {
	    y->second += f*aw*bw;
	  }
	}
      }
    }
  }

  std::map<bp_t, float>::iterator y;
  for (y=c.begin(); y!=c.end(); ++y) {
    bp_freq.push_back(*y);
  }
}
#endif

#if 0
template < class E >
template < class Seq, class BPM >
Node<E>::
Node(uint first, uint last, const Seq& seq, const BPM& bpm, uint n_edges)
  : first_(first), last_(last), edges_(n_edges), weight_(1.0), bp_freq_()
{
  make_cnt(bp_freq_, first_, last_, seq, bpm);
}

template < class E >
template < class Seq, class BPM >
Node<E>::
Node(const Pos& pos, const Seq& seq, const BPM& bpm, uint n_edges)
  : first_(pos.first), last_(pos.second),
    edges_(n_edges), weight_(1.0), bp_freq_()
{
  make_cnt(bp_freq_, first_, last_, seq, bpm);
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

#if 0
template
Node<Edge>::
Node(const Pos& pos, const char& c_first, const char& c_last,
     float weight, uint n_edges);

template
Node<Edge>::
Node(const Pos& pos, const Column<char>& c_first, const Column<char>& c_last,
     float weight, uint n_edges);

#else

template
Node<Edge>::
Node(const Pos& pos, const std::string& seq,
     const BPMatrix& bpm, uint n_edges);

template
Node<Edge>::
Node(const Pos& pos, const std::list<std::string>& seq,
     const BPMatrix& bpm, uint n_edges);
#endif

#endif
