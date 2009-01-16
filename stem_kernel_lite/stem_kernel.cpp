// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <vector>
#include <list>
#include <cassert>
#include "dag.h"
#include "dptable.h"
#include "stem_kernel.h"

template < class ST, class D >
typename StemKernel<ST,D>::value_type
StemKernel<ST,D>::
operator()(const Data& xx, const Data& yy) const
{
  typedef typename Data::Node Node;
  typename ScoreTable::WorkArea wa;
  st_.initialize(wa, xx, yy);

  const std::vector<Node>& x(xx.tree);
  const std::vector<Node>& y(yy.tree);
  const std::vector<uint> root_x(xx.root);
  const std::vector<uint> root_y(yy.root);
  const std::vector<uint>& max_pa(xx.max_pa);
  DPTable<value_type> K0(x.size(), y.size());
  DPTable<value_type> G0(x.size(), y.size());
  std::vector<value_type> K1(y.size());
  std::vector<value_type> G1(y.size());
  typename Node::const_iterator ix, iy;

  for (uint i=0; i!=x.size(); ++i) {
    K0.allocate(i);
    G0.allocate(i);
    for (uint j=0; j!=y.size(); ++j) {
      // initialize
      if (x[i].empty() && y[j].empty()) {
	K0[i][j] = G0[i][j] = 1.0;
	continue;
      }

      // MATCH
      K1[j] = G1[j] = 0.0;
      if (!x[i].empty() && !y[j].empty() &&
	  (len_band_==0 ||
	   static_cast<uint>(abs(x[i].length()-y[j].length()))<=len_band_) ) {
	// score for this node
	value_type v_s = st_.node_score(xx, yy, i, j);
	for (ix=x[i].begin(); ix!=x[i].end(); ++ix) {
	  for (iy=y[j].begin(); iy!=y[j].end(); ++iy) {
	    value_type e_s = st_.edge_score(wa, xx, yy, *ix, *iy, i, j);
	    value_type v = G0[ix->to()][iy->to()]*v_s*e_s;
	    K1[j] += v;
	    G1[j] += v;
	  }
	}
      }

      // IY
      for (iy=y[j].begin(); iy!=y[j].end(); ++iy) {
	value_type v_s = st_.node_score(yy, j);
	value_type e_s = st_.edge_score(wa, yy, *iy, j);
	K1[j] += K1[iy->to()];
	G1[j] += G1[iy->to()]*v_s*e_s;
      }

      // IX
      K0[i][j] = K1[j];
      G0[i][j] = G1[j];
      for (ix=x[i].begin(); ix!=x[i].end(); ++ix) {
	value_type v_s = st_.node_score(xx, i);
	value_type e_s = st_.edge_score(wa, xx, *ix, i);
	K0[i][j] += K0[ix->to()][j];
	G0[i][j] += G0[ix->to()][j]*v_s*e_s;
      }
    }

    for (ix=x[i].begin(); ix!=x[i].end(); ++ix) {
      if (max_pa[ix->to()]<=i) {
	K0.deallocate(ix->to());
	G0.deallocate(ix->to());
      }
    }
  }

  value_type ret=0.0;
  for (uint i=0; i!=root_x.size(); ++i) {
    for (uint j=0; j!=root_y.size(); ++j) {
      ret += K0[root_x[i]][root_y[j]];
    }
  }
  return ret;
}

// instantiation
#include <string>
#include "score_table.h"
#include "data.h"
#include "../common/rna.h"

typedef double ValueType;

#if 0
template
class StemKernel<SimpleScoreTable<ValueType,SData>, SData>;

template
class StemKernel<SubstScoreTable<ValueType,SData>, SData>;
#endif

template
class StemKernel<SimpleScoreTable<ValueType,MData>, MData>;

template
class StemKernel<SubstScoreTable<ValueType,MData>, MData>;


