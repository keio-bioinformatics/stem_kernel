// $Id$

#ifndef __INC_SCORE_TABLE_H__
#define __INC_SCORE_TABLE_H__

#include <boost/multi_array.hpp>
#include "dag.h"

template < class D, class V >
class SimpleNodeScore
{
public:
  typedef D Data;
  typedef V value_type;
  typedef typename Data::Seq Seq;
  typedef typename Data::Node Node;
  typedef typename Data::Edge Edge;

  SimpleNodeScore(value_type gap, value_type stack, value_type covar)
    : gap_(gap), stack_(stack), covar_(covar)
  {
  }

  value_type node_score(const Data& xx, const Data& yy,	uint i, uint j) const;

  value_type node_score(const Data& xx, uint i) const
  {
    return gap_*gap_;
  }

private:
  value_type gap_;
  value_type stack_;
  value_type covar_;
};

template < class D, class V >
class SubstNodeScore
{
public:
  typedef D Data;
  typedef V value_type;
  typedef typename Data::Seq Seq;
  typedef typename Data::Node Node;
  typedef typename Data::Edge Edge;

  SubstNodeScore(value_type gap, value_type beta);

  value_type node_score(const Data& xx, const Data& yy,	uint i, uint j) const;

  value_type node_score(const Data& xx, uint i) const
  {
    return gap_*gap_;
  }

private:
  value_type gap_;
  boost::multi_array<value_type,4> co_subst_;
  value_type beta_;
};

template < class D, class V >
class SimpleEdgeScore
{
public:
  typedef D Data;
  typedef V value_type;
  typedef typename Data::Seq Seq;
  typedef typename Data::Node Node;
  typedef typename Data::Edge Edge;
  struct WorkArea {
    //std::vector<value_type> gap;
  };

  SimpleEdgeScore(value_type gap) : gap_(gap) { }

  void initialize(WorkArea& wa, const Data& xx, const Data& yy) const;

  value_type edge_score(const WorkArea& wa,
			const Data& xx, const Data& yy,
			const Edge& ix, const Edge& iy,
			uint i, uint j) const;

  value_type edge_score(const WorkArea& wa,
			const Data& xx, const Edge& ix, uint i) const;

  value_type edge_ext_score(const WorkArea& wa,
			    const Data& xx, const Data& yy,
			    const Edge& ix, const Edge& iy,
			    uint i, uint j) const;

private:
  value_type gap_;

  static std::vector<value_type> gap_vec_;
};

template < class D, class V >
class SubstEdgeScore
{
public:
  typedef D Data;
  typedef V value_type;
  typedef typename Data::Seq Seq;
  typedef typename Data::Node Node;
  typedef typename Data::Edge Edge;
  struct WorkArea {
    std::vector<value_type> gap;
  };

  SubstEdgeScore(value_type gap, value_type beta);

  void initialize(WorkArea& wa, const Data& xx, const Data& yy) const;

  value_type edge_score(const WorkArea& wa,
			const Data& xx, const Data& yy,
			const Edge& ix, const Edge& iy,
			uint i, uint j) const;

  value_type edge_score(const WorkArea& wa,
			const Data& xx, const Edge& ix, uint i) const;

  value_type edge_ext_score(const WorkArea& wa,
			    const Data& xx, const Data& yy,
			    const Edge& ix, const Edge& iy,
			    uint i, uint j) const;

private:
  value_type gap_;
  value_type beta_;
  boost::multi_array<value_type,2> si_subst_;
};


template < class D, class V >
class SimpleScoreTable
  : public SimpleNodeScore<D,V>, public SimpleEdgeScore<D,V>
{
public:
  typedef V value_type;
  typedef D Data;
  typedef typename Data::Seq Seq;
  typedef typename SimpleEdgeScore<D,V>::WorkArea WorkArea;

  SimpleScoreTable(value_type gap, value_type stack, value_type covar,
		   value_type loop_gap)
    : SimpleNodeScore<D,V>(gap, stack, covar),
      SimpleEdgeScore<D,V>(loop_gap)
  {
  }
};

template < class D, class V >
class SubstScoreTable
  : public SubstNodeScore<D,V>, public SimpleEdgeScore<D,V>
{
public:
  typedef V value_type;
  typedef D Data;
  typedef typename Data::Seq Seq;
  typedef typename Data::Node Node;
  typedef typename Data::Edge Edge;
  typedef typename SimpleEdgeScore<D,V>::WorkArea WorkArea;

  SubstScoreTable(value_type gap, value_type beta, value_type loop_gap)
    : SubstNodeScore<D,V>(gap, beta),
      SimpleEdgeScore<D,V>(loop_gap)
  {
  }
};

#endif	// __INC_SCORE_TABLE_H__

// Local Variables:
// mode: C++
// End:
