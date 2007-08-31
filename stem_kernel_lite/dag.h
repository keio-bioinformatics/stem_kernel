// $Id$

#ifndef __INC_DAG_H__
#define __INC_DAG_H__

#include <vector>
#include <list>
#include <iostream>
#include "../common/cyktable.h"
#include "../common/rna.h"

namespace DAG {
  typedef std::pair<rna_t, rna_t> bp_t;
  typedef std::pair<bp_t, float> bp_freq_t;
  typedef std::list<bp_freq_t>::const_iterator bp_freq_iterator;

  class Edge
  {
  public:
    Edge() : to_(static_cast<uint>(-1)), gaps_(0), weight_(1.0) { }

    Edge(uint to, const Pos& p_pos, const Pos& c_pos, float weight=1.0)
      : to_(to), gaps_(0), p_pos_(p_pos), c_pos_(c_pos), weight_(weight)
    {
      gaps_  = c_pos_.first - p_pos_.first - 1;
      gaps_ += p_pos_.second - c_pos_.second - 1;
    }

    Edge(uint to, const Pos& p_pos, float weight=1.0)
      : to_(to), gaps_(0), p_pos_(p_pos), c_pos_(p_pos), weight_(weight)
    {
      gaps_ = p_pos_.second - p_pos_.first - 1;
    }

    Edge(const Edge& e)
      : to_(e.to_), gaps_(e.gaps_),
	p_pos_(e.p_pos_), c_pos_(e.c_pos_),
	weight_(e.weight_)
    { }

    Edge& operator=(const Edge& e)
    {
      if (this != &e) {
	to_ = e.to_;
	gaps_ = e.gaps_;
	p_pos_ = e.p_pos_;
	c_pos_ = e.c_pos_;
	weight_ = e.weight_;
      }
      return *this;
    }

    uint gaps() const { return gaps_; }
    uint to() const { return to_; }
    const Pos& p_pos() const { return p_pos_; }
    const Pos& c_pos() const { return c_pos_; }
    float weight() const { return weight_; }
      
  private:
    uint to_;
    uint gaps_;
    Pos p_pos_;
    Pos c_pos_;
    float weight_;
  };

  template < class E = Edge >
  class Node
  {
  public:
    typedef E Edge;
    typedef typename std::vector<Edge>::iterator iterator;
    typedef typename std::vector<Edge>::const_iterator const_iterator;
      
  public:
    Node()
      : first_(), last_(), edges_(), weight_(), bp_freq_()
    {
    }

    Node(uint first, uint last, const std::list<bp_freq_t>& bp_freq,
	 uint n_edges=0, float weight=1.0)
      : first_(first), last_(last),
	edges_(n_edges), weight_(weight), bp_freq_(bp_freq)
    {
    }

    Node(const Pos& pos, const std::list<bp_freq_t>& bp_freq,
	 uint n_edges=0, float weight=1.0)
      : first_(pos.first), last_(pos.second),
	edges_(n_edges), weight_(weight), bp_freq_(bp_freq)
    {
    }

    Node(uint first, uint last, float weight=1.0)
      : first_(first), last_(last),
	edges_(), weight_(weight), bp_freq_()
    {
    }

    Node(const Pos& pos, float weight=1.0)
      : first_(pos.first), last_(pos.second),
	edges_(), weight_(weight), bp_freq_()
    {
    }

    Node(const Node& n)
      : first_(n.first_), last_(n.last_),
	edges_(n.edges_), weight_(n.weight_), bp_freq_(n.bp_freq_)
    {
    }

    Node& operator=(const Node& n)
    {
      if (this != &n) {
	first_ = n.first_;
	last_ = n.last_;
	edges_ = n.edges_;
	weight_ = n.weight_;
	bp_freq_ = n.bp_freq_;
      }
      return *this;
    }

    void push_back(const Edge& edge)
    {
      edges_.push_back(edge);
    }

    const_iterator begin() const { return edges_.begin(); }
    const_iterator end() const { return edges_.end(); }
    iterator begin() { return edges_.begin(); }
    iterator end() { return edges_.end(); }
    bool empty() const { return edges_.empty(); }
    uint size() const { return edges_.size(); }
    const Edge& operator[](uint i) const { return edges_[i]; }
    Edge& operator[](uint i) { return edges_[i]; }

    uint first() const { return first_; }
    uint last() const { return last_; }
    float weight() const { return weight_; }
    //const std::list<bp_freq_t>& cnt() const { return bp_freq_; }
    bp_freq_iterator bp_freq_begin() const { return bp_freq_.begin(); }
    bp_freq_iterator bp_freq_end() const { return bp_freq_.end(); }
    uint length() const { return last_-first_; }

  private:
    uint first_;
    uint last_;
    std::vector<Edge> edges_;
    float weight_;
    std::list<bp_freq_t> bp_freq_;
  };
};

#endif	// __INC_DAG_H__

// Local Variables:
// mode: C++
// End:
