// $Id$

#ifndef __INC_PROFILE_H__
#define __INC_PROFILE_H__

#include <vector>
#include <string>
#include <utility>
#include <map>
#include "rna.h"

class ProfileSequence
{
public:
  typedef float value_type;
  typedef std::vector<value_type> Column;
  typedef std::vector<Column>::const_iterator const_iterator;

public:
  ProfileSequence()
    : n_seqs_(0), profile_()
  {
  }

  ProfileSequence(const ProfileSequence& s)
    : n_seqs_(s.n_seqs_), profile_(s.profile_)
  {
  }

  ProfileSequence(const std::string& seq)
    : n_seqs_(0), profile_()
  {
    initialize(seq.size());
    add_sequence(seq);
  }

  template < class It >
  ProfileSequence(It b, It e)
    : n_seqs_(0), profile_()
  {
    initialize(b->size());
    for (It x=b; x!=e; ++x) add_sequence(*x);
  }

  template < class T >
  ProfileSequence(const T& ma)
    : n_seqs_(0), profile_()
  {
    initialize(ma.begin()->size());
    typename T::const_iterator x;
    for (x=ma.begin(); x!=ma.end(); ++x) add_sequence(*x);
  }

  ProfileSequence& operator=(const ProfileSequence& s);

  value_type n_seqs() const { return n_seqs_; }

  uint size() const { return profile_.size(); }

  const Column& operator[](uint i) const { return profile_[i]; }

  const_iterator begin() const { return profile_.begin(); }

  const_iterator end() const { return profile_.end(); }

private:
  void initialize(uint sz);

  void add_sequence(const std::string& seq, value_type w=1.0);

private:  
  value_type n_seqs_;
  std::vector<Column> profile_;
};

class BPProfileMaker
{
public:
  typedef float value_type;
  typedef std::pair<rna_t, rna_t> bp_t;

public:
  BPProfileMaker(const std::string& seq)
    : n_seqs_(1), ma_(seq.size())
  {
    for (uint i=0; i!=ma_.size(); ++i) ma_[i].resize(n_seqs_);
    add_sequence(seq, 0);
  }

  template < class It >
  BPProfileMaker(It b, It e)
    : n_seqs_(std::distance(b,e)), ma_(b->size())
  {
    for (uint i=0; i!=ma_.size(); ++i) ma_[i].resize(n_seqs_);
    uint n=0;
    for (It x=b; x!=e; ++x) add_sequence(*x, n++);
  }

  template < class T >
  BPProfileMaker(const T& ma)
    : n_seqs_(ma.size()), ma_(ma.begin()->size())
  {
    for (uint i=0; i!=ma_.size(); ++i) ma_[i].resize(n_seqs_);
    uint n=0;
    typename T::const_iterator x;
    for (x=ma.begin(); x!=ma.end(); ++x) add_sequence(*x, n++);
  }

  template < class T >
  BPProfileMaker(const MASequence<T>& ma);

  void make_profile(value_type p, uint i, uint j,
		    std::map<bp_t, value_type>& v) const;

  void make_profile(const std::vector<value_type>& p, uint i, uint j,
		    std::map<bp_t, value_type>& v) const;

  uint n_seqs() const { return n_seqs_; }

private:
  void add_sequence(const std::string& seq, uint n);
  
private:
  uint n_seqs_;
  std::vector< std::vector<rna_t> > ma_;
};

#endif	// __INC_PROFILE_H__

// Local Variables:
// mode: C++
// End:
