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

  void add_sequence(const std::string& seq, value_type w=1.0);

  void add_sequence(const ProfileSequence& seq, value_type w=1.0);

private:
  void initialize(uint sz);

private:  
  value_type n_seqs_;
  std::vector<Column> profile_;
};

#endif	// __INC_PROFILE_H__

// Local Variables:
// mode: C++
// End:
