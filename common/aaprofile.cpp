// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cctype>
#include <cassert>
#include <cstring>
#include "aaprofile.h"

//static
int 
AAProfileSequence::
char2aa(int c)
{
  static const char* aa_char = "ARNDCQEGHILKMFPSTWYVBZX";
  assert(strlen(aa_char)==N_AA);
  const char* p;
  for (p=aa_char; *p!=0; ++p)
    if (*p==toupper(c)) return p-aa_char;
  return p-aa_char;
}

// public
AAProfileSequence&
AAProfileSequence::
operator=(const AAProfileSequence& s)
{
  if (this != &s) {
    n_seqs_ = s.n_seqs_;
    profile_ = s.profile_;
  }
  return *this;
}

// private
void
AAProfileSequence::
initialize(unsigned int sz)
{
  profile_.resize(sz);
  for (unsigned int i=0; i!=profile_.size(); ++i)
    profile_[i].resize(N_AA+1, 0.0);
  n_seqs_ = 0;
}

// public
void
AAProfileSequence::
add_sequence(const std::string& seq, value_type w /*=1.0*/)
{
  if (profile_.size() != seq.size())
    throw "alignment error";
    
  for (unsigned int i=0; i!=profile_.size(); ++i) {
    profile_[i][char2aa(seq[i])] += w;
  }
  n_seqs_ += w;
}

// public
void
AAProfileSequence::
add_sequence(const AAProfileSequence& seq, value_type w /*=1.0*/)
{
  if (profile_.size() != seq.size())
    throw "alignment error";

  for (unsigned int i=0; i!=profile_.size(); ++i) {
    for (unsigned int j=0; j!=N_AA+1; ++j) {
      profile_[i][j] += seq[i][j]*w;
    }
  }
  n_seqs_ += w;
}

