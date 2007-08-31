// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "profile.h"

static
float
iupac_weight[N_IUPAC][N_RNA] = {
  // A      C      G      T(U)
  {  1.0,   0.0,   0.0,   0.0  }, // A
  {  0.0,   1.0,   0.0,   0.0  }, // C
  {  0.0,   0.0,   1.0,   0.0  }, // G
  {  0.0,   0.0,   0.0,   1.0  }, // T(U)
  {  0.0,   0.0,   0.0,   0.0  }, // GAP
  { 1.0/2,  0.0,  1.0/2,  0.0  }, // R = A or G
  {  0.0,  1.0/2,  0.0,  1.0/2 }, // Y = C or T(U)
  { 1.0/2, 1.0/2,  0.0,   0.0  }, // M = A or C
  {  0.0,   0.0,  1.0/2, 1.0/2 }, // K = G or T(U)
  {  0.0,  1.0/2, 1.0/2,  0.0  }, // S = C or G
  { 1.0/2,  0.0,   0.0,  1.0/2 }, // W = A or T(U)
  {  0.0,  1.0/3, 1.0/3, 1.0/3 }, // B = C or G or T(U)
  { 1.0/3,  0.0,  1.0/3, 1.0/3 }, // D = A or G or T(U)
  { 1.0/3, 1.0/3,  0.0,  1.0/3 }, // H = A or C or T(U)
  { 1.0/3, 1.0/3, 1.0/3,  0.0  }, // V = A or C or G
  { 1.0/4, 1.0/4, 1.0/4, 1.0/4 }, // N = A or C or G or T(U)
};

// public
ProfileSequence&
ProfileSequence::
operator=(const ProfileSequence& s)
{
  if (this != &s) {
    n_seqs_ = s.n_seqs_;
    profile_ = s.profile_;
  }
  return *this;
}

// private
void
ProfileSequence::
initialize(uint sz)
{
  profile_.resize(sz);
  for (uint i=0; i!=profile_.size(); ++i)
    profile_[i].resize(N_RNA, 0.0);
  n_seqs_ = 0;
}

// public
void
ProfileSequence::
add_sequence(const std::string& seq, value_type w /*=1.0*/)
{
  if (profile_.size() != seq.size())
    throw "alignment error";
    
  for (uint i=0; i!=profile_.size(); ++i) {
    rna_t r = char2rna(seq[i]);
    if (r != RNASymbol<rna_t>::GAP) {
      for (uint j=0; j!=N_RNA; ++j) {
	profile_[i][j] += iupac_weight[r][j]*w;
      }
    }
  }
  n_seqs_ += w;
}

// public
void
ProfileSequence::
add_sequence(const ProfileSequence& seq, value_type w /*=1.0*/)
{
  if (profile_.size() != seq.size())
    throw "alignment error";

  for (uint i=0; i!=profile_.size(); ++i) {
    for (uint j=0; j!=N_RNA; ++j) {
      profile_[i][j] += seq[i][j]*w;
    }
  }
  n_seqs_ += w;
}

