// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "profile.h"

#if 0
static
bool
iupac_symbol[N_IUPAC][N_RNA] = {
  // A      C      G      T(U)
  { true,  false, false, false }, // A
  { false, true,  false, false }, // C
  { false, false, true,  false }, // G
  { false, false, false, true  }, // T(U)
  { false, false, false, false }, // GAP
  { true,  false, true,  false }, // R = A or G
  { false, true,  false, true  }, // Y = C or T(U)
  { true,  true,  false, false }, // M = A or C
  { false, false, true,  true  }, // K = G or T(U)
  { false, true,  true,  false }, // S = C or G
  { true,  false, false, true  }, // W = A or T(U)
  { false, true,  true,  true  }, // B = C or G or T(U)
  { true,  false, true,  true  }, // D = A or G or T(U)
  { true,  true,  false, true  }, // H = A or C or T(U)
  { true,  true,  true,  false }, // V = A or C or G
  { true,  true,  true,  true  }, // N = A or C or G or T(U)
};
#endif

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

// private
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
	profile_[i][j] = iupac_weight[r][j];
      }
    }
  }
  n_seqs_ += w;
}


// public
void
BPProfileMaker::
make_profile(value_type p, uint i, uint j,
	     std::map<bp_t, value_type>& v) const
{
  for (uint n=0; n!=ma_[0].size(); ++n) {
    for (uint a=0; a!=N_RNA; ++a) {
      float aw=iupac_weight[ma_[i][n]][a];
      if (aw==0.0) continue;
      for (uint b=0; b!=N_RNA; ++b) {
	float bw=iupac_weight[ma_[j][n]][b];
	if (bw==0.0) continue;
	bp_t k = std::make_pair(a,b);
	std::map<bp_t, value_type>::iterator x = v.find(k);
	if (x==v.end()) {
	  v.insert(std::make_pair(k, p*aw*bw));
	} else {
	  x->second += p*aw*bw;
	}
      }
    }
  }
}

// public
void
BPProfileMaker::
make_profile(const std::vector<value_type>& p, uint i, uint j,
	     std::map<bp_t, value_type>& v) const
{
  if (p.size()!=ma_[0].size()) throw "size error";

  for (uint n=0; n!=ma_[0].size(); ++n) {
    for (uint a=0; a!=N_RNA; ++a) {
      float aw=iupac_weight[ma_[i][n]][a];
      if (aw==0.0) continue;
      for (uint b=0; b!=N_RNA; ++b) {
	float bw=iupac_weight[ma_[j][n]][b];
	if (bw==0.0) continue;
	bp_t k = std::make_pair(a,b);
	std::map<bp_t, value_type>::iterator x = v.find(k);
	if (x==v.end()) {
	  v.insert(std::make_pair(k, p[n]*aw*bw));
	} else {
	  x->second += p[n]*aw*bw;
	}
      }
    }
  }
}

// private
void
BPProfileMaker::
add_sequence(const std::string& seq, uint n)
{
  if (ma_.size() != seq.size())
    throw "alignment error";

  for (uint i=0; i!=seq.size(); ++i)
    ma_[i][n] = char2rna(seq[i]);
}
