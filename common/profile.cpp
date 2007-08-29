// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "profile.h"

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
    ma_[n][i] = char2rna(seq[i]);
}
