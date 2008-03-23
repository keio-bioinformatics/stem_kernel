// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cmath>
#include <vector>
#include <list>
#include <sstream>
#include "data.h"
#include "bpmatrix.h"
#include "../common/profile.h"
#include "../common/cyktable.h"
#include "../common/fa.h"
#include "../common/maf.h"
#include "../common/aln.h"

template < class IS >
static
void
fill_weight(const IS& s,
	    std::vector<float>& p_l,
	    std::vector<float>& p_r,
	    std::vector<float>& p_u,
	    const BPMatrix::Options& opts)
{
  BPMatrix bp(s, opts);
  assert(bp.size()==p_l.size());
  assert(bp.size()==p_r.size());
  assert(bp.size()==p_u.size());
  for (uint i=0; i<bp.size(); ++i) {
    p_l[i]=p_r[i]=0.0;
    for (uint j=0; j<i; ++j)
      p_l[i] += bp(j+1, i+1);
    for (uint j=i+1; j<bp.size(); ++j) 
      p_r[i] += bp(i+1, j+1);
    p_u[i] = 1.0-(p_l[i]+p_r[i]);
    if (p_u[i] < 0.0) p_u[i] = 0.0; // error around
  }
}

template < class S, class IS >
Data<S,IS>::
Data(const IS& s, const BPMatrix::Options& opts)
  : seq(s),
    p_left(s.begin()->size(), 0.0),
    p_right(s.begin()->size(), 0.0),
    p_unpair(s.begin()->size(), 1.0)
{
  fill_weight(s, p_left, p_right, p_unpair, opts);
}
 
template < class S, class IS >
Data<S,IS>::
Data(const IS& s)
  : seq(s), p_left(), p_right(), p_unpair()
{
}

enum { TP_UNKNOWN, TP_FA, TP_ALN, TP_MAF };

static
uint
check_filetype(const char* f)
{
#if 0  
  boost::spirit::file_iterator<> fi(f);
  std::string s;
  MASequence<std::string> ma;
  if (load_fa(s, fi)) return TP_FA;
  if (load_aln(ma, fi)) return TP_ALN;
  if (load_maf(ma, fi)) return TP_MAF;
  return TP_UNKNOWN;
#else
  std::ifstream in(f);
  std::string l;
  while (std::getline(in, l)) {
    if (l[0]=='>') return TP_FA;
    else if (l.compare(0, 7, "CLUSTAL")==0) return TP_ALN;
    else if (l.compare(0, 2, "a ")==0) return TP_MAF;
  }
  return TP_UNKNOWN;
#endif
}

#if 0
DataLoader<SData>::
DataLoader(const char* filename,
	   const BPMatrix::Options& bp_opts, bool use_bp)
  : bp_opts_(bp_opts),
    use_bp_(use_bp),
    filename_(filename),
    type_(check_filetype(filename)),
    fi_(filename)
{
  if (!fi_) {
    std::ostringstream os;
    os << filename_ << ": no such file";
    throw os.str().c_str();
  }
}

SData*
DataLoader<SData>::
get()
{
  std::string s;
  switch (type_) {
  case TP_FA:
    if (load_fa(s, fi_)) {
      if (use_bp_)
	return new SData(s, bp_opts_);
      else
	return new SData(s);
    } else {
      return NULL;
    }
    break;
  default:
    return NULL;
    break;
  }
  return NULL;
}
#endif

DataLoader<MData>::
DataLoader(const char* filename,
	   const BPMatrix::Options& bp_opts, bool use_bp)
  : bp_opts_(bp_opts),
    use_bp_(use_bp),
    filename_(filename),
    type_(check_filetype(filename)),
    fi_()
{
  switch (type_) {
  case TP_FA:
  case TP_ALN:
  case TP_MAF:
    fi_ = boost::spirit::file_iterator<>(filename_);
    break;
  default:
    break;
  }
  if (!fi_) {
    std::ostringstream os;
    os << filename_ << ": no such file";
    throw os.str().c_str();
    //return false;
  }
}

MData*
DataLoader<MData>::
get()
{
  bool ret=false;
  std::list<std::string> ma;
  switch (type_) {
  case TP_FA:
    ret=load_fa(ma, fi_);
    break;
  case TP_ALN:
    ret=load_aln(ma, fi_);
    break;
  case TP_MAF:
    ret=load_maf(ma, fi_);
    break;
  default:
    break;
  }

  if (ret) {
    std::list<std::string>::const_iterator x;
    uint l=ma.begin()->size();
    for (x=ma.begin(); x!=ma.end(); ++x) {
      if (l!=x->size()) throw "wrong alignment";
    }
  
    if (use_bp_)
      return new MData(ma, bp_opts_);
    else
      return new MData(ma);
  }
  return NULL;
}

#if 0
DataLoader<SData>*
DataLoaderFactory< DataLoader<SData> >::
get_loader(const char* filename) const
{
  switch (check_filetype(filename)) {
  case TP_FA:
    return new DataLoader<SData>(filename, bp_opts_, use_bp_);
    break;
  default:
    return NULL;
    break;
  }
  return NULL;
}
#endif

// instantiation

#include <string>
#include "bpmatrix.h"
#include "../common/rna.h"

#if 0
template
class Data<ProfileSequence, std::string>;
#endif

template
class Data<ProfileSequence, std::list<std::string> >;
