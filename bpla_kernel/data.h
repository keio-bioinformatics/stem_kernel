// $Id$

#ifndef __INC_DATA_H__
#define __INC_DATA_H__

#include <vector>
#include <boost/spirit/iterator/file_iterator.hpp>
#include "../common/bpmatrix.h" // for definitions folding methods
#include "../common/profile.h"
namespace Vienna {
  extern "C" {
#include <ViennaRNA/utils.h>
  };
};

template < class S, class IS >
struct Data
{
  typedef S Seq;
  
  Seq seq;
  std::vector<float> p_left;
  std::vector<float> p_right;
  std::vector<float> p_unpair;

  Data() : seq(), p_left(), p_right(), p_unpair() { }

  Data(const IS& seq, float pf_scale, const BPMatrix::Options& opts);

  Data(const IS& seq);

  Data(const Data& x)
    : seq(x.seq), p_left(x.p_left), p_right(x.p_right), p_unpair(x.p_unpair)
  {
  }

  uint size() const { return seq.size(); }
};

typedef Data< ProfileSequence, std::string> SData;
typedef Data< ProfileSequence, std::list<std::string> > MData;

template < class D >
class DataLoader
{
public:
  typedef D Data;

public:
  DataLoader(const char* filename,
	     const BPMatrix::Options& bp_opts, bool use_bp)
  {
  }
  ~DataLoader();

  Data* get() { return NULL; }

private:
};

template < >
class DataLoader<SData>
{
public:
  typedef SData Data;

public:
  DataLoader(const char* filename,
	     const BPMatrix::Options& bp_opts, bool use_bp);
  DataLoader(const char* filename, const char* pf_scales,
	     const BPMatrix::Options& bp_opts, bool use_bp);
  ~DataLoader();
  Data* get();

private:
  const BPMatrix::Options& bp_opts_;
  bool use_bp_;
  std::string filename_;
  uint type_;
  boost::spirit::file_iterator<> fi_;
  std::ifstream* pf_is_;
};

template < >
class DataLoader<MData>
{
public:
  typedef MData Data;

public:
  DataLoader(const char* filename,
	     const BPMatrix::Options& bp_opts, bool use_bp);
  DataLoader(const char* filename, const char* pf_scapes,
	     const BPMatrix::Options& bp_opts, bool use_bp);
  ~DataLoader();
  Data* get();

private:
  const BPMatrix::Options& bp_opts_;
  bool use_bp_;
  std::string filename_;
  uint type_;
  boost::spirit::file_iterator<> fi_;
  std::ifstream* pf_is_;
};

template < class LD >
class DataLoaderFactory
{
public:
  typedef LD Loader;
  typedef typename Loader::Data Data;

public:
  DataLoaderFactory(const BPMatrix::Options& bp_opts)
    : bp_opts_(bp_opts), use_bp_(true)
  {
    Vienna::init_rand();
  }

  DataLoaderFactory() : use_bp_(false)
  {
  }
  
  Loader* get_loader(const char* filename) const
  {
    return new Loader(filename, bp_opts_, use_bp_);
  }

  Loader* get_loader(const char* filename, const char* pf_scales) const
  {
    return new Loader(filename, pf_scales, bp_opts_, use_bp_);
  }

private:
  BPMatrix::Options bp_opts_;
  bool use_bp_;
};

template < >
class DataLoaderFactory< DataLoader<SData> >
{
public:
  typedef DataLoader<SData> Loader;
  typedef Loader::Data Data;

public:
  DataLoaderFactory(const BPMatrix::Options& bp_opts)
    : bp_opts_(bp_opts),  use_bp_(true)
  {
    Vienna::init_rand();
  }

  DataLoaderFactory() : use_bp_(false)
  {
  }
  
  Loader* get_loader(const char* filename) const;

private:
  BPMatrix::Options bp_opts_;
  bool use_bp_;
};

#endif	// __INC_DATA_H__

// Local Variables:
// mode: C++
// End:
