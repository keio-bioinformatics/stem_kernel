// $Id$

#ifndef __INC_DATA_H__
#define __INC_DATA_H__

#include <vector>
#include <boost/version.hpp>
#if BOOST_VERSION >= 103800
#include <boost/spirit/include/classic.hpp>
#else
#include <boost/spirit.hpp>
#endif
#include "dag.h"
#include "../common/bpmatrix.h" // for definitions folding methods
#include "../common/profile.h"
namespace Vienna {
  extern "C" {
#include <ViennaRNA/utils.h>
  };
};

#ifndef BOOST_SPIRIT_CLASSIC_NS
#define BOOST_SPIRIT_CLASSIC_NS boost::spirit
#endif

template < class S, class IS, class N = DAG::Node<DAG::Edge> >
struct Data
{
  typedef S Seq;
  typedef N Node;
  typedef typename Node::Edge Edge;
  
  std::vector<Node> tree;
  Seq seq;
  std::vector<uint> root;
  std::vector<uint> max_pa;
  std::vector<float> weight;

  Data() : tree(), seq(), root(), max_pa(), weight() { }

  Data(const IS& seq, float th, float pf_scale, const BPMatrix::Options& opts);

  Data(const IS& seq);

  Data(const Data& x)
    : tree(x.tree), seq(x.seq), root(x.root),
      max_pa(x.max_pa), weight(x.weight)
  {
  }

  uint used_memory_size() const;
  uint max_node_size() const;
};

typedef Data< ProfileSequence, std::list<std::string> > MData;

template < class D >
class DataLoader
{
public:
  typedef D Data;

public:
  DataLoader(const char* filename, float th,
	     const BPMatrix::Options& bp_opts, bool use_bp)
  {
  }

  DataLoader(const char* filename, float th, const char* pf_scales,
	     const BPMatrix::Options& bp_opts, bool use_bp)
  {
  }

  Data* get() { return NULL; }

private:
};

template < >
class DataLoader<MData>
{
public:
  typedef MData Data;

public:
  DataLoader(const char* filename, float th,
	     const BPMatrix::Options& bp_opts, bool use_bp);
  DataLoader(const char* filename, float th, const char* pf_scales,
	     const BPMatrix::Options& bp_opts, bool use_bp);
  ~DataLoader();

  Data* get();

private:
  float th_;
  const BPMatrix::Options& bp_opts_;
  bool use_bp_;
  std::string filename_;
  uint type_;
  BOOST_SPIRIT_CLASSIC_NS::file_iterator<> fi_;
  std::ifstream* pf_is_;
};

template < class LD >
class DataLoaderFactory
{
public:
  typedef LD Loader;
  typedef typename Loader::Data Data;

public:
  DataLoaderFactory(float th, const BPMatrix::Options& bp_opts)
    : th_(th), bp_opts_(bp_opts), use_bp_(true)
  {
    Vienna::init_rand();
  }

  DataLoaderFactory() : th_(0.0), use_bp_(false)
  {
  }
  
  Loader* get_loader(const char* filename) const
  {
    return new Loader(filename, th_, bp_opts_, use_bp_);
  }

  Loader* get_loader(const char* filename, const char* pf_scales) const
  {
    return new Loader(filename, th_, pf_scales, bp_opts_, use_bp_);
  }

private:
  float th_;
  BPMatrix::Options bp_opts_;
  bool use_bp_;
};

#endif	// __INC_DATA_H__

// Local Variables:
// mode: C++
// End:
