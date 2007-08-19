// $Id$

#ifndef __INC_DATA_H__
#define __INC_DATA_H__

#include <vector>
#include <boost/spirit/iterator/file_iterator.hpp>
#include "dag.h"

template < class S, class IS, class N = DAG::Node<DAG::Edge> >
struct Data {
  typedef S Seq;
  typedef N Node;
  typedef typename Node::Edge Edge;
  
  std::vector<Node> tree;
  Seq seq;
  std::vector<uint> root;
  std::vector<float> weight;
  std::vector<uint> max_pa;

  Data() : tree(), seq(), root(), weight(), max_pa() { }

  Data(const IS& seq, uint method, float th);

  Data(const Data& x)
    : tree(x.tree), seq(x.seq), root(x.root),
      weight(x.weight), max_pa(x.max_pa)
  {
  }
};

struct Col
{
  std::vector<unsigned char> cnt_;

  Col() : cnt_(N_RNA, 0) { }

  unsigned char cnt(uint i) const { return cnt_[i]; }
  unsigned char& cnt(uint i) { return cnt_[i]; }
};

enum { ALIFOLD, FOLD, LFOLD, NO_BPMATRIX };	// available methods
  

#if 0
typedef Data<std::string> SData;
typedef Data< MASequence<std::string> > MData;
#else
typedef Data<RNASequence,std::string> SData;
//typedef Data< MASequence<RNASequence>,MASequence<std::string> > MData;
typedef Data< std::vector<Col>,MASequence<std::string> > MData;
#endif

template < class D >
class DataLoader
{
public:
  typedef D Data;

public:
  DataLoader(const char* filename, uint method, float th) { };

  Data* get() { return NULL; }

private:
  uint method_;
  float th_;
};

template < >
class DataLoader<SData>
{
public:
  typedef SData Data;

public:
  DataLoader(const char* filename, uint method, float th);
  Data* get();

private:
  uint method_;
  float th_;
  std::string filename_;
  uint type_;
  boost::spirit::file_iterator<> fi_;
};

template < >
class DataLoader<MData>
{
public:
  typedef MData Data;

public:
  DataLoader(const char* filename, uint method, float th);

  Data* get();

private:
  uint method_;
  float th_;
  std::string filename_;
  uint type_;
  boost::spirit::file_iterator<> fi_;
};

template < class LD >
class DataLoaderFactory
{
public:
  typedef LD Loader;
  typedef typename Loader::Data Data;

public:
  DataLoaderFactory(uint method, float th)
    : method_(method), th_(th)
  {
  }
  
  Loader* get_loader(const char* filename) const
  {
    return new Loader(filename, method_, th_);
  }

private:
  uint method_;
  float th_;
};

template < >
class DataLoaderFactory< DataLoader<SData> >
{
public:
  typedef DataLoader<SData> Loader;
  typedef Loader::Data Data;

public:
  DataLoaderFactory(uint method, float th)
    : method_(method), th_(th)
  {
  }
  
  Loader* get_loader(const char* filename) const;

private:
  uint method_;
  float th_;
};



#endif	// __INC_DATA_H__

// Local Variables:
// mode: C++
// End:
