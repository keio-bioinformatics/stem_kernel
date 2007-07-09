// $Id$

#ifndef __INC_DATA_H__
#define __INC_DATA_H__

#include <vector>
#include <boost/spirit/iterator/file_iterator.hpp>
#include "dag.h"
#include "maf.h"

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

  Data(const IS& seq);

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
void set_folding_method(uint method);

void set_bp_threshold(float th);

void set_window_size(uint win_sz, uint pair_sz);
  
template < class Data >
bool
load_examples(const std::string& label,
	      const char* filename,
	      std::vector< std::pair<std::string, Data> >& ex,
	      uint n_th=1, uint th_no=0);

template < class Data >
bool
load_examples(const std::string& label,
	      const char* filename,
	      std::vector< std::pair<std::string, Data> >& ex,
	      const std::vector< uint >& sv_index,
	      uint n_th=1, uint th_no=0);

template < class Data >
struct MakeData
{
  MakeData(const char* f) { }
  bool operator()(Data& data)
  {
    return false;
  }
};

#if 0
typedef Data<std::string> SData;
typedef Data< MASequence<std::string> > MData;
#else
typedef Data<RNASequence,std::string> SData;
//typedef Data< MASequence<RNASequence>,MASequence<std::string> > MData;
typedef Data< std::vector<Col>,MASequence<std::string> > MData;
#endif

template < >
struct MakeData<SData>
{
  std::string filename;
  uint type;
  boost::spirit::file_iterator<> fi;
  
  MakeData(const char* f);
  bool operator()(SData& data, bool skip=false);
};

template < >
struct MakeData<MData>
{
  std::string filename;
  uint type;
  boost::spirit::file_iterator<> fi;
  MAF maf;
  bool available;

  MakeData(const char* f);
  bool operator()(MData& data, bool skip=false);
};

#endif	// __INC_DATA_H__

// Local Variables:
// mode: C++
// End:
