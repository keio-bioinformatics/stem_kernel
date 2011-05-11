// $Id$
#ifndef __INC_PHMM_H__
#define __INC_PHMM_H__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <deque>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>
#include "log_value.h"

typedef unsigned int uint;

class Ribosum
{
public:
  typedef PHMMTS::LogValue<double> score_type;

public:
  static boost::shared_ptr<Ribosum> get_instance();
  score_type operator()(uint s, uint t, char x) const;
  score_type operator()(uint s, uint t, uint x) const;
  score_type operator()(uint s, uint t, char x, char y) const;
  score_type operator()(uint s, uint t, uint x, uint y) const;

private:
  score_type trans_[3][3];
  score_type emit_[4][4];
  static boost::shared_ptr<Ribosum> instance_;

private:  
  Ribosum();
};

template < class ScoreTable = Ribosum >
class PairHMM
{
public:
  enum { M=0, IX=1, IY=2, N_STATES=3 };
  struct Pos {
    uint s;
    uint x;
    uint y;
    Pos() {}
    Pos(uint s, uint x, uint y) : s(s), x(x), y(y) { }
  };
  typedef ScoreTable score_table;
  typedef typename ScoreTable::score_type score_type;
  typedef boost::multi_array<score_type, 3> FwTable;
  typedef boost::multi_array<score_type, 3> BkTable;
  typedef boost::multi_array<double, 3> FBTable;
  typedef boost::multi_array<uint, 3> TrTable;
  typedef std::deque<Pos> DPPath;
  
public:
  PairHMM(const ScoreTable& score = *ScoreTable::get_instance())
    : score_(score)
  {
  }

  template < class Seq >
  score_type forward(const Seq& x, const Seq& y, FwTable& fw) const;
  template < class Seq >
  score_type backward(const Seq& x, const Seq& y, BkTable& bk) const;
  template < class Seq >
  void forward_backward(const Seq& x, const Seq& y,
			const FwTable& fw, const BkTable& bk,
			FBTable& fb) const;
  template < class Seq >
  void forward_backward(const Seq& x, const Seq& y,
			FBTable& fb) const
  {
    FwTable fw;
    BkTable bk;
    forward(x, y, fw);
    backward(x, y, bk);
    forward_backward(x, y, fw, bk, fb);
  }

  double
  forward(const FBTable& fb, TrTable& tr) const;

  void
  traceback(const FBTable& fb, const TrTable& tr, DPPath& res) const;

  void
  map_path(const FBTable& fb, DPPath& res) const
  {
    TrTable tr;
    forward(fb, tr);
    traceback(fb, tr, res);
  }

private:
  const ScoreTable& score_;
};

#endif // __INC_PHMM_H__

// Local Variables:
// mode: C++
// End:
