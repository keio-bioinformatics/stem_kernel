// $Id:

#include "phmm.h"
#include <cassert>
#include <string>
#include <boost/shared_ptr.hpp>
#include "log_value.h"

using PHMMTS::LogValue;

template < class ScoreTable >
template < class Seq >
typename PairHMM<ScoreTable>::score_type
PairHMM<ScoreTable>::
forward(const Seq& x, const Seq& y, FwTable& fw) const
{
  uint st[] = {M, IX, IY};
  const uint N=3;
  fw.resize(boost::extents[N][x.size()+1][y.size()+1]);
  for (uint i=0; i!=x.size()+1; ++i)
    for (uint j=0; j!=y.size()+1; ++j)
      for (uint* s=st; s!=st+N; ++s) 
	fw[*s][i][j] = static_cast<score_type>(0.0);
  fw[M][0][0]  = static_cast<score_type>(1.0);

  for (uint i=1; i!=x.size()+1; ++i) {
    fw[M][i][0] = fw[IY][i][0] = static_cast<score_type>(0.0);
    for (uint* s=st; s!=st+N; ++s) {
      fw[IX][i][0] += fw[*s][i-1][0] * score_(*s, IX, x[i-1]);
    }
  }

  for (uint j=1; j!=y.size()+1; ++j) {
    fw[M][0][j] = fw[IX][0][j] = static_cast<score_type>(0.0);
    for (uint* s=st; s!=st+N; ++s) {
      fw[IY][0][j] += fw[*s][0][j-1] * score_(*s, IY, y[j-1]);
    }
  }

  for (uint i=1; i!=x.size()+1; ++i) {
    for (uint j=1; j!=y.size()+1; ++j) {
      for (uint* s=st; s!=st+N; ++s) {
	fw[M][i][j] += fw[*s][i-1][j-1] * score_(*s, M, x[i-1], y[j-1]);
	fw[IX][i][j] += fw[*s][i-1][j] * score_(*s, IX, x[i-1]);
	fw[IY][i][j] += fw[*s][i][j-1] * score_(*s, IY, y[j-1]);
      }
    }
  }

  return fw[M][x.size()][y.size()];
}

template < class ScoreTable >
template < class Seq >
typename PairHMM<ScoreTable>::score_type
PairHMM<ScoreTable>::
backward(const Seq& x, const Seq& y, BkTable& bk) const
{
  uint st[] = {M, IX, IY};
  uint N=3;
  bk.resize(boost::extents[N][x.size()+1][y.size()+1]);
  for (uint i=0; i!=x.size()+1; ++i)
    for (uint j=0; j!=y.size()+1; ++j)
      for (uint* s=st; s!=st+N; ++s) 
	bk[*s][i][j] = static_cast<score_type>(0.0);
  bk[M][x.size()][y.size()]  = static_cast<score_type>(1.0);

  for (uint i=x.size(); i!=0; --i) {
    for (uint j=y.size(); j!=0; --j) {
      for (uint* s=st; s!=st+N; ++s) {
	bk[*s][i-1][j-1] += bk[M][i][j] * score_(*s, M, x[i-1], y[j-1]);
	bk[*s][i-1][j] += bk[IX][i][j] * score_(*s, IX, x[i-1]);
	bk[*s][i][j-1] += bk[IY][i][j] * score_(*s, IY, y[j-1]);
      }
    }
  }

  for (uint j=y.size(); j!=0; --j) {
    bk[M][0][j] = bk[IX][0][j] = static_cast<score_type>(0.0);
    for (uint* s=st; s!=st+N; ++s) {
      bk[*s][0][j-1] += bk[IY][0][j] * score_(*s, IY, y[j-1]);
    }
  }

  for (uint i=x.size(); i!=0; --i) {
    bk[M][i][0] = bk[IY][i][0] = static_cast<score_type>(0.0);
    for (uint* s=st; s!=st+N; ++s) {
      bk[*s][i-1][0] += bk[IX][i][0] * score_(*s, IX, x[i-1]);
    }
  }

  return bk[M][0][0];
}

template < class ScoreTable >
template < class Seq >
void
PairHMM<ScoreTable>::
forward_backward(const Seq& x, const Seq& y,
		 const FwTable& fw, const BkTable& bk, FBTable& fb) const
{
  uint st[] = {M, IX, IY};
  uint N=3;
  fb.resize(boost::extents[N][x.size()+1][y.size()+1]);
  score_type w = fw[M][x.size()][y.size()];

  for (uint i=0; i!=x.size()+1; ++i) {
    for (uint j=0; j!=y.size()+1; ++j) {
      for (uint* s=st; s!=st+N; ++s) {
	fb[*s][i][j] = fw[*s][i][j] * bk[*s][i][j] / w;
      }
    }
  }
}

template < class ScoreTable >
double
PairHMM<ScoreTable>::
forward(const FBTable& fb, TrTable& tr) const
{
  uint st[] = {M, IX, IY};
  const uint N=3;
  uint x_size=fb[0].size();
  uint y_size=fb[0][0].size();
  FBTable fw(boost::extents[N][x_size][y_size]);
  tr.resize(boost::extents[N][x_size][y_size]);

  for (uint i=0; i!=x_size; ++i)
    for (uint j=0; j!=y_size; ++j)
      for (uint* s=st; s!=st+N; ++s) {
	fw[*s][i][j] = 0.0;
	tr[*s][i][j] = static_cast<uint>(-1);
      }
  for (uint* s=st; s!=st+N; ++s)
    fw[*s][0][0] = fb[*s][0][0];

  for (uint i=1; i!=x_size; ++i) {
    fw[M][i][0] = fw[IY][i][0] = static_cast<score_type>(0.0);
    for (uint* s=st; s!=st+N; ++s) {
      double v=fw[*s][i-1][0]+fb[IX][i][0];
      if (tr[IX][i][0]==static_cast<uint>(-1) || fw[IX][i][0]<v) {
	fw[IX][i][0] = v;
	tr[IX][i][0] = *s;
      }
    }
  }

  for (uint j=1; j!=y_size; ++j) {
    fw[M][0][j] = fw[IX][0][j] = static_cast<score_type>(0.0);
    for (uint* s=st; s!=st+N; ++s) {
      double v=fw[*s][0][j-1]+fb[IY][0][j];
      if (tr[IY][0][j]==static_cast<uint>(-1) || fw[IY][0][j]<v) {
	fw[IY][0][j] = v;
	tr[IY][0][j] = *s;
      }
    }
  }

  for (uint i=1; i!=x_size; ++i) {
    for (uint j=1; j!=y_size; ++j) {
      for (uint* s=st; s!=st+N; ++s) {
	double v=fw[*s][i-1][j-1]+fb[M][i][j];
	if (tr[M][i][j]==static_cast<uint>(-1) || fw[M][i][j]<v) {
	  fw[M][i][j] = v;
	  tr[M][i][j] = *s;
	}

	v=fw[*s][i-1][j]+fb[IX][i][j];
	if (tr[IX][i][j]==static_cast<uint>(-1) || fw[IX][i][j]<v) {
	  fw[IX][i][j] = v;
	  tr[IX][i][j] = *s;
	}

	v=fw[*s][i][j-1]+fb[IY][i][j];
	if (tr[IY][i][j]==static_cast<uint>(-1) || fw[IY][i][j]<v) {
	  fw[IY][i][j] = v;
	  tr[IY][i][j] = *s;
	}
      }
    }
  }

  return fw[M][x_size-1][y_size-1];
}

template < class ScoreTable >
void
PairHMM<ScoreTable>::
traceback(const FBTable& fb, const TrTable& tr, DPPath& res) const
{
  //uint st[] = {M, IX, IY};
  //const uint N=N;
  uint x_size=fb[0].size();
  uint y_size=fb[0][0].size();

  Pos pos(M, x_size-1, y_size-1);
  res.push_front(pos);
  while (pos.x!=0 && pos.y!=0) {
    switch (pos.s) {
    case M:
      pos=Pos(tr[pos.s][pos.x][pos.y], pos.x-1, pos.y-1);
      break;
    case IX:
      pos=Pos(tr[pos.s][pos.x][pos.y], pos.x-1, pos.y);
      break;
    case IY:
      pos=Pos(tr[pos.s][pos.x][pos.y], pos.x, pos.y-1);
      break;
    default:
      assert(!"no such state");
      break;
    }
    res.push_front(pos);
  }
}

static
uint
char2rna(char x)
{
  switch(x) {
  case 'a': case 'A': return 0;
  case 'c': case 'C': return 1;
  case 'g': case 'G': return 2;
  case 'u': case 'U': return 3;    
  }
  assert(!"unknown nucleotide");
  return 4;
}

static
double
ribosum_trans[3][3] = {
  {   0.0,  -5.0,  -5.0 }, // M
  { -10.0,  -5.0, -15.0 }, // IX
  { -10.0,  -5.0, -15.0 }, // IY
};

static
double
ribosum_emit[4][4] = {
  {  2.22, -1.86, -1.46, -1.39 }, // A
  { -1.86,  1.16, -2.48, -1.05 }, // C
  { -1.46, -2.48,  1.03, -1.74 }, // G
  { -1.39, -1.05, -1.74,  1.65 }, // U
};

Ribosum::score_type
Ribosum::
operator()(uint s, uint t, char x) const
{
  return trans_[s][t];
}

Ribosum::score_type
Ribosum::
operator()(uint s, uint t, uint x) const
{
  return trans_[s][t];
}

Ribosum::score_type
Ribosum::
operator()(uint s, uint t, char x, char y) const
{
  return (*this)(s,t,char2rna(x),char2rna(y));
}

Ribosum::score_type
Ribosum::
operator()(uint s, uint t, uint x, uint y) const
{
  assert(t==0);
  return trans_[s][t]*emit_[x][y];
}

Ribosum::
Ribosum()
{
  for (uint i=0; i!=3; ++i) {
    for (uint j=0; j!=3; ++j) {
      trans_[i][j] =
	static_cast<score_type>(score_type::ExpOf(ribosum_trans[i][j]));
    }
  }
  for (uint i=0; i!=4; ++i) {
    for (uint j=0; j!=4; ++j) {
      emit_[i][j] =
	static_cast<score_type>(score_type::ExpOf(ribosum_emit[i][j]));
    }
  }
}

//static
boost::shared_ptr<Ribosum>
Ribosum::
get_instance()
{
  if (instance_==NULL)
    instance_ = boost::shared_ptr<Ribosum>(new Ribosum);
  return instance_;
}

//static
boost::shared_ptr<Ribosum> Ribosum::instance_;

// instantiation
template class PairHMM<Ribosum>;

template
PairHMM<Ribosum>::score_type
PairHMM<Ribosum>::
forward(const std::string& x, const std::string& y,
	FwTable& fw) const;

template
PairHMM<Ribosum>::score_type
PairHMM<Ribosum>::
backward(const std::string& x, const std::string& y,
	 BkTable& bk) const;

template
void
PairHMM<Ribosum>::
forward_backward(const std::string& x, const std::string& y,
		 const FwTable& fw, const BkTable& bk, FBTable& fb) const;

#if 0
#include <iostream>

int
main()
{
  boost::shared_ptr<Ribosum> st=Ribosum::get_instance();
  PairHMM<Ribosum> phmm(*st);
  PairHMM<Ribosum>::DPTable fw, bk, fb;
  std::string x, y;
  while (1) {
    std::getline(std::cin, x);
    if (x.empty()) break;
    std::getline(std::cin, y);
    if (y.empty()) break;
    std::cout << phmm.forward(x, y, fw) << " "
	      << phmm.backward(x, y, bk) << std::endl;;
    phmm.forward_backward(x, y, fw, bk, fb);
    for (uint i=0; i!=x_size; ++i) {
      for (uint j=0; j!=y_size; ++j) {
	std::cout << fb[0][i][j] << " ";
      }
      std::cout << std::endl;
    }
  }
  return 0;
}
#endif
