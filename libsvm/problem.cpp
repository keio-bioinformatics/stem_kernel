// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "problem.h"
#include <deque>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstring>
#include <cerrno>
#ifdef HAVE_BOOST_IOSTREAMS
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef HAVE_BOOST_XPRESSIVE_XPRESSIVE_HPP
#include <boost/xpressive/xpressive.hpp>
#endif

#ifdef HAVE_BOOST_IOSTREAMS
namespace io = boost::iostreams;
#endif

Problem::
Problem() : x_(), y_()
{
  prob_.x = NULL;
  prob_.y = NULL;
  prob_.l = 0;
}

Problem::
Problem(const std::vector<node_ary>& x,
	const std::vector<double>& y)
  : x_(x), y_(y)
{
  prob_.x = NULL;
  prob_.y = NULL;
  prob_.l = 0;
  set_problem();
}

Problem::
Problem(const Problem& prob)
  : x_(prob.x_), y_(prob.y_)
{
  prob_.x = NULL;
  prob_.y = NULL;
  prob_.l = 0;
  set_problem();
}

Problem::
~Problem()
{
  if (prob_.x!=NULL) delete[] prob_.x;
}

Problem&
Problem::
operator=(const Problem& prob)
{
  if (this!=&prob) {
    x_ = prob.x_;
    y_ = prob.y_;
    set_problem();
  }
  return *this;
}

std::ostream&
Problem::
print(std::ostream& os) const
{
  for (int i=0; i!=prob_.l; ++i) {
    os << prob_.y[i] << " ";
    for (int j=0; prob_.x[i][j].index!=-1; ++j) {
      os << prob_.x[i][j].index << ":"
	 << prob_.x[i][j].value << " ";
    }
    os << std::endl;
  }
  return os;
}

void
Problem::
split(uint fold, uint n, Problem& tr, Problem& ts) const
{
  uint tr_n=0, ts_n=0;
  for (uint i=0; i!=x_.size(); ++i) {
    if (i%fold==n) { ts_n++; } else { tr_n++; }
  }

  std::vector<node_ary> tr_x(tr_n);
  std::vector<double> tr_y(tr_n);
  std::vector<node_ary> ts_x(ts_n);
  std::vector<double> ts_y(ts_n);

  tr_n=0; ts_n=0;
  for (uint i=0; i!=x_.size(); ++i) {
    if (i%fold==n) {
      ts_y[ts_n]=y_[i];
      ts_x[ts_n]=x_[i];
      ts_n++;
    } else {
      tr_y[tr_n]=y_[i];
      tr_x[tr_n]=x_[i];
      tr_n++;
    }
  }

  tr = Problem(tr_x, tr_y);
  ts = Problem(ts_x, ts_y);
}

// private
void
Problem::
set_problem()
{
  assert(x_.size()==y_.size());
  prob_.l = x_.size();
  prob_.y = &y_[0];
  if (prob_.x!=NULL) delete[] prob_.x;
  prob_.x = new struct svm_node*[prob_.l];
  for (uint i=0; i!=x_.size(); ++i)
    prob_.x[i] = x_[i].get();
}

#ifdef HAVE_BOOST_IOSTREAMS
static
void
open_dat(const std::string& f, std::ifstream& fis,
	 io::filtering_stream<io::input>& is)
{
  fis.open(f.c_str());
  if (!fis) {
    std::string e;
    e += strerror(errno);
    e += ": ";
    e += f;
    throw e.c_str();
  }
  if (f.rfind(".gz")+3==f.size())
    is.push(io::gzip_decompressor());
  else if (f.rfind(".bz2")+4==f.size())
    is.push(io::bzip2_decompressor());
  is.push(fis);
}
#else
static
void
open_dat(const std::string& f, std::ifstream& fs)
{
  is.open(f.c_str());
  if (!is) {
    std::string e;
    e += strerror(errno);
    e += ": ";
    e += f;
    throw e.c_str();
  }
}
#endif

void
Data::
read(const std::string& f)
{
#ifdef HAVE_BOOST_IOSTREAMS
  std::ifstream fis;
  io::filtering_stream<io::input> is;
  open_dat(f, fis, is);
#else
  std::ifstream is;
  open_dat(f, is);
#endif
  read(is);
}

void
Data::
read(std::istream& is)
{
  std::deque<std::string> labels;
  std::deque<node_ary> vec;
  std::string l;
  while (std::getline(is, l)) {
    std::istringstream ss(l);
    std::string v;
    ss >> v;
    labels.push_back(v);
    
    uint n=0;
    while (ss >> v) n++;
    node_ary na=node_ary(new svm_node[n+1]);
    vec.push_back(na);

    ss.clear(); ss.seekg(0);
    ss >> v;
    n=0;
    while (ss >> v) {
      sscanf(v.c_str(), "%d:%lf",
	     &na[n].index, &na[n].value);
      n++;
    }
    na[n].index=-1;
  }

  labels_.resize(labels.size());
  std::copy(labels.begin(), labels.end(), labels_.begin());
  vec_.resize(vec.size());
  std::copy(vec.begin(), vec.end(), vec_.begin());
}

bool
Data::
add(const std::string& f)
{
#ifdef HAVE_BOOST_IOSTREAMS
  std::ifstream fis;
  io::filtering_stream<io::input> is;
  open_dat(f, fis, is);
#else
  std::ifstream is;
  open_dat(f, is);
#endif
  return add(is);
}

bool
Data::
add(std::istream& is)
{
  std::string l;
  uint i=0;
  while (std::getline(is, l)) {
    std::istringstream ss(l);
    std::string v;
    ss >> v;
    if (v!=labels_[i]) return false;
    
    uint n=0;
    while (ss >> v) {
      int idx;
      double val;
      sscanf(v.c_str(), "%d:%lf", &idx, &val);
      if (vec_[i][n].index!=idx) return false;
      if (vec_[i][n].index==0) {
	if (vec_[i][n].value!=val) return false;
      } else {
	vec_[i][n].value += val;
      }
      n++;
    }
    if (vec_[i][n].index!=-1) return false;

    i++;
  }

  return true;
}

Problem
Data::
select() const
{
  std::vector<double> y(vec_.size(), 0);
  return Problem(vec_, y);
}

Problem
Data::
select(const std::unordered_map<std::string,int>& pn_map) const
{
  uint n=0;
  for (uint i=0; i!=labels_.size(); ++i) {
    std::unordered_map<std::string,int>::const_iterator f;
    f=pn_map.find(labels_[i]);
    if (f!=pn_map.end() && f->second!=0) {
      n++;
    }
  }
  if (n==0) return Problem();

  std::vector<node_ary> x(n);
  std::vector<double> y(n);
  n=0;
  for (uint i=0; i!=labels_.size(); ++i) {
    std::unordered_map<std::string,int>::const_iterator f;
    f=pn_map.find(labels_[i]);
    if (f!=pn_map.end()) {
      if (f->second>0) {
	y[n] = +1;
	x[n] = vec_[i];
	n++;
      } else if (f->second<0) {
	y[n] = -1;
	x[n] = vec_[i];
	n++;
      }
    }
  }
  return Problem(x, y);
}

#ifdef HAVE_BOOST_XPRESSIVE_XPRESSIVE_HPP
Problem
Data::
select(const std::string& pos, const std::string& neg) const
{
  using namespace boost::xpressive;

  sregex pos_reg = sregex::compile(pos);
  sregex neg_reg = sregex::compile(neg);
  sregex pos_neg_reg = pos_reg | neg_reg;

  uint n=0;
  for (uint i=0; i!=labels_.size(); ++i) {
    if (regex_match(labels_[i], pos_neg_reg)) {
      n++;
    }
  }
  if (n==0) return Problem();

  std::vector<node_ary> x(n);
  std::vector<double> y(n);
  n=0;
  for (uint i=0; i!=labels_.size(); ++i) {
    if (regex_match(labels_[i], pos_reg)) {
      y[n] = +1;
      x[n] = vec_[i];
      n++;
    } else if (regex_match(labels_[i], neg_reg)) {
      y[n] = -1;
      x[n] = vec_[i];
      n++;
    }
  }
  return Problem(x, y);
}
#endif

Problem
Data::
select(const std::vector<std::string>& pos,
       const std::vector<std::string>& neg) const
{
  if (pos.empty() && neg.empty())
    return select();
  
  std::unordered_map<std::string,int> pn_map;
  for (uint i=0; i!=pos.size(); ++i)
    pn_map[pos[i]]=+1;
  for (uint i=0; i!=neg.size(); ++i)
    pn_map[neg[i]]=-1;
  return select(pn_map);
}

void
Data::
normalize()
{
  uint sz = labels_.size();
  std::vector<double> diag(sz);
  std::vector<node_ary>::iterator v;
  for (v=vec_.begin(); v!=vec_.end(); ++v) {
    assert((*v)[0].index==0);
    int index=(*v)[0].value;
    for (uint i=1; (*v)[i].index!=-1; ++i) {
      if (index==(*v)[i].index) {
	diag[index-1] = sqrt((*v)[i].value);
	break;
      }
    }
  }

  for (v=vec_.begin(); v!=vec_.end(); ++v) {
    assert((*v)[0].index==0);
    uint index=(*v)[0].value;
    for (uint i=1; (*v)[i].index!=-1; ++i) {
      (*v)[i].value /= diag[index-1] * diag[(*v)[i].index-1];
    }
  }
}

std::ostream&
Data::
print(std::ostream& os) const
{
  for (uint i=0; i!=labels_.size(); ++i) {
    os << labels_[i] << " ";
    for (uint j=0; vec_[i][j].index!=-1; ++j) {
      os << vec_[i][j].index << ":"
	 << vec_[i][j].value << " ";
    }
    os << std::endl;
  }
  return os;
}

