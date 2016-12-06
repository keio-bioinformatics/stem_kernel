// $Id$

#ifndef __INC_PROBLEM_H__
#define __INC_PROBLEM_H__

#include <iosfwd>
#include <vector>
#include <string>
#include <unordered_map>
#include <svm.h>
#include <boost/shared_array.hpp>

typedef unsigned int uint;

typedef boost::shared_array<struct svm_node> node_ary;

class Problem
{
public:
  Problem();
  Problem(const std::vector<node_ary>& x,
	  const std::vector<double>& y);
  Problem(const Problem& prob);
  ~Problem();
  Problem& operator=(const Problem& prob);

  const struct svm_problem* get() const { return &prob_; }
  uint size() const { return prob_.l; }
  const struct svm_node* x(uint i) const { return prob_.x[i]; }
  const std::vector<node_ary>& x() const { return x_; }
  const std::vector<double>& y() const { return y_; }

  std::ostream& print(std::ostream& os) const;
  void split(uint fold, uint n, Problem& tr, Problem& ts) const;
    

private:
  void set_problem();

private:
  struct svm_problem prob_;
  std::vector<node_ary> x_;
  std::vector<double> y_;
};

class Data
{
public:
  Data() : labels_(), vec_() { }

  Data(std::istream& is) : labels_(), vec_()
  {
    read(is);
  }

  void read(std::istream& is);
  void read(const std::string& f);
  bool add(std::istream& is);
  bool add(const std::string& f);
  Problem select(const std::unordered_map<std::string,int>& pn_map) const;
#ifdef HAVE_BOOST_XPRESSIVE_XPRESSIVE_HPP
  Problem select(const std::string& pos, const std::string& neg) const;
#endif
  Problem select(const std::vector<std::string>& pos,
		 const std::vector<std::string>& neg) const;
  Problem select() const;

  void normalize();

  std::ostream& print(std::ostream& os) const;

private:
  std::vector<std::string> labels_;
  std::vector<node_ary> vec_;
};

#endif	// __INC_PROBLEM_H__

// Local Variables:
// mode: C++
// End:
