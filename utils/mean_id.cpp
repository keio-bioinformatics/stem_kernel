#include <iostream>
#include <list>
#include <string>
#include <boost/multi_array.hpp>
#include <boost/spirit.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost::spirit;

struct fa_parser : public grammar< fa_parser >
{
  fa_parser(std::list<std::string>& x) : fa(x) {}

  std::list<std::string>& fa;

  struct append_seq
  {
    std::list<std::string>& fa;

    append_seq(std::list<std::string>& x) : fa(x) { }

    template < class Ite >
    void operator()(Ite i1, Ite i2) const
    {
      std::string s(i1, i2);
      std::string::const_iterator x;
      boost::erase_all(s, "-");
      for (x=s.begin(); x!=s.end(); ++x)
	fa.back().push_back(*x);
    }
  };

  struct add_seq
  {
    std::list<std::string>& fa;

    add_seq(std::list<std::string>& x) : fa(x) { }

    template < class Ite >
    void operator()(Ite i1, Ite i2) const
    {
      fa.push_back(std::string());
    }
  };

  template <class ScannerT>
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t fa;
    rule_t head;
    rule_t seq_l;
    rule_t seq;
    rule_t paren;

    definition(const fa_parser& self)
    {
      fa = *(head >> seq >> !paren) >> end_p;
      head = ch_p('>') >> (*(blank_p | graph_p))[add_seq(self.fa)] >> eol_p;
      seq_l = *(graph_p - chset<>("<().") - eol_p);
      seq = +(seq_l[append_seq(self.fa)] >> eol_p);
      paren = +chset<>("<().") >> eol_p;
    }

    const rule_t& start() const { return fa; }
  };
};

bool
load_fa(std::list<std::string>& fa, const char* filename)
{
  fa_parser parser(fa);
  file_iterator<> fi(filename);
  if (!fi) {
    std::cout << "unable to open file: " << filename << std::endl;
    return false;
  }
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  return info.full;
}

int
dp_match(const std::string& x, const std::string& y)
{
  typedef boost::multi_array<int,2> dp_type;
  dp_type dp(boost::extents[x.size()+1][y.size()+1]);

  dp[0][0]=0;
  for (uint i=1; i!=x.size()+1; ++i) dp[i][0]=0;
  for (uint j=1; j!=y.size()+1; ++j) dp[0][j]=0;

  for (uint j=1; j!=y.size()+1; ++j) {
    for (uint i=1; i!=x.size()+1; ++i) {
      if (x[i-1]==y[j-1])
	dp[i][j]=dp[i-1][j-1]+1;
      else
	dp[i][j]=dp[i-1][j-1]-1;
      if (dp[i][j]<dp[i-1][j])
	dp[i][j]=dp[i-1][j];
      if (dp[i][j]<dp[i][j-1])
	dp[i][j]=dp[i][j-1];
    }
  }
  
  return dp[x.size()][y.size()];
}

int
main(int argc, char** argv)
{
  std::list<std::string> fa;
  load_fa(fa, argv[1]);
  std::list<std::string>::const_iterator x,y;
  float min_v=100, max_v=0;
  float sum_v=0.0, sum=0.0;
  uint n=0;
  for (x=fa.begin(); x!=fa.end(); ++x) {
    for (y=x, ++y; y!=fa.end(); ++y) {
      float v=dp_match(*x,*y)*2;
      v /= x->size()+y->size();
      if (v<min_v) min_v=v;
      if (v>max_v) max_v=v;
      sum += v;
      sum_v += v*v;
      n++;
    }
  }
  std::cout << sum/n << ", "
	    << sqrt(sum_v/n - sum/n*sum/n) << ", "
	    << min_v << ", "
	    << max_v << std::endl;
  return 0;
}
