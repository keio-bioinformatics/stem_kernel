#include <iostream>
#include <list>
#include <string>
#include <cmath>
#include <boost/multi_array.hpp>
#include "fa.h"

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
  std::list<Fasta> fa;
  std::cout << "load " 
	    << Fasta::load(fa, argv[1])
	    << " seqs" << std::endl;
  std::list<Fasta>::const_iterator x,y;
  float min_v=100, max_v=0;
  float sum_v=0.0, sum=0.0;
  uint n=0;
  for (x=fa.begin(); x!=fa.end(); ++x) {
    for (y=x, ++y; y!=fa.end(); ++y) {
      float v=dp_match(x->seq(), y->seq())*2;
      v /= x->seq().size()+y->seq().size();
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
