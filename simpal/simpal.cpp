// $Id:$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "../common/pf_wrapper.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <ext/hash_map>
#include <cmath>
#include <boost/algorithm/string.hpp>


namespace __gnu_cxx {
  template <>
  struct hash<std::string> : public hash<const char*>
  {
    size_t operator()(const std::string& s) const
    {
      return hash<const char*>::operator()(s.c_str());
    }
  };
};


class KernelFunc;

class Pals
{
  friend class KernelFunc;
public:
  Pals() : pals_() { }  
  Pals(const std::string& seq, uint seed_length, uint min_loop, uint max_dist)
    : pals_()
  {
    make_pal_map(seq, seed_length, min_loop, max_dist);
  }

  Pals(const Pals& pals)
    : pals_(pals.pals_)
  {
  }

  Pals& operator=(const Pals& pals)
  {
    if (this!=&pals) {
      pals_=pals.pals_;
    }
    return *this;
  }

private:
  struct Pal 
  {
    Pal(const std::string& key, uint dist)
      : key_(key), dist_(dist)
    {
    }

    bool operator<(const Pal& p) const
    {
      return key_!=p.key_ ? key_<p.key_ : dist_<p.dist_;
    }
    
    std::string key_;
    int dist_;
  };

  struct RevComp
  {
    char operator()(char x) const
    {
      switch (x) {
      case 'a': case 'A':
	return 'u';
	break;
      case 'c': case 'C':
	return 'g';
	break;
      case 'g': case 'G':
	return 'c';
	break;
      case 'u': case 'U': case 't': case 'T':
	return 'a';
	break;
      default:
	break;
      }
      return 0;
    }
  };

  typedef std::map<Pal,float> PalMap;
  typedef __gnu_cxx::hash_multimap<std::string, int> Hash;

private:
  void
  make_pal_map(const std::string& seq,
	       uint seed_length, uint min_loop, uint max_dist)
  {
    Hash f,r;
    fill_hash(seq, f, seed_length, +1);
    std::string rev(seq.size(), ' ');
    std::transform(seq.rbegin(), seq.rend(), rev.begin(), RevComp());
    fill_hash(rev, r, seed_length, -1);
    find_pals(seq, f, r, seq.size(), seed_length, min_loop, max_dist);
  }

  void
  fill_hash(const std::string& seq, Hash& h, uint seed_length, int mult)
  {
    for (uint i=0; i!=seq.size()-seed_length; ++i) {
      std::string k(&seq[i], seed_length);
      h.insert(Hash::value_type(k, (i+1)*mult));
    }
  }

  void
  find_pals(const std::string& seq, const Hash& fwd, const Hash& rev,
	    uint seq_len, uint seed_length, uint min_loop, uint max_dist)
  {
   
    std::string st;
    PFWrapper pf(seq, true);
    //float f=pf.fold(st);
    pf.fold(st);
    
    //    std::cout << seq << std::endl;
    /*
    std::cout << seq << std::endl << st << std::endl << f << std::endl;
    
    for (uint j=0; j!=pf.size(); ++j){
      for (uint i=j;; --i){
	std::cout << pf(i,j) << " ";
	if (i==0) break;
      }
      std::cout << std::endl;
    }
    */
    
     
    Hash::const_iterator x;
    for (x=fwd.begin(); x!=fwd.end(); ++x) {
      std::pair<Hash::const_iterator,Hash::const_iterator> r;
      r=rev.equal_range(x->first);
      Hash::const_iterator y;
      for (y=r.first; y!=r.second; ++y) {
	uint d=seq_len-(x->second-y->second+2*seed_length-2);
	if (d>=min_loop && d<=max_dist) {
	  
	  //std::cout << x->second << " " << y->second << " " << d << std::endl;
	  float weight = 1;
	  
	  Pal p(x->first, d);
	  /*
	    std::cout << "seq_len\n";
	    std::cout << seq_len << std::endl;
	    std::cout << "x->second\ty->second\n";
	    std::cout << x->second << std::endl;
	    std::cout << y->second << std::endl;
	  */
	  
	  for(uint i=0; i < seed_length; ++i)
	    {	      
	      std::cout << "m\tn\n";
	      uint m = x->second + i ;
	      uint n = y->second + seq_len -i + 1;
	      
	      /*
		std::cout << "pf(m,n)\n";
		std::cout << pf(m, n) << std::endl;
	      */

	      /*
		std::cout << "1/(1-log(pf(m,n))) \n";
		std::cout << 1/(1-log(pf(m, n))) << std::endl;
	      */
	      
	      /*
		if(pf(m,n) >= 1.0e-1){
		
		weight = 1 * weight;
		}
	      */
	      /*
		else if(pf(m,n) < 1.0e-1)
	      {
	      
	      weight = 0 * weight;
	      
	      //weight = (1/(1-log(pf(m,n)))) * weight;
	      }
	      */
	      
	      weight = pf(m,n) * weight; 
	      
	      
	      
	    }
	  
	  //std::cout << "weight\n";
	  //std::cout << weight << std::endl;
	  
	  
	  
	  if (pals_.find(p)==pals_.end()) {
	    pals_.insert(std::make_pair(p,weight));
	  } else {
	    pals_[p]+= weight;
	  }
	}
      }
    }
  }
  

private:
  PalMap pals_;
};

class KernelFunc
{
public:
  typedef double value_type;

public:
  KernelFunc(int tolerance=1)
    : tolerance_(tolerance)
  {
  }
  
  double
  operator()(const Pals& apal, const Pals& bpal) const
  {
    
    
    //std::cout << "-----\n";
    
    double ret=0.0;
    
    Pals::PalMap::const_iterator a,b;

    for (a=apal.pals_.begin(); a!=apal.pals_.end(); ++a) {
      //std::cout << "--\n";
      for (b=bpal.pals_.begin(); b!=bpal.pals_.end(); ++b) {
	int dist=std::abs(a->first.dist_-b->first.dist_);
	if (tolerance_<0 ||
	    mismatches(a->first.key_,b->first.key_)<=tolerance_) {
#if 0
	  std::cout << a->first.key_ << " " << a->second << " " << a->first.dist_ << " "
		    << b->first.key_ << " " << b->second << " " << b->first.dist_ << " "
		    << dist << " "
		    << (1/dist)*a->second*b->second << std::endl;
	            << exp(-dist)*a->second*b->second << std::endl;
#endif
	  ret += exp(-dist)*a->second*b->second;
	  //ret += (1/dist)*a->second*b->second;
	}
      }
    }
    return ret;
  }

private:
  int
  mismatches(const std::string& x, const std::string& y) const
  {
    assert(x.size()==y.size());
    uint c=0;
    for (uint i=0; i!=x.size(); ++i) if (x[i]!=y[i]) c++;
    return c;
  }

private:
  int tolerance_;
};

#include "../common/kernel_matrix.h"
#include "../common/fasta.h"
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

typedef std::pair<std::string,Pals> Example;
typedef std::vector<Example> ExampleSet;

uint
load_examples(const std::string& label, Fasta& f, ExampleSet& ex,
	      uint seed_length, uint min_loop, uint max_dist)
{
  while (!f.IsEOF() && f.GetNextSeq()) {
    f.ToLower();
    Pals pal(f.Seq(), seed_length, min_loop, max_dist);
    ex.push_back(Example(label,pal));
  }
  return ex.size();
}

namespace po = boost::program_options;

int
main(int argc, char** argv)
{
  
  namespace po = boost::program_options;
  
  typedef double value_type;
  uint seed_length=3;
  uint min_loop=3;
  uint tolerance=1;
  uint max_dist=300;
  uint n_th=1;
  std::string output_file;
  std::string test_norm_output;

#ifdef HAVE_MPI
  MPI::Init(argc, const_cast<char**&>(argv));
  int my_rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("seed-length,s", po::value<uint>(&seed_length)->default_value(3),
     "seed length of palindrome")
    ("min-loop,l", po::value<uint>(&min_loop)->default_value(3),
     "minimum loop length")
    ("tolerance,t", po::value<uint>(&tolerance)->default_value(1),
     "allowable mismatches")
    ("max-distance,m", po::value<uint>(&max_dist)->default_value(300),
     "maximum distance of palindrome")
#if !defined (HAVE_MPI) && defined (HAVE_BOOST_THREAD)
    ("threads,t", po::value<uint>(&n_th)->default_value(1),
     "set the number of threads")
#endif
    ("norm", po::value<std::string>(&test_norm_output),
     "set the filename for norms of test examples")
    ("normalize,n", "normalize the kernel matrix");
  po::variables_map vm;
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).
    options(desc).allow_unregistered().run();
  std::vector<std::string> extra_args =
    collect_unrecognized(parsed.options, po::include_positional);
  std::vector<po::option>::iterator new_end=
    std::remove_if(parsed.options.begin(), parsed.options.end(),
		   bind(&po::option::unregistered, _1) );
  parsed.options.erase(new_end,  parsed.options.end());
  po::store(parsed, vm);
  po::notify(vm);


  if (vm.count("help") || extra_args.size()<3) {
    std::cout << "Kernel Gram Matrix Calculator for Simple Palindrome Kernel"
	      << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] output [label1 training-data1] ... [--test] [label1] [test-data1] ...\n\n"
	      << desc << std::endl;
    return 1;
  }

  output_file = extra_args[0];
  
  // load examples
  ExampleSet train, test;
  bool train_flag=true;
  for (uint i=1; i<extra_args.size(); i+=2) {
    if (extra_args[i]=="--test") {
      train_flag=false;
      i++;
    }
    std::ifstream in(extra_args[i+1].c_str());
    if (!in.is_open()) {
      perror(extra_args[i+1].c_str());
      return 1;
    }
    Fasta fasta(in);
    if (train_flag) {
      load_examples(extra_args[i], fasta, train,
		    seed_length, min_loop, max_dist);
    } else {
      load_examples(extra_args[i], fasta, test,
		    seed_length, min_loop, max_dist);
    }
  }

  KernelFunc kernel(tolerance);
  KernelMatrix<value_type> matrix;
  if (test.empty()) {
    matrix.calculate(train, kernel, vm.count("normalize"), n_th);
  } else {
    matrix.calculate(test, train, kernel, !test_norm_output.empty(), 
		     vm.count("normalize"), n_th);
  }
  //matrix.calculate(test, train, kernel, vm.count("normalize"));

#ifdef HAVE_MPI
  if (my_rank==0) {
#endif
    std::ofstream out(output_file.c_str());
    matrix.print(out);
    if (!test_norm_output.empty()) {
      std::ofstream tout(test_norm_output.c_str());
      for (uint i=0; i!=test.size(); ++i)
	tout << matrix(i) << std::endl;
    }
#ifdef HAVE_MPI
  }
#endif
 



  return 0;
}

