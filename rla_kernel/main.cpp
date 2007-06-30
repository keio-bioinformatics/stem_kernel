// $Id: main.cpp 100 2006-11-29 08:15:50Z satoken $

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/bind.hpp>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "string_kernel.h"
#include "../common/kernel_matrix.h"
#include "../common/input.h"
#include "../common/pf_wrapper.h"
#include "../common/fasta.h"

namespace Vienna{
  extern "C" {
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/PS_dot.h>
    extern int pfl_fold(char *sequence, int winSize, int pairdist,
		      float cutoff, struct plist **pl);
    extern void init_pf_foldLP(int length);
    extern void free_pf_arraysLP(void);
  };
};

namespace po = boost::program_options;
using namespace boost::lambda;

typedef double value_type;

typedef std::pair<std::string, StringKernel<value_type>::Seq > Example;
typedef std::vector< Example > ExampleSet;

uint char2rna(char c) 
{
  switch(c){
  case 'a':
    return 0;
    break;
  case 'c':
    return 1;
    break;
  case 'g':
    return 2;
    break;
  case 'u':
    return 3;
    break;
  case 't':
    return 3;
    break;
  default:
    return 4;
    break;
  }
}

int
load_examples(const std::string& label, std::ifstream& in, ExampleSet& ex, 
	      int win_sz_, int pair_sz_, double th_)
{
  typedef double value_type;
  typedef boost::multi_array<value_type,2> BPMAT;
  uint ret=0;
  InputSeq is(in);
  value_type tmp, tmp2, tmp3;
  while (is.get_next_seq()) {
    std::string s(is.seq());
    boost::algorithm::to_lower(s);
    
    //caluculate probability
    Vienna::plist *pl;
    Vienna::init_pf_foldLP(s.size()+1);
    //fold
    Vienna::pfl_fold(const_cast<char*>(s.c_str()),
		     win_sz_, pair_sz_, th_, &pl);
    BPMAT bp(boost::extents[s.size()+1][s.size()+1]);
    for (uint i=0; i!=s.size()+1; ++i) {
      for (uint j=0; j!=s.size()+1; ++j) {
	bp[i][j] = 0.0;
      }
    }
    for (uint k=0; pl[k].i!=0; ++k){
      bp[pl[k].i][pl[k].j] = pl[k].p;
    }
    free(pl);
    Vienna::free_pf_arraysLP();
    
    std::vector<uint> intseq;
    for(uint j=0; j < s.size(); j++){
      intseq.push_back(char2rna(s[j]));
    }
    std::vector<value_type> right;
    std::vector<value_type> left;
    std::vector<value_type> unpair;
   
 
    for(uint j=0; j!=s.size()+1; ++j){
      tmp =0;
      for(uint i=j;;--i){
	tmp+=bp[i+1][j+1];
	if(i==0) break;
      }
      right.push_back(tmp);
      tmp2 = 0;
      for(uint k=j; k!=s.size()+1; ++k){
	tmp2+=bp[j+1][k+1];
      }
      left.push_back(tmp2);
      tmp3 = 0;
      tmp3 = 1 - (tmp + tmp2);
      if(tmp3 < 0) tmp3 = 0;
      
      unpair.push_back(tmp3);
    }
    /*
    for(uint i=0; i != intseq.size();i++){
	std::cout << s[i]<< " "<< intseq[i] <<"\t"<< left[i] <<"\t\t"<< right[i] << "\t\t"<< unpair[i] << std::endl;  
      }
    
      std::cout << std::endl;
    */
    //StringKernel<value_type>::Seq example(s, right, left, unpair);
    StringKernel<value_type>::Seq seq(intseq, right, left, unpair);
    ex.push_back(std::make_pair(label,seq));
    ret++;
  }
  return ret;
}

int
load_examples(const std::string& label,
              std::ifstream& in, ExampleSet& ex)
{
  typedef double value_type;
  uint ret=0;
  InputSeq is(in);
  value_type tmp, tmp2, tmp3;

  while (is.get_next_seq()) {
    std::vector<value_type> right;
    std::vector<value_type> left;
    std::vector<value_type> unpair;
    std::string s(is.seq());
    boost::algorithm::to_lower(s);
    
    //convert from string to uint
    std::vector<uint> intseq;
    for(uint j=0; j < s.size(); j++){
      intseq.push_back(char2rna(s[j]));
    }
        
    //caluculate probability
    PFWrapper pf(s,true);
    float f=pf.fold();
    for(uint j=0; j!=pf.size(); ++j){
      tmp =0;
      for(uint i=j;;--i){
	tmp+=pf(i+1,j+1);
	if(i==0) break;
      }
      right.push_back(tmp);
      tmp2 = 0;
      for(uint k=j; k!=pf.size(); ++k){
	tmp2+=pf(j+1,k+1);
      }
      left.push_back(tmp2);
      tmp3 = 0;
      tmp3 = 1 - (tmp + tmp2);
      if(tmp3 < 0) tmp3 = 0;
      unpair.push_back(tmp3);
    }
    StringKernel<value_type>::Seq seq(intseq, right, left, unpair);
    ex.push_back(std::make_pair(label,seq));
    ret++;
  }
  return ret;
}

int
main(int argc, char** argv)
{
  typedef double value_type;
  float gap;
  float alpha;
  float beta;
  float ext;
  int win_sz;
  int pair_sz;
  double th=10E-100; 
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
    ("gap,g", po::value<float>(&gap)->default_value(-8.0), "set gap weight")
    ("ext,e", po::value<float>(&ext)->default_value(-0.75), "set extension weight")
    ("alpha,a", po::value<float>(&alpha)->default_value(4.5), "set alpha")
    ("beta,b", po::value<float>(&beta)->default_value(0.11), "set beta")
    ("win_sz,w", po::value<int>(&win_sz)->default_value(15), "set win_sz")
    ("pair_sz,p", po::value<int>(&pair_sz)->default_value(15), "set pair_sz")
#if !defined (HAVE_MPI) && defined (HAVE_BOOST_THREAD)
    ("threads,t", po::value<uint>(&n_th)->default_value(1),
     "set the number of threads")
#endif
    ("norm,x", po::value<std::string>(&test_norm_output),
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
  parsed.options.erase(new_end, parsed.options.end());
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help") || extra_args.size()<3) {
    std::cout << "Kernel Gram Matrix Calculator for String Kernel" << std::endl
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
    if (train_flag) {
      load_examples(extra_args[i], in, train);
      //load_examples(extra_args[i], in, train, win_sz, pair_sz, th);
    } else {
      load_examples(extra_args[i], in, test);
      //load_examples(extra_args[i], in, test, win_sz, pair_sz, th);
    }
  }

  StringKernel<value_type> kernel(gap, ext, alpha, beta);
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
