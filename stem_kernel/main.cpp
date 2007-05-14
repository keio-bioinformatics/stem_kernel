// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <memory>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "stem_kernel.h"
#include "../common/kernel_matrix.h"
#include "../common/example.h"
#include "../common/fasta.h"

namespace po = boost::program_options;
using namespace boost::lambda;

int
main(int argc, char* argv[])
{
  typedef double value_type;
  float gap;
  float stack;
  float subst;
  uint loop;
  uint n_th=1;
  uint band;
  float bp_bound=1.0;
  float ali_bound;
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
    ("enable-wobble-pair,w", "allow wobble base pairs")
    ("gap,g", po::value<float>(&gap)->default_value(0.8),
     "set the gap weight")
    ("stack,s", po::value<float>(&stack)->default_value(1.0),
     "set the stacking weight")
    ("loop,l", po::value<uint>(&loop)->default_value(3),
     "set minimum loop length")
    ("substitution,v", po::value<float>(&subst)->default_value(0.5),
     "set substitution weight for base pairs")
#ifdef HAVE_LIBRNA
    ("basepair-probability,p", po::value<float>(&bp_bound)->default_value(0.0),
     "weight all basepairs by its basepair probability calculated by pf_fold()")
#endif
    ("alignment-constraint,a", po::value<float>(&ali_bound)->default_value(0.0),
     "use alignment constraints with posteior probability over this value")
#if !defined (HAVE_MPI) && defined (HAVE_BOOST_THREAD)
    ("threads,t", po::value<uint>(&n_th)->default_value(1),
     "set the number of threads")
#endif
    ("band-width,b", po::value<uint>(&band)->default_value(0),
     "set band width")
    ("normalize,n", "normalize the kernel matrix")
    ("norm,x", po::value<std::string>(&test_norm_output),
     "set the filename for norms of test examples");
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
    std::cout << "Kernel Gram Matrix Calculator for Stem Kernel" << std::endl
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
      load_examples(extra_args[i], fasta, train);
    } else {
      load_examples(extra_args[i], fasta, test);
    }
  }

  bool use_GU = vm.count("enable-wobble-pair");
  KernelMatrix<value_type> matrix;
#ifdef HAVE_LIBRNA
  if (bp_bound<1.0) {
    typedef StemKernel<value_type,BPMatrix> Kernel;
    Kernel kernel(use_GU, loop, gap, stack, subst, band, ali_bound, bp_bound);
    if (test.empty()) {
      matrix.calculate(train, kernel, vm.count("normalize"), n_th);
    } else {
      matrix.calculate(test, train, kernel, !test_norm_output.empty(), 
		       vm.count("normalize"), n_th);
    }
  } else if (use_GU) {
#else
  if (use_GU) {
#endif
    typedef StemKernel<value_type,WobbleBasePair> Kernel;    
    Kernel kernel(use_GU, loop, gap, stack, subst, band, ali_bound);
    if (test.empty()) {
      matrix.calculate(train, kernel, vm.count("normalize"), n_th);
    } else {
      matrix.calculate(test, train, kernel, !test_norm_output.empty(), 
		       vm.count("normalize"), n_th);
    }
  } else {
    typedef StemKernel<value_type,NormalBasePair> Kernel;    
    Kernel kernel(use_GU, loop, gap, stack, subst, band, ali_bound);
    if (test.empty()) {
      matrix.calculate(train, kernel, vm.count("normalize"), n_th);
    } else {
      matrix.calculate(test, train, kernel, !test_norm_output.empty(), 
		       vm.count("normalize"), n_th);
    }
  }

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

#ifdef HAVE_MPI
  MPI::Finalize();
#endif

  return 0;
}
