// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "string_kernel.h"
#include "../common/kernel_matrix.h"
#include "../common/example.h"
#include "../common/fasta.h"

namespace po = boost::program_options;

int
main(int argc, char** argv)
{
  typedef double value_type;
  float gap;
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
    ("gap,g", po::value<float>(&gap)->default_value(1.0), "set gap weight")
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
  parsed.options.erase(new_end,  parsed.options.end());
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
    Fasta fasta(in);
    if (train_flag) {
      load_examples(extra_args[i], fasta, train);
    } else {
      load_examples(extra_args[i], fasta, test);
    }
  }

  StringKernel<value_type> kernel(gap);
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
