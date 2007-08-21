// $Id: main.cpp 209 2007-08-20 07:49:21Z satoken $

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <string>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "../common/framework.h"
#include "string_kernel.h"
#include "data.h"
#include "../common/rna.h"

using namespace boost::lambda;
namespace po = boost::program_options;

int
main(int argc, char** argv)
{
  typedef DataLoaderFactory<DataLoader<MData> > LDF;
  typedef std::pair<std::string, MData> Example;
  typedef std::vector< Example > ExampleSet;
  typedef StringKernel<double,MData> SuStringKernel;

#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
#endif

  Options opts;
  double alpha;
  double loop_gap;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()("help,h", "show this message");
  opts.add_options(desc);
  
  po::options_description k_desc("Kernel Options");
  k_desc.add_options()
    ("no-ribosum", "do not use the RIBOSUM substitution matrix")
    ("alpha,a", po::value<double>(&alpha)->default_value(0.2),
     "set the loop weight of the RIBOSUM for the string kernel")
    ("loop-gap,G", po::value<double>(&loop_gap)->default_value(0.6),
     "set the gap weight for loop regions");

  desc.add(k_desc);
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
    std::cout << "Kernel Matrix Calculator for Stem Kernels" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] output [label1 training-data1] ... [--test] [label1] [test-data1] ...\n\n"
	      << desc << std::endl;
    return 1;
  }

  opts.parse_extra_args(extra_args);

  bool res = false;
  try {
    if (!res) {
      LDF ldf;
      SuStringKernel kernel(loop_gap, alpha);
      App<SuStringKernel,LDF> app(kernel, ldf, opts);
      res = app.execute();
    }
  } catch (const char* str) {
#ifdef HAVE_MPI
    if (/*MPI::COMM_WORLD.Get_rank()==0*/ 1) {
#endif
      std::cout << str << std::endl;
#ifdef HAVE_MPI
    }
#endif    
  }

  return res ? 0 : 1;
}
