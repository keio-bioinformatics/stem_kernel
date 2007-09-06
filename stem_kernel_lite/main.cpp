// $Id$

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
#include "data.h"
#include "../common/rna.h"
#include "def_kernel.h"

using namespace boost::lambda;
namespace po = boost::program_options;

int
main(int argc, char** argv)
{

#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
#endif

  Options opts;
  BPMatrix::Options bp_opts;
  double gap;
  double stack;
  double covar;
  double alpha;
  double beta;
  double loop_gap;
  bool use_string=false;
  bool use_ribosum=false;
  bool use_log=false;
  uint len_band=0;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()("help,h", "show this message");
  opts.add_options(desc);

  po::options_description f_desc("Folding Options");
  bp_opts.add_options(f_desc);

  po::options_description k_desc("Kernel Options");
  k_desc.add_options()
    ("gap,g",
     po::value<double>(&gap)->default_value(0.5),
     "set the gap weight")
    ("stack,s",
     po::value<double>(&stack)->default_value(1.3),
     "set the weight for stacking base pairs")
    ("covariant,v",
     po::value<double>(&covar)->default_value(0.8),
     "set substitution (covariant) weight for base pairs")
    ("length-band",
     po::value<uint>(&len_band)->default_value(0),
     "set the band of difference of the length between bases")
    ("no-ribosum", "do not use the RIBOSUM substitution matrix")
    ("no-string",
     "do not convolute the string kernel with base pair probabilities")
    ("alpha,a",
     po::value<double>(&alpha)->default_value(0.2),
     "set the loop weight of the RIBOSUM for the string kernel")
    ("beta,b",
     po::value<double>(&beta)->default_value(0.3),
     "set the base pair weight of the RIBOSUM for the stem kernel")
    ("loop-gap,G",
     po::value<double>(&loop_gap)->default_value(0.6),
     "set the gap weight for loop regions")
    ("log",
     po::value<bool>(&use_log)->zero_tokens()->default_value(false),
     "use the logarithm of the kernel")
    ;

  desc.add(k_desc).add(f_desc);
  po::variables_map vm;
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).
    options(desc).allow_unregistered().run();
  std::vector<std::string> extra_args =
    collect_unrecognized(parsed.options, po::include_positional);
  std::vector<po::option>::iterator new_end =
    std::remove_if(parsed.options.begin(), parsed.options.end(),
		   bind(&po::option::unregistered, _1) );
  parsed.options.erase(new_end, parsed.options.end());
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

  use_string = !vm.count("no-string");
  use_ribosum = !vm.count("no-ribosum");

  opts.parse_extra_args(extra_args);

  bool res = false;
  try {
    typedef DataLoaderFactory<DataLoader<MData> > LDF;
    LDF ldf(bp_opts);
    if (use_string && use_ribosum) {
      if (!use_log) {
	SuStemStrKernel<double,MData>
	  kernel(alpha, beta, gap, loop_gap, len_band);
	App<SuStemStrKernel<double,MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else {
	LSuStemStrKernel<double, MData>
	  kernel(alpha, beta, gap, loop_gap, len_band);
	App<LSuStemStrKernel<double, MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      }
    } else if (!use_string && use_ribosum) {
      if (!use_log) {
	SuStemKernel<double, MData>
	  kernel(gap, beta, len_band);
	App<SuStemKernel<double, MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else {
	LSuStemKernel<double, MData>
	  kernel(gap, beta, len_band);
	App<LSuStemKernel<double, MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      }
    } else if (use_string && !use_ribosum) {
      SiStemStrKernel<double,MData>
	kernel(gap, stack, covar, loop_gap, stack, covar, len_band);
      App<SiStemStrKernel<double, MData>, LDF> app(kernel, ldf, opts);
      res = app.execute();
    } else /*if (!use_string && !use_ribosum)*/ {
      SiStemKernel<double, MData>
	kernel(gap, stack, covar, len_band);
      App<SiStemKernel<double, MData>, LDF> app(kernel, ldf, opts);
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
