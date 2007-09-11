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
  bool no_string;
  bool no_ribosum;
  bool use_log;
  double beta, loop_gap, stack, covar;
  uint len_band=0;
  double alpha, gap, str_match, str_mismatch;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()("help,h", "show this message");
  opts.add_options(desc);

  po::options_description f_desc("Folding Options");
  bp_opts.add_options(f_desc);

  po::options_description k_desc("Kernel Options");
  k_desc.add_options()
    ("no-ribosum",
     po::value<bool>(&no_ribosum)->zero_tokens()->default_value(false),
     "do not use the RIBOSUM substitution matrix")
    ("no-string",
     po::value<bool>(&no_string)->zero_tokens()->default_value(false),
     "do not convolute the string kernel")
    ("log",
     po::value<bool>(&use_log)->zero_tokens()->default_value(false),
     "use the logarithm of the kernel");
  
  po::options_description stem_desc("Options for the stem kernel");
  stem_desc.add_options()
    ("beta,b",
     po::value<double>(&beta)->default_value(0.3),
     "weight of the RIBOSUM for the stem kernel")
    ("loop-gap,g",
     po::value<double>(&loop_gap)->default_value(0.5),
     "gap weight for loop regions")
    ("stack,s",
     po::value<double>(&stack)->default_value(1.3),
     "match weight for stacking base pairs (with --no-ribosum)")
    ("covariant,v",
     po::value<double>(&covar)->default_value(0.8),
     "substitution (covariant) weight for base pairs (with --no-ribosum)")
    ("length-band",
     po::value<uint>(&len_band)->default_value(0),
     "the band of difference of the length between bases");

  po::options_description str_desc("Options for the string kernel");
  str_desc.add_options()
    ("alpha,a",
     po::value<double>(&alpha)->default_value(0.2),
     "weight of the RIBOSUM for the string kernel")
    ("gap,G",
     po::value<double>(&gap)->default_value(0.6),
     "gap weight for the string kernel")
    ("match",
     po::value<double>(&str_match)->default_value(1.0),
     "match weight for the string kernel (with --no-ribosum)")
    ("mismatch",
     po::value<double>(&str_mismatch)->default_value(0.8),
     "substitution (mismatch) weight for the string kernel (with --no-ribosum)");

  desc.add(k_desc).add(stem_desc).add(str_desc).add(f_desc);
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

  opts.parse_extra_args(extra_args);

  bool res = false;
  try {
    typedef DataLoaderFactory<DataLoader<MData> > LDF;
    LDF ldf(bp_opts);
    if (!no_string && !no_ribosum) {
      if (!use_log) {
	SuStemStrKernel<double, MData>
	  kernel(alpha, beta, loop_gap, gap, len_band);
	App<SuStemStrKernel<double,MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else {
	LSuStemStrKernel<double, MData>
	  kernel(alpha, beta, loop_gap, gap, len_band);
	App<LSuStemStrKernel<double, MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      }
    } else if (no_string && !no_ribosum) {
      if (!use_log) {
	SuStemKernel<double, MData>
	  kernel(loop_gap, beta, len_band);
	App<SuStemKernel<double, MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else {
	LSuStemKernel<double, MData>
	  kernel(loop_gap, beta, len_band);
	App<LSuStemKernel<double, MData>, LDF> app(kernel, ldf, opts);
	res = app.execute();
      }
    } else if (!no_string && no_ribosum) {
      SiStemStrKernel<double,MData>
	kernel(loop_gap, stack, covar, gap, str_match, str_mismatch, len_band);
      App<SiStemStrKernel<double, MData>, LDF> app(kernel, ldf, opts);
      res = app.execute();
    } else /*if (no_string && no_ribosum)*/ {
      SiStemKernel<double, MData>
	kernel(loop_gap, stack, covar, len_band);
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
