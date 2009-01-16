// $Id: main.cpp 209 2007-08-20 07:49:21Z satoken $

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <string>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include "../common/framework.h"
#include "data.h"
#include "../common/rna.h"
#include "def_kernel.h"

namespace po = boost::program_options;

int
main(int argc, char** argv)
{

#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
#endif

  Options opts;
  BPMatrix::Options bp_opts;
  bool no_ribosum;
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
    ("use-bp",
     "use base-paring probability weight")
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

  desc.add(k_desc).add(f_desc);
  po::variables_map vm;
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).
    options(desc).allow_unregistered().run();
  std::vector<std::string> extra_args =
    collect_unrecognized(parsed.options, po::include_positional);
  std::vector<po::option>::iterator new_end =
    std::remove_if(parsed.options.begin(), parsed.options.end(),
		   boost::bind(&po::option::unregistered, _1) );
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
    if (!vm.count("use-bp")) {
      if (!no_ribosum) {
	typedef DataLoaderFactory<DataLoader<MData> > LDF;
	typedef StringKernel<double,MData> SuStringKernel;
	LDF ldf;
	SuStringKernel kernel(gap, alpha);
	App<SuStringKernel,LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else {
	typedef DataLoaderFactory<DataLoader<MData> > LDF;
	typedef StringKernel<double,MData> SiStringKernel;
	LDF ldf;
	SiStringKernel kernel(gap, str_match, str_mismatch);
	App<SiStringKernel,LDF> app(kernel, ldf, opts);
	res = app.execute();
      }
    } else {
      if (!no_ribosum) {
	typedef DataLoaderFactory<DataLoader<MData> > LDF;
	typedef StringKernel<double,MData> SuStringKernel;
	LDF ldf(0.0, bp_opts);
	SuStringKernel kernel(gap, alpha);
	App<SuStringKernel,LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else {
	typedef DataLoaderFactory<DataLoader<MData> > LDF;
	typedef StringKernel<double,MData> SiStringKernel;
	LDF ldf(0.0, bp_opts);
	SiStringKernel kernel(gap, str_match, str_mismatch);
	App<SiStringKernel,LDF> app(kernel, ldf, opts);
	res = app.execute();
      }
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
