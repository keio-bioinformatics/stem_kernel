// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <string>
//#include <cstring>
#include <cassert>
//#include <sys/types.h>
//#include <unistd.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "../common/framework.h"
#include "score_table.h"
#include "stem_kernel.h"
#include "string_kernel.h"
#include "ss_kernel.h"
#include "data.h"
//#include "bpmatrix.h"
#include "../common/rna.h"

using namespace boost::lambda;
namespace po = boost::program_options;

int
main(int argc, char** argv)
{
  typedef DataLoaderFactory<DataLoader<SData> > LDF;
  typedef std::pair<std::string, SData> Example;
  typedef std::vector< Example > ExampleSet;
  typedef SimpleScoreTable<SData,double> SiScoreTable;
  typedef StemKernel<SiScoreTable, SData> SiKernel;
  typedef SubstScoreTable<SData,double> SuScoreTable;
  typedef StemKernel<SuScoreTable, SData> SuKernel;
  typedef StemStrKernel<SuScoreTable, SData> SSKernel;

#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
#endif

  Options opts;
  double th=0.001;
  uint folding_method=FOLD;
  //uint win_sz, pair_sz;
  double gap;
  double stack;
  double covar;
  double alpha;
  double beta;
  double loop_gap;
  bool use_string=false;
  bool use_ribosum=false;
  uint len_band=0;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message");
  opts.add_options(desc);
  desc.add_options()
    ("basepair,p", po::value<double>(&th)->default_value(0.01),
     "set the threshold of basepairing probability")
    ("gap,g", po::value<double>(&gap)->default_value(0.5),
     "set the gap weight")
    ("stack,s", po::value<double>(&stack)->default_value(1.3),
     "set the weight for stacking base pairs")
    ("covariant,v", po::value<double>(&covar)->default_value(0.8),
     "set substitution (covariant) weight for base pairs")
#if 0
    ("window-size,w", po::value<uint>(&win_sz)->default_value(0),
     "set the window size for folding RNAs")
    ("pair-width", po::value<uint>(&pair_sz)->default_value(0),
     "set the pair width for pairing bases")
#endif
    ("length-band", po::value<uint>(&len_band)->default_value(0),
     "set the band of difference of the length between bases")
    ("no-ribosum", "do not use the RIBOSUM substitution matrix")
    ("no-string", "do not convolute the string kernel with base pair probabilities")
    ("alpha,a", po::value<double>(&alpha)->default_value(0.2),
     "set the loop weight of the RIBOSUM for the string kernel")
    ("beta,b", po::value<double>(&beta)->default_value(0.3),
     "set the base pair weight of the RIBOSUM for the stem kernel")
    ("loop-gap,G", po::value<double>(&loop_gap)->default_value(0.6),
     "set the gap weight for loop regions")
    ("--log", "use log")
    ;
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

  use_string = !vm.count("no-string");
  use_ribosum = !vm.count("no-ribosum");
  use_log = !vm.count("log");

  opts.parse_extra_args(extra_args);

  bool res = false;
  try {
    if (!res) {
      LDF ldf(folding_method, th);
      if (use_string) {
	SuScoreTable st(gap, beta, loop_gap);
	SSKernel kernel(st, loop_gap, alpha, len_band);
	App<SSKernel,LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else if (use_ribosum) {
	SuScoreTable st(gap, beta, loop_gap);
	SuKernel kernel(st, len_band);
	App<SuKernel,LDF> app(kernel, ldf, opts);
	res = app.execute();
      } else {
	SiScoreTable st(gap, stack, covar, loop_gap);
	SiKernel kernel(st, len_band);
	App<SiKernel,LDF> app(kernel, ldf, opts);
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
