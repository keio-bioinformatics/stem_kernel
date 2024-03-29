// $Id: main.cpp 242 2007-10-01 07:29:51Z satoken $

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include "../common/framework.h"
#include "data.h"
#include "../common/rna.h"
#include "bpla_kernel.h"

namespace po = boost::program_options;

float default_score_table[4][4]={
  /* A         C          G          U       */
  { 5.846613, -1.860000, -1.460000, -1.390000}, // A
  {-1.860000,  4.786613, -2.480000, -1.050000}, // C
  {-1.460000, -2.480000,  4.656613, -1.740000}, // G
  {-1.390000, -1.050000, -1.740000,  5.276613}, // U
};

template <class T>
bool
read_score_table(boost::multi_array<T,2>& score_table, const std::string& filename)
{
  
  std::ifstream in(filename.c_str());
  if (in.fail()) { perror(filename.c_str()); return false; }
  std::string a, b;
  float v;
  while (in >> a >> b >> v) {
    score_table[char2rna(a[0])][char2rna(b[0])] = v;
  }
  return true;
}

int
main(int argc, char** argv)
{

#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
#endif

  Options opts;
  BPMatrix::Options bp_opts;
  float gap;
  float alpha;
  float beta;
  float ext;
  std::string score_file;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()("help,h", "show this message");
  opts.add_options(desc);

  po::options_description f_desc("Folding Options");
  bp_opts.add_options(f_desc);

  po::options_description k_desc("Kernel Options");
  k_desc.add_options()
    ("noBP", "do not use base-pairing profiles, i.e. run as local alignment kernels.")
    ("SW", "Smith-Waterman kernel instead of local alignment kernels.")
    ("gap,g", po::value<float>(&gap)->default_value(-8.0), "set gap weight")
    ("ext,e", po::value<float>(&ext)->default_value(-0.75), "set extension weight")
    ("alpha,a", po::value<float>(&alpha)->default_value(4.5), "set alpha")
    ("beta,b", po::value<float>(&beta)->default_value(0.11), "set beta")
    ("score", po::value<std::string>(&score_file), "specify the score table file");
  
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
    std::cout << "Kernel Matrix Calculator for BPLA Kernels" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] output [label1 training-data1] ... [--test] [label1] [test-data1] ...\n\n"
	      << desc << std::endl;
    return 1;
  }

  opts.parse_extra_args(extra_args);

  boost::multi_array<double,2> score_table(boost::extents[N_RNA][N_RNA]);
  if (!score_file.empty()) {
    read_score_table(score_table, score_file);
  } else {
    for (uint i=0; i!=N_RNA; ++i) {
      for (uint j=0; j!=N_RNA; ++j) {
        score_table[i][j] = default_score_table[i][j];
      }
    }
  }

  bool res = false;
  try {
    if (vm.count("noBP")) {
      typedef DataLoaderFactory<DataLoader<MData> > LDF;
      LDF ldf;
      BPLAKernel<double,MData> kernel(score_table, vm.count("noBP"), vm.count("SW"), gap, ext, alpha, beta);
      App<BPLAKernel<double,MData>, LDF> app(kernel, ldf, opts);
      res = app.execute();
    } else {
      typedef DataLoaderFactory<DataLoader<MData> > LDF;
      LDF ldf(bp_opts);
      BPLAKernel<double,MData> kernel(score_table, vm.count("noBP"), vm.count("SW"), gap, ext, alpha, beta);
      App<BPLAKernel<double,MData>, LDF> app(kernel, ldf, opts);
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
