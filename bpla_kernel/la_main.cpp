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
#include "bpla_kernel.h"

namespace po = boost::program_options;

const uint N_AA=AAProfileSequence::N_AA;

float blosum62[N_AA][N_AA] = {
/*  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X  */
  { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0}, // A 
  {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1}, // R 
  {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1}, // N 
  {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1}, // D 
  { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2}, // C 
  {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1}, // Q 
  {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1}, // E 
  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1}, // G 
  {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1}, // H 
  {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1}, // I 
  {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1}, // L 
  {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1}, // K 
  {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1}, // M 
  {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1}, // F 
  {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2}, // P 
  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0}, // S 
  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0}, // T 
  {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2}, // W 
  {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1}, // Y 
  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1}, // V 
  {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1}, // B 
  {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1}, // Z 
  { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1}, // X 
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
  float alpha=1.0;
  float beta;
  float ext;
  std::string score_file;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()("help,h", "show this message");
  opts.add_options(desc);

  po::options_description k_desc("Kernel Options");
  k_desc.add_options()
    ("gap,g", po::value<float>(&gap)->default_value(-10.0), "set gap weight")
    ("ext,e", po::value<float>(&ext)->default_value(-1.0), "set extension weight")
    //("alpha,a", po::value<float>(&alpha)->default_value(4.5), "set alpha")
    ("beta,b", po::value<float>(&beta)->default_value(0.11), "set beta")
    ("score", po::value<std::string>(&score_file), "specify the score table file");
  
  desc.add(k_desc);
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
    std::cout << "Kernel Matrix Calculator for LA Kernels" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0]
	      << " [options] output [label1 training-data1] ... [--test] [label1] [test-data1] ...\n\n"
	      << desc << std::endl;
    return 1;
  }

  opts.parse_extra_args(extra_args);

  boost::multi_array<double,2> score_table(boost::extents[N_AA][N_AA]);
  if (!score_file.empty()) {
    read_score_table(score_table, score_file);
  } else {
    for (uint i=0; i!=N_AA; ++i) {
      for (uint j=0; j!=N_AA; ++j) {
        score_table[i][j] = blosum62[i][j];
      }
    }
  }

  bool res = false;
  try {
    typedef DataLoaderFactory<DataLoader<AAMData> > LDF;
    LDF ldf;
    BPLAKernel<double,AAMData> kernel(score_table, true, gap, ext, alpha, beta);
    App<BPLAKernel<double,AAMData>, LDF> app(kernel, ldf, opts);
    res = app.execute();
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
