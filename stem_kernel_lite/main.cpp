// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cassert>
#include <sys/types.h>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#ifdef HAVE_BOOST_IOSTREAMS
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "../common/kernel_matrix.h"
#include "fa.h"
#include "maf.h"
#include "aln.h"
#include "score_table.h"
#include "stem_kernel.h"
#include "string_kernel.h"
#include "ss_kernel.h"
#include "data.h"
#include "model.h"
#include "bpmatrix.h"
#include "../common/rna.h"
//#include "../libsvm/svm_util.h"

using namespace boost::lambda;
namespace po = boost::program_options;
#ifdef HAVE_BOOST_IOSTREAMS
namespace io = boost::iostreams;
#endif

// options
static double gap;
static double stack;
static double covar;
static double alpha;
static double beta;
static double loop_gap;
static uint n_th;
static bool enable_normalize=false;
static bool enable_test_normalize=false;
static bool use_string_only=false;
static bool use_string=false;
static bool use_ribosum=false;
static std::string test_norm_output;
static std::vector<std::string> trained_model_file;
static std::vector<std::string> predict_output;
static bool predict_only=false;
static uint skip=0;
static uint len_band=0;

template <class Kernel, class ExSet>
static
bool
do_train(const std::string& output_file,
	 const Kernel& kernel, const ExSet& train)
{
  KernelMatrix<double> matrix;
  matrix.calculate(train, kernel, enable_normalize, n_th);

#ifdef HAVE_MPI
  if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
    try {
      std::ofstream out(output_file.c_str());
      if (!out) throw output_file.c_str();
#ifdef HAVE_BOOST_IOSTREAMS
      io::filtering_stream<io::output> fout;
      if (output_file.rfind(".gz")+3==output_file.size())
	fout.push(io::gzip_compressor());
      else if (output_file.rfind(".bz2")+4==output_file.size())
	fout.push(io::bzip2_compressor());
      fout.push(out);
      matrix.print(fout);
#else
      matrix.print(out);
#endif
      if (!test_norm_output.empty()) {
	std::ofstream tout(test_norm_output.c_str());
	if (!tout) throw test_norm_output.c_str();
	for (uint i=0; i!=matrix.self().size(); ++i)
	  tout << matrix.self()[i] << std::endl;
      }
    } catch (const char* f) {
#ifdef HAVE_MPI
      MPI::COMM_WORLD.Abort(1);
#endif
      std::ostringstream os;
      os << f << ": cannot open for writing";
      throw os.str().c_str();
      return false;
    }
#ifdef HAVE_MPI
  }
#endif
  return true;
}

template <class Kernel, class ExSet>
static
bool
do_predict(const std::string& output_file,
	   const std::vector<std::string>& extra_args,
	   const Kernel& kernel, const ExSet& train,
	   const std::vector<uint>& sv_index)
{
  typedef typename Kernel::value_type value_type;
  typedef typename ExSet::value_type Example;
  typedef typename Example::second_type Data;
  
#ifdef HAVE_MPI
  uint rank=MPI::COMM_WORLD.Get_rank();
#endif

  bool do_svm_predict=false;
  if (trained_model_file.size()==predict_output.size())
    do_svm_predict=true;
  else
    predict_only=false;

#if 0
  std::vector<SVMPredict*> pout; 
  std::vector<svm_node> x;
#endif
  std::ofstream out, tout;
  std::vector<value_type> diag;
  std::vector<value_type> vec;
  
#ifdef HAVE_MPI
  if (rank==0) {
#endif
    diag.resize(train.size());
    vec.resize(train.size());

    if (!predict_only) {
      try {
	out.open(output_file.c_str());
	if (!out) throw output_file.c_str();
	if (!test_norm_output.empty()) {
	  tout.open(test_norm_output.c_str());
	  if (!tout) throw test_norm_output.c_str();
	}
      } catch (const char* f) {
#ifdef HAVE_MPI
	MPI::COMM_WORLD.Abort(1);
#endif
	std::ostringstream os;
	os << f << ": cannot open for writing";
	throw os.str().c_str();
	return false;
      }
    }
#if 0
    if (do_svm_predict) {
      x.resize(train.size()+2);
      pout.resize(predict_output.size(), NULL);
      try {
	for (uint k=0; k!=predict_output.size(); ++k) {
	  pout[k] = new SVMPredict(predict_output[k].c_str(),
				   trained_model_file[k].c_str(), true);
	  if (!pout[k]->is_open()) throw predict_output[k].c_str();
	}
      } catch (const char* f) {
	for (uint k=0; k!=predict_output.size(); ++k) 
	  if (pout[k]) delete pout[k];
	std::cout << f << ": cannot open for writing" << std::endl;
#ifdef HAVE_MPI
	MPI::COMM_WORLD.Abort(1);
#endif
	std::ostringstream os;
	os << f << ": cannot open for writing";
	throw os.str().c_str();
	return false;
      }
    }
#endif
#ifdef HAVE_MPI
  }
#endif
  
  if (enable_normalize)
    KernelMatrix<value_type>::diagonal(diag, train, sv_index, kernel, n_th);

  uint cnt=0;
  bool test_flag=false;
  for (uint i=1; i<extra_args.size(); i+=2) {
    if (extra_args[i]=="--test") {
      test_flag=true;
      i++;
    }
    if (!test_flag) continue;

#ifdef HAVE_MPI
    if (rank==0) {
#endif
      std::cout << "predicting " << extra_args[i+1] << std::endl;
#ifdef HAVE_MPI
    }
#endif

    bool norm = enable_normalize || !test_norm_output.empty();
    MakeData<Data> mkdata(extra_args[i+1].c_str());
    std::string label = extra_args[i];
    //double target = atof(label.c_str());
    while (true) {
      Data data;
      if (mkdata(data, cnt<skip)) {
	if (cnt++<skip) continue;

	Example ex(label, data);
	value_type self;
	KernelMatrix<value_type>::calculate(vec, ex, train, sv_index,
					    kernel, n_th, norm ? &self : NULL);
#ifdef HAVE_MPI
	if (rank==0) {
#endif
	  if (enable_normalize) {
	    for (uint j=0; j!=vec.size(); ++j)
	      vec[j] /= sqrt(diag[j]*self);
	  }
	  if (!predict_only) {
	    out << label << " 0:" << cnt << " ";
	    for (uint j=0; j!=vec.size(); ++j) {
	      out << (j+1) << ":" << vec[j] << " ";
	    }
	    out << std::endl;
	    if (!test_norm_output.empty()) tout << self << std::endl;
	  }
#if 0
	  if (do_svm_predict) {
	    SVMPredict::make_svm_node(cnt, vec, x);
	    for (uint k=0; k!=predict_output.size(); ++k) {
	      pout[k]->do_svm_predict(target, x);
	    }
	  }
#endif
#ifdef HAVE_MPI
	}
#endif
      } else {
	break;
      }
    }
  }

  if (do_svm_predict) {
#ifdef HAVE_MPI
    if (rank==0) {
#endif
#if 0
      for (uint k=0; k!=predict_output.size(); ++k) {
	delete pout[k];
      }
#endif
#ifdef HAVE_MPI
    }
#endif
  }

  return true;
}

template < class Seq, class Data >
static 
bool
do_it(const std::vector<std::string>& extra_args)
{
  typedef std::pair<std::string, Seq> Example;
  typedef std::vector< Example > ExampleSet;
  typedef SimpleScoreTable<Data,double> SiScoreTable;
  typedef StemKernel<SiScoreTable, Data> SiKernel;
  typedef SubstScoreTable<Data,double> SuScoreTable;
  typedef StemKernel<SuScoreTable, Data> SuKernel;
  typedef StemStrKernel<SuScoreTable, Data> SSKernel;
  typedef StringKernel<double,Data> SuStringKernel;

  bool predict_mode=
    extra_args.end()!=std::find(extra_args.begin(), extra_args.end(), "--test");
  
  if (!predict_mode) {
    // train mode
    ExampleSet train;
    for (uint i=1; i<extra_args.size(); i+=2) {
      bool res=false;
      res=load_examples(extra_args[i], extra_args[i+1].c_str(), train);
      if (!res) {
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	  std::cout << "fail to load " << extra_args[i+1] << std::endl;
#ifdef HAVE_MPI
	}
#endif
	return false;
      }
#ifdef HAVE_MPI
      if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	std::cout << "load " << extra_args[i+1]
		  << " as label " << extra_args[i] << std::endl;
#ifdef HAVE_MPI
      }
#endif
    }
    // calculate the matrix
    if (use_string_only) {
      SuStringKernel kernel(loop_gap, alpha);
      return do_train(extra_args[0], kernel, train);
    } else if (use_string) {
      SuScoreTable st(gap, beta, loop_gap);
      SSKernel kernel(st, loop_gap, alpha, len_band);
      return do_train(extra_args[0], kernel, train);
    } else if (use_ribosum) {
      SuScoreTable st(gap, beta, loop_gap);
      SuKernel kernel(st, len_band);
      return do_train(extra_args[0], kernel, train);
    } else {
      SiScoreTable st(gap, stack, covar, loop_gap);
      SiKernel kernel(st, len_band);
      return do_train(extra_args[0], kernel, train);
    }
  } else {
    // predict mode
    std::vector<uint> sv_index;
    if (!trained_model_file.empty()) {
      if (!load_sv_index(sv_index, trained_model_file))
	return false;
    }

    // load examples
    ExampleSet train;
    for (uint i=1; i<extra_args.size(); i+=2) {
      if (extra_args[i]=="--test") {
	break;
      }
      bool res=false;

#ifdef HAVE_MPI
      if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	std::cout << "loading " << extra_args[i+1]
		  << " as label " << extra_args[i] << std::flush;
#ifdef HAVE_MPI
      }
#endif

#ifdef HAVE_MPI
      uint rank = MPI::COMM_WORLD.Get_rank();
      uint num_procs = MPI::COMM_WORLD.Get_size();
      if (sv_index.empty())
	res=load_examples(extra_args[i], extra_args[i+1].c_str(), train,
			  num_procs, rank);
      else
	res=load_examples(extra_args[i], extra_args[i+1].c_str(), train,
			  sv_index, num_procs, rank);
#else
      if (sv_index.empty())
	res=load_examples(extra_args[i], extra_args[i+1].c_str(), train);
      else
	res=load_examples(extra_args[i], extra_args[i+1].c_str(), train,
			  sv_index);
#endif

#ifdef HAVE_MPI
      if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	if (res)
	  std::cout << " done." << std::endl;
	else
	  std::cout << " fail." << std::endl;
#ifdef HAVE_MPI
      }
#endif
      if (!res) return false;
    }

    // calculate the matrix
    if (use_string_only) {
      SuStringKernel kernel(loop_gap, alpha);
      return do_predict(extra_args[0], extra_args, kernel, train, sv_index);
    } else if (use_string) {
      SuScoreTable st(gap, beta, loop_gap);
      SSKernel kernel(st, loop_gap, alpha, len_band);
      return do_predict(extra_args[0], extra_args, kernel, train, sv_index);
    } else if (use_ribosum) {
      SuScoreTable st(gap, beta, loop_gap);
      SuKernel kernel(st, len_band);
      return do_predict(extra_args[0], extra_args, kernel, train, sv_index);
    } else {
      SiScoreTable st(gap, stack, covar, loop_gap);
      SiKernel kernel(st, len_band);
      return do_predict(extra_args[0], extra_args, kernel, train, sv_index);
    }
  }

  return false;
}

template < class Data >
static 
bool
do_it(const std::vector<std::string>& extra_args)
{
  return do_it<Data, Data>(extra_args);
}

#ifdef HAVE_MPI
class MPIState
{
public:
  MPIState(int& argc, char**& argv)
  {
    MPI::Init(argc, const_cast<char**&>(argv));
  }

  ~MPIState()
  {
    if (MPI::Is_initialized())
      MPI::Finalize();
  }
};
#endif

int
main(int argc, char** argv)
{
  double th;
  //uint win_sz, pair_sz;
#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
#endif

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("basepair,p", po::value<double>(&th)->default_value(0.01),
     "set the threshold of basepairing probability")
    ("beta,b", po::value<double>(&beta)->default_value(0.1),
     "set the base pair weight of the RIBOSUM for the stem kernel")
    ("gap,g", po::value<double>(&gap)->default_value(0.4),
     "set the gap weight")
    ("alpha,a", po::value<double>(&alpha)->default_value(0.2),
     "set the loop weight of the RIBOSUM for the string kernel")
    ("loop-gap,G", po::value<double>(&loop_gap)->default_value(0.8),
     "set the gap weight for loop regions")
#if 0
    ("stack,s", po::value<double>(&stack)->default_value(1.3),
     "set the weight for stacking base pairs")
    ("covariant,v", po::value<double>(&covar)->default_value(0.8),
     "set substitution (covariant) weight for base pairs")
    ("use-alifold", "use basepairing probablity matrices by pf version of the alifold instead of averaged basepairing probability matrices")
    ("window-size,w", po::value<uint>(&win_sz)->default_value(0),
     "set the window size for folding RNAs")
    ("pair-width", po::value<uint>(&pair_sz)->default_value(0),
     "set the pair width for pairing bases")
#endif
    ("pf-scale", "calculate appropriate pf_scales by MFE")
    ("length-band", po::value<uint>(&len_band)->default_value(0),
     "set the band of difference of the length between bases")
    //("no-ribosum", "do not use the RIBOSUM substitution matrix")
    ("no-string", "do not convolute the la kernel with base pair probabilities")
    //("string-only", "use only the string kernel with base pair probabilities")
    ("la-kernel", "run as local alignment kernel")
#if !defined (HAVE_MPI) && defined (HAVE_BOOST_THREAD)
    ("threads,t", po::value<uint>(&n_th)->default_value(1),
     "set the number of threads")
#endif
    ("normalize,n", "normalize the kernel matrix")
    ("norm,x", po::value<std::string>(&test_norm_output),
     "set the filename for norms of test examples")
#if 0
    ("no-matrix", "do not output matrix")
    ("model", po::value< std::vector<std::string> >(&trained_model_file),
     "the model file trained by svm-train if you already have")
    ("predict", po::value< std::vector<std::string> >(&predict_output),
     "output file name of prediction results")
    ("skip", po::value<uint>(&skip), "skip lines")
#endif
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
  enable_normalize = vm.count("normalize");
  enable_test_normalize = !test_norm_output.empty();
  //use_string_only = vm.count("string-only");
  use_string_only = vm.count("la-kernel");
  use_string = !vm.count("no-string");
  use_ribosum = !vm.count("no-ribosum");
  set_folding_method(vm.count("use-alifold") ? ALIFOLD : FOLD);
  set_use_pf_scale(vm.count("pf-scale"));
  //if (win_sz>0) set_folding_method(LFOLD);
  if (use_string_only) set_folding_method(NO_BPMATRIX);
  set_bp_threshold(th);
  //set_window_size(win_sz, pair_sz);
  predict_only = vm.count("no-matrix");

  bool res = false;
  try {
    if (!res) res = do_it<SData>(extra_args);
    if (!res) res = do_it<MData>(extra_args);
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
