// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cassert>
#include <sys/types.h>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "../common/kernel_matrix.h"
#include "../common/glob.h"
#include "score_table.h"
#include "stem_kernel.h"
#include "string_kernel.h"
#include "ss_kernel.h"
#include "data.h"
#include "model.h"
#include "bpmatrix.h"
#include "../common/rna.h"
#include "../libsvm/svm_util.h"

using namespace boost::lambda;
namespace po = boost::program_options;

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

class Output
{
  typedef double value_type;
  std::ofstream* out_;	// output stream for the kernel matrix
  std::ofstream* tout_;	// output stream for the norm (k(x,x))
  std::vector<SVMPredict*> pout_;
  std::ostringstream s_out_;
  std::ostringstream s_tout_;
  const static int MAX = 10*1024*1024; // 10MB

public:
  Output(const char* output, const char* norm_output,
	 const std::vector<std::string>& models,
	 const std::vector<std::string>& pout_files)
    : out_(open(output)), tout_(open(norm_output)), pout_(),
      s_out_(), s_tout_()
  {
    pout_.resize(pout_files.size(), NULL);
    for (uint k=0; k!=pout_files.size(); ++k) {
      pout_[k] = new SVMPredict(pout_files[k].c_str(),
				models[k].c_str(), true);
      if (!pout_[k]->is_open()) throw pout_files[k].c_str();
    }
  }
    
  ~Output()
  {
    if (out_) {
      *out_ << s_out_.str();
      delete out_;
    }
    if (tout_) {
      *tout_ << s_tout_.str();
      delete tout_;
    }
    for (uint i=0; i!=pout_.size(); ++i)
      if (pout_[i]) delete pout_[i];
  }
    
  void output(uint cnt, const std::string& label,
	      const std::vector<value_type>& vec, value_type self)
  {
    kernel_output(cnt, label, vec);
    norm_output(self);
    prob_output(cnt, label, vec);
  }

private:
  std::ofstream* open(const char* file)
  {
    if (file!=NULL) {
      std::ofstream* os = new std::ofstream(file);
      if (!os->is_open()) throw file;
      return os;
    }
    return NULL;
  }
    
  void kernel_output(uint cnt, const std::string& label,
		     const std::vector<value_type>& vec)
  {
    if (out_) {
      s_out_ << label << " 0:" << cnt << " ";
      for (uint j=0; j!=vec.size(); ++j) {
	s_out_ << (j+1) << ":" << vec[j] << " ";
      }
      s_out_ << std::endl;
      if (s_out_.str().size()>MAX) {
	*out_ << s_out_.str();
	s_out_.str("");
      }
    }
  }

  void norm_output(value_type self)
  {
    if (tout_) {
      s_tout_ << self << std::endl;
      if (s_tout_.str().size()>MAX) {
	*tout_ << s_tout_.str();
	s_tout_.str("");
      }
    }
  }

  void prob_output(uint cnt, const std::string& label,
		   const std::vector<double>& vec)
  {
    std::vector<svm_node> x(vec.size()+2);
    SVMPredict::make_svm_node(cnt, vec, x);
    for (uint k=0; k!=pout_.size(); ++k) {
      pout_[k]->do_svm_predict(atof(label.c_str()), x);
    }
  }
};

template < class K, class LDF >
class App
{
public:
  typedef typename LDF::Data Data;
  typedef std::pair<std::string, Data> Example;
  typedef std::vector<Example> ExampleSet;
  typedef double value_type;

public:
  App(const K& kernel, const LDF& ldf)
    : kernel_(kernel), ldf_(ldf)
  {
  }

  bool train(const std::vector<std::string>& labels,
	     const std::vector<std::string>& files,
	     const std::string& output_file, bool normalize,
	     uint n_th) const
  {
    bool res=false;
    ExampleSet ex;
    res=load_examples(ex, labels, files);
    if (!res) return false;
    
    KernelMatrix<value_type> matrix;
    matrix.calculate(ex, kernel_, normalize, n_th);

#ifdef HAVE_MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
      try {
	std::ofstream out(output_file.c_str());
	if (!out) throw output_file.c_str();
	matrix.print(out);
      } catch (const char* f) {
#ifdef HAVE_MPI
	MPI::COMM_WORLD.Abort(1);
#endif
	std::ostringstream os;
	os << f << ": cannot open for writing";
	throw os.str().c_str();
      }
#ifdef HAVE_MPI
    }
#endif
    return true;
  }

  bool predict(const std::vector<std::string>& labels,
	       const std::vector<std::string>& files,
	       const std::vector<std::string>& ts_labels,
	       const std::vector<std::string>& ts_files,
	       const std::vector<uint>& sv_index,
	       const std::string& output,
	       const std::string& norm_output,
	       const std::vector<std::string>& trained_model_file,
	       const std::vector<std::string>& predict_output,
	       bool normalize, bool predict_only, uint n_th) const
  {
    bool res=false;
    ExampleSet ex;
    res=load_examples(ex, labels, files);
    if (!res) return false;
    
#ifdef HAVE_MPI
    uint rank=MPI::COMM_WORLD.Get_rank();
#endif
    
    Output* out=NULL;
    std::vector<value_type> diag;
    std::vector<value_type> vec;
  
#ifdef HAVE_MPI
    if (rank==0) {
#endif
      diag.resize(ex.size());
      vec.resize(ex.size());

      try {
	out =
	  new Output(!predict_only ? output.c_str() : NULL,
		     !norm_output.empty() ? norm_output.c_str() : NULL,
		     trained_model_file, predict_output);
      } catch (const char* f) {
#ifdef HAVE_MPI
	MPI::COMM_WORLD.Abort(1);
#endif
	std::ostringstream os;
	os << f << ": cannot open for writing";
	throw os.str().c_str();
	if (out) delete out;
	return false;
      }
#ifdef HAVE_MPI
    }
#endif
  
    if (normalize)
      KernelMatrix<value_type>::diagonal(diag, ex, sv_index, kernel_, n_th);

    uint cnt=0;
    for (uint i=0; i!=ts_files.size(); ++i) {
      bool norm = normalize || !norm_output.empty();

      Glob glob(ts_files[i].c_str());
      if (glob.empty()) {
	std::ostringstream os;
	os << ts_files[i] << ": no matches found";
	throw os.str().c_str();
	//return false;
      }
      Glob::const_iterator p;
      for (p=glob.begin(); p!=glob.end(); ++p) {
#ifdef HAVE_MPI
	if (rank==0) {
#endif
	  std::cout << "predicting " << *p << std::endl;
#ifdef HAVE_MPI
	}
#endif
	typename LDF::Loader* loader=ldf_.get_loader(p->c_str());
	if (loader==NULL) return false;
	
	Data* data;
	while (data=loader->get()) {
	  value_type self;
	  KernelMatrix<value_type>::
	    calculate(vec, std::make_pair(ts_labels[i], *data),
		      ex, sv_index, kernel_, n_th, norm ? &self : NULL);
	  delete data;
#ifdef HAVE_MPI
	  if (rank==0) {
#endif
	    if (normalize) {
	      for (uint j=0; j!=vec.size(); ++j)
		vec[j] /= sqrt(diag[j]*self);
	    }
	    ++cnt;
	    out->output(cnt, ts_labels[i], vec, self);
#ifdef HAVE_MPI
	  }
#endif
	}
	delete loader;
      }
    }
    if (out) delete out;
    return true;
  }

private:
  bool load_examples(ExampleSet& ex,
		     const std::vector<std::string>& labels,
		     const std::vector<std::string>& files) const
  {
    assert(labels.size()==files.size());
    for (uint i=0; i!=files.size(); ++i) {
      Glob glob(files[i].c_str());
      if (glob.empty()) {
	std::ostringstream os;
	os << files[i] << ": no matches found";
	throw os.str().c_str();
	//return false;
      }
      Glob::const_iterator p;
      for (p=glob.begin(); p!=glob.end(); ++p) {
	typename LDF::Loader* loader=ldf_.get_loader(p->c_str());
	if (loader==NULL) return false;
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	  std::cout << "loading " << *p
		    << " as label " << labels[i] << std::flush;
#ifdef HAVE_MPI
	}
#endif
	Data* d;
	while (d=loader->get()) {
	  ex.push_back(std::make_pair(labels[i], *d));
	  delete d;
	}
	delete loader;
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	  std::cout << " done." << std::endl;
#ifdef HAVE_MPI
	}
#endif
      }
    }
    return true;
  }

private:
  const K& kernel_;
  const LDF& ldf_;
};

template < class LDF >
static
bool
do_it(const std::vector<std::string>& extra_args,
      const LDF& ldf)
{
  typedef typename LDF::Data Data;
  typedef std::pair<std::string, Data> Example;
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
    uint n=(extra_args.size()-1)/2;
    std::vector<std::string> labels(n);
    std::vector<std::string> files(n);
    n=0;
    for (uint i=1; i<extra_args.size(); i+=2) {
      labels[n]=extra_args[i];
      files[n]=extra_args[i+1];
      ++n;
    }

    if (use_string_only) {
      SuStringKernel kernel(loop_gap, alpha);
      App<SuStringKernel,LDF> app(kernel, ldf);
      return app.train(labels, files, extra_args[0], enable_normalize, n_th);
    } else if (use_string) {
      SuScoreTable st(gap, beta, loop_gap);
      SSKernel kernel(st, loop_gap, alpha, len_band);
      App<SSKernel,LDF> app(kernel, ldf);
      return app.train(labels, files, extra_args[0], enable_normalize, n_th);
    } else if (use_ribosum) {
      SuScoreTable st(gap, beta, loop_gap);
      SuKernel kernel(st, len_band);
      App<SuKernel,LDF> app(kernel, ldf);
      return app.train(labels, files, extra_args[0], enable_normalize, n_th);
    } else {
      SiScoreTable st(gap, stack, covar, loop_gap);
      SiKernel kernel(st, len_band);
      App<SiKernel,LDF> app(kernel, ldf);
      return app.train(labels, files, extra_args[0], enable_normalize, n_th);
    }

  } else {
    uint x = std::find(extra_args.begin(),
		       extra_args.end(), "--test") - extra_args.begin();
    uint n=(x-1)/2;
    uint m=(extra_args.size()-2)/2-n;
    std::vector<std::string> labels(n);
    std::vector<std::string> files(n);
    std::vector<std::string> ts_labels(m);
    std::vector<std::string> ts_files(m);
    n=0; m=0;
    for (uint i=1; i<x; i+=2) {
      labels[n]=extra_args[i];
      files[n]=extra_args[i+1];
      ++n;
    }
    for (uint i=x+1; i<extra_args.size(); i+=2) {
      ts_labels[m]=extra_args[i];
      ts_files[m]=extra_args[i+1];
      ++m;
    }

    // predict mode
    std::vector<uint> sv_index;
    if (!trained_model_file.empty()) {
      if (!load_sv_index(sv_index, trained_model_file))
	return false;
    }

    // calculate the matrix
    if (use_string_only) {
      SuStringKernel kernel(loop_gap, alpha);
      App<SuStringKernel,LDF> app(kernel, ldf);
      return app.predict(labels, files, ts_labels, ts_files,
			 sv_index, extra_args[0], test_norm_output,
			 trained_model_file, predict_output,
			 enable_normalize, predict_only, n_th);
    } else if (use_string) {
      SuScoreTable st(gap, beta, loop_gap);
      SSKernel kernel(st, loop_gap, alpha, len_band);
      App<SSKernel,LDF> app(kernel, ldf);
      return app.predict(labels, files, ts_labels, ts_files,
			 sv_index, extra_args[0], test_norm_output,
			 trained_model_file, predict_output,
			 enable_normalize, predict_only, n_th);
    } else if (use_ribosum) {
      SuScoreTable st(gap, beta, loop_gap);
      SuKernel kernel(st, len_band);
      App<SuKernel,LDF> app(kernel, ldf);
      return app.predict(labels, files, ts_labels, ts_files,
			 sv_index, extra_args[0], test_norm_output,
			 trained_model_file, predict_output,
			 enable_normalize, predict_only, n_th);
    } else {
      SiScoreTable st(gap, stack, covar, loop_gap);
      SiKernel kernel(st, len_band);
      App<SiKernel,LDF> app(kernel, ldf);
      return app.predict(labels, files, ts_labels, ts_files,
			 sv_index, extra_args[0], test_norm_output,
			 trained_model_file, predict_output,
			 enable_normalize, predict_only, n_th);
    }
  }
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
  uint win_sz, pair_sz;
#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
#endif

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("basepair,p", po::value<double>(&th)->default_value(0.01),
     "set the threshold of basepairing probability")
    ("gap,g", po::value<double>(&gap)->default_value(0.5),
     "set the gap weight")
    ("stack,s", po::value<double>(&stack)->default_value(1.3),
     "set the weight for stacking base pairs")
    ("covariant,v", po::value<double>(&covar)->default_value(0.8),
     "set substitution (covariant) weight for base pairs")
    ("use-alifold", "use basepairing probablity matrices by pf version of the alifold instead of averaged basepairing probability matrices")
    ("window-size,w", po::value<uint>(&win_sz)->default_value(0),
     "set the window size for folding RNAs")
    ("pair-width", po::value<uint>(&pair_sz)->default_value(0),
     "set the pair width for pairing bases")
    ("length-band", po::value<uint>(&len_band)->default_value(0),
     "set the band of difference of the length between bases")
    ("no-ribosum", "do not use the RIBOSUM substitution matrix")
    ("no-string", "do not convolute the string kernel with base pair probabilities")
    //("string-only", "use only the string kernel with base pair probabilities")
    ("la-kernel", "run as local alignment kernel")
    ("alpha,a", po::value<double>(&alpha)->default_value(0.2),
     "set the loop weight of the RIBOSUM for the string kernel")
    ("beta,b", po::value<double>(&beta)->default_value(0.3),
     "set the base pair weight of the RIBOSUM for the stem kernel")
    ("loop-gap,G", po::value<double>(&loop_gap)->default_value(0.6),
     "set the gap weight for loop regions")
#if !defined (HAVE_MPI) && defined (HAVE_BOOST_THREAD)
    ("threads,t", po::value<uint>(&n_th)->default_value(1),
     "set the number of threads")
#endif
    ("normalize,n", "normalize the kernel matrix")
    ("norm,x", po::value<std::string>(&test_norm_output),
     "set the filename for norms of test examples")
    ("no-matrix", "do not output matrix")
    ("model", po::value< std::vector<std::string> >(&trained_model_file),
     "the model file trained by svm-train if you already have")
    ("predict", po::value< std::vector<std::string> >(&predict_output),
     "output file name of prediction results")
    ("skip", po::value<uint>(&skip), "skip lines");
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
  predict_only = vm.count("no-matrix");
  uint folding_method=vm.count("use-alifold") ? ALIFOLD : FOLD;
  if (use_string_only) folding_method=NO_BPMATRIX;

  bool res = false;
  try {
    if (!res) {
      DataLoaderFactory<DataLoader<SData> > ldf(folding_method, th);
      res=do_it(extra_args, ldf);
    }
    if (!res) {
      DataLoaderFactory<DataLoader<MData> > ldf(folding_method, th);
      res=do_it(extra_args, ldf);
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
