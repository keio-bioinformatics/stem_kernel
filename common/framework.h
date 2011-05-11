// $Id$

#ifndef __INC_FRAMEWORK_H__
#define __INC_FRAMEWORK_H__

#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/timer.hpp>
#ifdef HAVE_BOOST_IOSTREAMS
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#include "../common/glob_wrapper.h"
#include "../common/kernel_matrix.h"
#include "../libsvm/svm_util.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

typedef unsigned int uint;

#ifdef HAVE_BOOST_IOSTREAMS
namespace io = boost::iostreams;
#endif

struct Options
{
  // inputs
  std::vector<std::string> labels;
  std::vector<std::string> files;
  std::vector<std::string> pf_files;
  std::vector<std::string> ts_labels;
  std::vector<std::string> ts_files;
  std::vector<std::string> pf_ts_files;
  std::vector<uint> sv_index;
  // outputs
  std::string output;
  std::string norm_output;
  std::vector<std::string> trained_model_file;
  std::vector<std::string> predict_output;
  // opts
  uint n_th;
  bool normalize;
  bool use_pf_scale_file;
  bool predict_only;
  uint skip;
  bool predict_mode;

  Options()
    : labels(), files(), pf_files(), ts_labels(), ts_files(), pf_ts_files(), sv_index(),
      output(), norm_output(), trained_model_file(), predict_output(),
      n_th(1), normalize(false), use_pf_scale_file(false),
      predict_only(false), skip(0), predict_mode(false)
  {}

  void
  add_options(boost::program_options::options_description& opts);

  void
  parse_extra_args(const std::vector<std::string>& extra_args);
};

class Output
{
  typedef double value_type;
  std::ofstream* out_;	// output stream for the kernel matrix
  std::ofstream* tout_;	// output stream for the norm (k(x,x))
  std::vector<SVMPredict*> pout_;
  std::ostringstream s_out_;
  std::ostringstream s_tout_;
  const static uint MAX = 10*1024*1024; // 10MB

public:
  Output(const char* output, const char* norm_output,
	 const std::vector<std::string>& models,
	 const std::vector<std::string>& pout_files);
    
  ~Output();
    
  void output(uint cnt, const std::string& label,
	      const std::vector<value_type>& vec, value_type self);

private:
  std::ofstream* open(const char* file);
    
  void kernel_output(uint cnt, const std::string& label,
		     const std::vector<value_type>& vec);

  void norm_output(value_type self);

  void prob_output(uint cnt, const std::string& label,
		   const std::vector<double>& vec);
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
  App(const K& kernel, const LDF& ldf, const Options& opts)
    : kernel_(kernel), ldf_(ldf), opts_(opts)
  {
  }

  bool execute() const
  {
    return opts_.predict_mode ? predict() : train();
  }

private:
  bool train() const
  {
    bool res=false;
    ExampleSet ex;
    if (opts_.use_pf_scale_file)
      res=load_examples(ex, opts_.labels, opts_.files, opts_.pf_files);
    else
      res=load_examples(ex, opts_.labels, opts_.files);
    if (!res) return false;
    
    KernelMatrix<value_type> matrix;
    double elapsed =
      matrix.calculate(ex, kernel_, opts_.normalize, opts_.n_th);

#ifdef HAVE_MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
      try {
	std::cout << "elapsed time: " << elapsed << "s" << std::endl;
	std::ofstream out(opts_.output.c_str());
	if (!out) throw opts_.output.c_str();
#ifdef HAVE_BOOST_IOSTREAMS
	io::filtering_stream<io::output> fout;
	if (opts_.output.rfind(".gz")+3==opts_.output.size())
	  fout.push(io::gzip_compressor());
	else if (opts_.output.rfind(".bz2")+4==opts_.output.size())
	  fout.push(io::bzip2_compressor());
	fout.push(out);
	matrix.print(fout);
#else
	matrix.print(out);
#endif
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

  bool predict() const
  {
    bool res=false;
    ExampleSet ex;
    if (opts_.use_pf_scale_file)
      res=load_examples(ex, opts_.labels, opts_.files, opts_.pf_files);
    else
      res=load_examples(ex, opts_.labels, opts_.files);
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
	  new Output(!opts_.predict_only ? opts_.output.c_str() : NULL,
		     !opts_.norm_output.empty() ? opts_.norm_output.c_str() : NULL,
		     opts_.trained_model_file, opts_.predict_output);
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
  
    if (opts_.normalize) {
      double e = KernelMatrix<value_type>::
	diagonal(diag, ex, opts_.sv_index, kernel_, opts_.n_th);
#ifdef HAVE_MPI
      if (rank==0) {
#endif
	std::cout << "elapsed time for diagonals: "
		  << e << "s" << std::endl;
#ifdef HAVE_MPI
      }
#endif
    }

    uint cnt=0;
    bool norm = opts_.normalize || !opts_.norm_output.empty();

    for (uint i=0; i!=opts_.ts_files.size(); ++i) {
      double elapsed = 0.0;
#ifdef HAVE_MPI
      if (rank==0) {
#endif
	std::cout << "predicting " << opts_.ts_files[i] << std::flush;
#ifdef HAVE_MPI
      }
#endif
      Glob glob(opts_.ts_files[i].c_str());
      if (glob.empty()) {
	std::ostringstream os;
	os << opts_.ts_files[i] << ": no matches found";
	throw os.str().c_str();
	//return false;
      }
      Glob pf_glob;
      if (opts_.use_pf_scale_file) {
	pf_glob.expand(opts_.pf_ts_files[i].c_str());
	if (pf_glob.empty()) {
	  std::ostringstream os;
	  os << opts_.pf_ts_files[i] << ": no matches found";
	  throw os.str().c_str();
	  //return false;
	}
	if (pf_glob.size() != glob.size()) {
	  std::ostringstream os;
	  os << opts_.pf_ts_files[i] << ": the number of matched files is inconsistent";
	  throw os.str().c_str();
	  //return false;
	}
      }
      Glob::const_iterator p;
      Glob::const_iterator q=pf_glob.begin();
      for (p=glob.begin(); p!=glob.end(); ++p) {
	typename LDF::Loader* loader = NULL;
	if (opts_.use_pf_scale_file)
	  loader = ldf_.get_loader(p->c_str(), (q++)->c_str());
	else
	  loader = ldf_.get_loader(p->c_str());
	if (loader==NULL) return false;
	
	while (true) {
	  boost::timer tm;
	  Data* data = loader->get();
	  elapsed += tm.elapsed();
	  if (data==NULL) break;
	  value_type self;
	  elapsed += KernelMatrix<value_type>::
	    calculate(vec, std::make_pair(opts_.ts_labels[i], *data),
		      ex, opts_.sv_index, kernel_, opts_.n_th,
		      norm ? &self : NULL);
	  delete data;
#ifdef HAVE_MPI
	  if (rank==0) {
#endif
	    if (opts_.normalize) {
	      boost::timer tm2;
	      for (uint j=0; j!=vec.size(); ++j)
		vec[j] /= sqrt(diag[j]*self);
	      elapsed += tm2.elapsed();
	    }
	    ++cnt;
	    out->output(cnt, opts_.ts_labels[i], vec, self);
#ifdef HAVE_MPI
	  }
#endif
	}
	delete loader;
      }
#ifdef HAVE_MPI
      if (rank==0) {
#endif
	std::cout << " (" << elapsed << "s) done." << std::endl;
#ifdef HAVE_MPI
      }
#endif
    }
    if (out) delete out;
    return true;
  }

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
	double elapsed = 0.0;
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
	while (true) {
	  boost::timer tm;
	  Data* d = loader->get();
	  elapsed += tm.elapsed();
	  if (d==NULL) break;
	  ex.push_back(std::make_pair(labels[i], *d));
	  delete d;
	}
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	  std::cout << " (" << elapsed << "s) done." << std::endl;
#ifdef HAVE_MPI
	}
#endif
	delete loader;
      }
    }
    return true;
  }

  bool load_examples(ExampleSet& ex,
		     const std::vector<std::string>& labels,
		     const std::vector<std::string>& files,
		     const std::vector<std::string>& pf_files) const
  {
    assert(labels.size()==files.size());
    assert(labels.size()==pf_files.size());
    for (uint i=0; i!=files.size(); ++i) {
      Glob glob(files[i].c_str());
      Glob pf_glob(pf_files[i].c_str());
      if (glob.empty()) {
	std::ostringstream os;
	os << files[i] << ": no matches found";
	throw os.str().c_str();
	//return false;
      }
      if (pf_glob.empty()) {
	std::ostringstream os;
	os << pf_files[i] << ": no matches found";
	throw os.str().c_str();
	//return false;
      }
      Glob::const_iterator p, q;
      for (p=glob.begin(), q=pf_glob.begin();
	   p!=glob.end() && q!=pf_glob.end(); ++p, ++q) {
	double elapsed = 0.0;
	typename LDF::Loader* loader=ldf_.get_loader(p->c_str(), q->c_str());
	if (loader==NULL) return false;
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	  std::cout << "loading " << *p
		    << " as label " << labels[i] << std::flush;
#ifdef HAVE_MPI
	}
#endif
	while (true) {
	  boost::timer tm;
	  Data* d = loader->get();
	  elapsed += tm.elapsed();
	  if (d==NULL) break;
	  ex.push_back(std::make_pair(labels[i], *d));
	  delete d;
	}
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
	  std::cout << " (" << elapsed << "s) done." << std::endl;
#ifdef HAVE_MPI
	}
#endif
	delete loader;
      }
    }
    return true;
  }

private:
  const K& kernel_;
  const LDF& ldf_;
  const Options& opts_;
};

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

#endif	//__INC_FRAMEWORK_H__

// Local Variables:
// mode: C++
// End:
