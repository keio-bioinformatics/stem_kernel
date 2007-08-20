// $Id:$

#ifndef __INC_FRAMEWORK_H__
#define __INC_FRAMEWORK_H__

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <glob.h>
#include <boost/program_options.hpp>
#include "../common/kernel_matrix.h"
#include "../libsvm/svm_util.h"

class Glob
{
public:
  typedef std::list<std::string>::const_iterator const_iterator;
  
  Glob(const char* pattern, uint flag=0) : path_()
  {
    glob_t buf;
    glob(pattern, flag, NULL, &buf);
    for (uint i=0; i!=buf.gl_pathc; ++i) {
      path_.push_back(std::string(buf.gl_pathv[i]));
    }
    globfree(&buf);
  }
  const_iterator begin() const { return path_.begin(); }
  const_iterator end() const { return path_.end(); }
  bool empty() const { return path_.empty(); }
  
private:
  std::list<std::string> path_;
};

struct Options
{
  // inputs
  std::vector<std::string> labels;
  std::vector<std::string> files;
  std::vector<std::string> ts_labels;
  std::vector<std::string> ts_files;
  std::vector<uint> sv_index;
  // outputs
  std::string output;
  std::string norm_output;
  std::vector<std::string> trained_model_file;
  std::vector<std::string> predict_output;
  // opts
  uint n_th;
  bool normalize;
  bool predict_only;
  uint skip;
  bool predict_mode;

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
    res=load_examples(ex, opts_.labels, opts_.files);
    if (!res) return false;
    
    KernelMatrix<value_type> matrix;
    matrix.calculate(ex, kernel_, opts_.normalize, opts_.n_th);

#ifdef HAVE_MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
      try {
	std::ofstream out(opts_.output.c_str());
	if (!out) throw opts_.output.c_str();
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

  bool predict() const
  {
    bool res=false;
    ExampleSet ex;
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
  
    if (opts_.normalize)
      KernelMatrix<value_type>::
	diagonal(diag, ex, opts_.sv_index, kernel_, opts_.n_th);

    uint cnt=0;
    for (uint i=0; i!=opts_.ts_files.size(); ++i) {
      bool norm = opts_.normalize || !opts_.norm_output.empty();

      Glob glob(opts_.ts_files[i].c_str());
      if (glob.empty()) {
	std::ostringstream os;
	os << opts_.ts_files[i] << ": no matches found";
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
	while ((data=loader->get())!=NULL) {
	  value_type self;
	  KernelMatrix<value_type>::
	    calculate(vec, std::make_pair(opts_.ts_labels[i], *data),
		      ex, opts_.sv_index, kernel_, opts_.n_th,
		      norm ? &self : NULL);
	  delete data;
#ifdef HAVE_MPI
	  if (rank==0) {
#endif
	    if (opts_.normalize) {
	      for (uint j=0; j!=vec.size(); ++j)
		vec[j] /= sqrt(diag[j]*self);
	    }
	    ++cnt;
	    out->output(cnt, opts_.ts_labels[i], vec, self);
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
	while ((d=loader->get())!=NULL) {
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
  const Options& opts_;
};

#endif	//__INC_FRAMEWORK_H__
