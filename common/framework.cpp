// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "framework.h"

namespace po = boost::program_options;

void
Options::  
add_options(boost::program_options::options_description& opts)
{
  opts.add_options()
#if !defined (HAVE_MPI) && defined (HAVE_BOOST_THREAD)
    ("threads,t",
     po::value<uint>(&n_th)->default_value(1),
     "set the number of threads")
#endif
    ("normalize,n",
     po::value<bool>(&normalize)->zero_tokens()->default_value(false),
     "normalize the kernel matrix")
    ("norm,x",
     po::value<std::string>(&norm_output),
     "set the filename for norms of test examples")
    ("no-matrix",
     po::value<bool>(&predict_only)->zero_tokens()->default_value(false),
     "do not output matrix")
    ("model",
     po::value< std::vector<std::string> >(&trained_model_file),
     "the model file trained by svm-train if you already have")
    ("predict",
     po::value< std::vector<std::string> >(&predict_output),
     "output file name of prediction results")
#if 0    
    ("skip",
     po::value<uint>(&skip),
     "skip lines")
#endif
    ;
}

void
Options::
parse_extra_args(const std::vector<std::string>& extra_args)
{
  output = extra_args[0];
  predict_mode =
    extra_args.end()!=std::find(extra_args.begin(),
				extra_args.end(), "--test");
  if (!predict_mode) {
    uint n=(extra_args.size()-1)/2;
    labels.resize(n);
    files.resize(n);
    n=0;
    for (uint i=1; i<extra_args.size(); i+=2) {
      labels[n]=extra_args[i];
      files[n]=extra_args[i+1];
      ++n;
    }
  } else {
    uint x = std::find(extra_args.begin(),
		       extra_args.end(), "--test") - extra_args.begin();
    uint n=(x-1)/2;
    uint m=(extra_args.size()-2)/2-n;
    labels.resize(n);
    files.resize(n);
    ts_labels.resize(m);
    ts_files.resize(m);
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
    if (!trained_model_file.empty()) {
      if (!load_sv_index(sv_index, trained_model_file))
	return ;
    }
  }
}

Output::
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

Output::
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

void
Output::
output(uint cnt, const std::string& label,
       const std::vector<value_type>& vec, value_type self)
{
  kernel_output(cnt, label, vec);
  norm_output(self);
  prob_output(cnt, label, vec);
}

std::ofstream*
Output::
open(const char* file)
{
  if (file!=NULL) {
    std::ofstream* os = new std::ofstream(file);
    if (!os->is_open()) throw file;
    return os;
  }
  return NULL;
}

void
Output::
kernel_output(uint cnt, const std::string& label,
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

void
Output::
prob_output(uint cnt, const std::string& label,
	    const std::vector<double>& vec)
{
  std::vector<svm_node> x(vec.size()+2);
  SVMPredict::make_svm_node(cnt, vec, x);
  for (uint k=0; k!=pout_.size(); ++k) {
    pout_[k]->do_svm_predict(atof(label.c_str()), x);
  }
}

void
Output::
norm_output(value_type self)
{
  if (tout_) {
    s_tout_ << self << std::endl;
    if (s_tout_.str().size()>MAX) {
      *tout_ << s_tout_.str();
      s_tout_.str("");
    }
  }
}
