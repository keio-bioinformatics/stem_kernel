// $Id:$
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <iterator>
#include "svm_util.h"

SVMPredict::
SVMPredict(const char* output_file, const char* model_file,
	   bool predict_probability)
  : out_(output_file),
    s_out_(),
    model_(NULL),
    predict_probability_(predict_probability),
    correct_(0), total_(0), error_(0),
    sumv_(0), sumy_(0), sumvv_(0), sumyy_(0), sumvy_(0)
{
  model_ = svm_load_model(model_file);
  if (predict_probability_) {
    int nr_class=svm_get_nr_class(model_);
    std::vector<int> labels(nr_class);
    svm_get_labels(model_, &labels[0]);
    out_ << "labels ";
    std::copy(labels.begin(), labels.end(),
	      std::ostream_iterator<int>(out_, " "));
    out_ << std::endl;
  }
}

SVMPredict::
~SVMPredict()
{
  out_ << s_out_.str();
  if (model_) svm_destroy_model(model_);
}

void
SVMPredict::
do_svm_predict(double target, const std::vector<svm_node>& x)
{
  int svm_type=svm_get_svm_type(model_);
  unsigned int nr_class=svm_get_nr_class(model_);
  double v;
  std::vector<double> prob_estimates(nr_class);
  if (predict_probability_ && (svm_type==C_SVC || svm_type==NU_SVC)) {
    v = svm_predict_probability(model_, &x[0], &prob_estimates[0]);
    s_out_ << v << " ";
    std::copy(prob_estimates.begin(), prob_estimates.end(),
	      std::ostream_iterator<double>(s_out_, " "));
    s_out_ << std::endl;
  } else {
    v = svm_predict(model_, &x[0]);
    std::vector<double> dec_values(nr_class*(nr_class-1)/2);
    svm_predict_values(model_, &x[0], &dec_values[0]);
    unsigned int pos=0;
    s_out_ << target << " ";
    for(unsigned int i=0; i!=nr_class; ++i)
      for(unsigned int j=i+1; j!=nr_class; ++j)
	s_out_ << dec_values[pos++] << " ";
    s_out_ << std::endl;
  }

  if (s_out_.str().size()>MAX) {
    out_ << s_out_.str();
    s_out_.str("");
  }

  if(v == target)
    ++correct_;
  error_ += (v-target)*(v-target);
  sumv_ += v;
  sumy_ += target;
  sumvv_ += v*v;
  sumyy_ += target*target;
  sumvy_ += v*target;
  ++total_;
}

//static
void
SVMPredict::
make_svm_node(unsigned int n, const std::vector<double>& v, std::vector<svm_node>& x)
{
  //x.resize(v.size()+2);
  x[0].index = 0;
  x[0].value = n;
  for (unsigned int i=0; i!=v.size(); ++i) {
    x[i+1].index = i+1;
    x[i+1].value = v[i];
  }
  x[x.size()-1].index = -1;
}
