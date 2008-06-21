// $Id:$

#ifndef __INC_SVM_UTIL_H__
#define __INC_SVM_UTIL_H__

#include <iosfwd>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "../libsvm/svm.h"

class SVMPredict
{
public:
  SVMPredict(const char* output_file, const char* model_file,
	     bool predict_probability);
  ~SVMPredict();

  bool is_open() const { return out_.is_open(); }

  void do_svm_predict(double target, const std::vector<svm_node>& x);

  static
  void make_svm_node(unsigned int n, const std::vector<double>& v,
		     std::vector<svm_node>& x);

private:
  std::ofstream out_;
  std::ostringstream s_out_;
  svm_model* model_;
  bool predict_probability_;
  int correct_;
  int total_;
  double error_;
  double sumv_;
  double sumy_;
  double sumvv_;
  double sumyy_;
  double sumvy_;

  const static unsigned int MAX = 10*1024*1024; // 10MB
};

bool
load_sv_index(std::vector<unsigned int>& sv_index,
	      const std::vector<std::string>& models);

#endif	//  __INC_SVM_UTIL_H__

// Local Variables:
// mode: C++
// End:
