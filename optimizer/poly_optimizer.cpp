// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "optimizer.h"
#include "poly_kernel.h"
#include "gradient.h"

using namespace boost::lambda;
namespace po = boost::program_options;

int
main(int argc, char* argv[])
{
  int degree;
  double gamma;
  double coef0;
  double C;
  double eps, factr, pgtol;
  uint fold;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("degree,d",
     po::value<int>(&degree)->default_value(3),
     "set degree in kernel function")
    ("gamma,g",
     po::value<double>(&gamma)->default_value(1.0),
     "set gamma in kernel function")
    ("coef0,r",
     po::value<double>(&coef0)->default_value(0),
     "set coef0 in kernel function")
    ("cost,c",
     po::value<double>(&C)->default_value(1.0),
     "set the parameter C of C-SVC")
    ("epsilon,e",
     po::value<double>(&eps)->default_value(0.001, "0.001"),
     "set tolerance of termination criterion for the SVM training")
    ("factr",
     po::value<double>(&factr)->default_value(1e10, "1e10"),
     "set the factor for machine epsilon used for the convergence tolerance of iterations")
    ("pgtol",
     po::value<double>(&pgtol)->default_value(1e-5, "1e-5"),
     "set the tolerance of the maximum grandient of the object function")
    ("nfold,v",
     po::value<uint>(&fold)->default_value(10),
     "n-fold cross validation mode");

  po::variables_map vm;
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).
    options(desc).
    allow_unregistered().
    run();
  std::vector<std::string> extra_args =
    collect_unrecognized(parsed.options, po::include_positional);
  std::vector<po::option>::iterator new_end =
    std::remove_if(parsed.options.begin(), parsed.options.end(),
		   bind(&po::option::unregistered, _1) );
  parsed.options.erase(new_end, parsed.options.end());
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help") || extra_args.empty()) {
    std::cout << "Hyperparameter optimizer for polynomial kernels" << std::endl
	      << "Usage:" << std::endl
	      << " " << argv[0] << " [options] training-data"
	      << std::endl << std::endl << desc << std::endl;
    return 1;
  }

  std::vector<int> labels;
  std::vector<PolyKernel::Data> data;

  std::cout << "loading \"" << extra_args[0] << "\"" << std::flush;
  std::ifstream is(extra_args[0].c_str());
  PolyKernel::read_data(is, labels, data);
  std::cout << " done."  << std::endl;

  std::vector<double> param(2);
  std::vector<double> lbd(2);
  std::vector<double> ubd(2);
  std::vector<long> nbd(2);
  param[0] = gamma; lbd[0] = 1e-5; nbd[0] = LBFGSB::LOWER_BOUND;
  param[1] = coef0; nbd[1] = LBFGSB::UNBOUND;
  Optimizer<PolyKernel,GradientComputationAUC> opt;
  opt.set_degree(degree);
  opt.optimize(labels, data, fold, param, C, lbd, ubd, nbd,
	       eps, factr, pgtol, &std::cout);

  std::cout << std::endl
	    << "Optimized Parameters:" << std::endl
	    << "  C=" << C << ", "
	    << "gamma=" << param[0]
	    << " coef0=" << param[1]
	    << std::endl;
    
  return 0;
}
