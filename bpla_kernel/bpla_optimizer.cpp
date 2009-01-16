// $Id$

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <iostream>
#include <string>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/timer.hpp>
#include "../common/glob_wrapper.h"
#include "data.h"
#include "../common/rna.h"
#include "../optimizer/optimizer.h"
#include "bpla_kernel.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace boost::lambda;
namespace po = boost::program_options;

static std::ostream* os=&std::cout;

template < class DataT >
class CalcMatrix
{
public:
  typedef DataT Data;
  
public:
  static void calculate_matrix(const std::vector<Data>& data,
			       const std::vector<double>& param,
			       Kmat& kmat, Gmat& gmat)
  {
#ifdef HAVE_MPI
    uint n_data = data.size();
    uint n_cell = (n_data+1)*n_data/2;
    uint n_proc = MPI::COMM_WORLD.Get_size();
    uint rank = MPI::COMM_WORLD.Get_rank();

    // calculation
    std::vector<double> k_l(n_cell/n_proc+1);
    boost::multi_array<double,2> g_l(boost::extents[param.size()][n_cell/n_proc+1]);
    for (uint i=0, n=0, c=0; i!=data.size(); ++i) {
      for (uint j=i; j!=data.size(); ++j) {
	if (c++ % n_proc == rank) {
	  std::vector<double> g(param.size());
	  k_l[n] = BPLAKernel<double,MData>::
	    compute_gradients(data[i], data[j], param, g);
	  for (uint l=0; l!=g.size(); ++l) g_l[l][n] = g[l];
	  ++n;
	}
      }
    }

    // synchronization
    std::vector<double> k_g(n_cell/n_proc+1);
    boost::multi_array<double,2> g_g(boost::extents[param.size()][n_cell/n_proc+1]);
    for (uint p=0; p!=n_proc; ++p) {
      if (p==rank) {
	std::copy(k_l.begin(), k_l.end(), k_g.begin());
	for (uint l=0; l!=g_g.size(); ++l) 
	  std::copy(g_l[l].begin(), g_l[l].end(), g_g[l].begin());
      }
      MPI::COMM_WORLD.Bcast(&k_g[0], k_g.size(), MPI::DOUBLE, p);
      for (uint l=0; l!=g_l.size(); ++l)
	MPI::COMM_WORLD.Bcast(&g_g[l][0], g_g[l].size(), MPI::DOUBLE, p);

      for (uint i=0, n=0, c=0; i!=data.size(); ++i) {
	for (uint j=i; j!=data.size(); ++j) {
	  if (c++ % n_proc == p) {
	    kmat[i][j] = kmat[j][i] = k_g[n];
	    for (uint l=0; l!=g_g.size(); ++l)
	      gmat[l][i][j] = gmat[l][j][i] = g_g[l][n];
	    ++n;
	  }
	}
      }
    }
#else
    double k;
    std::vector<double> g(param.size());
    for (uint i=0; i!=data.size(); ++i) {
      for (uint j=i; j!=data.size(); ++j) {
	k = BPLAKernel<double,MData>::
	  compute_gradients(data[i], data[j], param, g);
	kmat[i][j] = kmat[j][i] = k;
	for (uint l=0; l!=g.size(); ++l)
	  gmat[l][i][j] = gmat[l][j][i] = g[l];
      }
    }
#endif
  }
};

template < class DataT >
class CalcMatrixN
{
public:
  typedef DataT Data;
  
public:
  static void calculate_matrix(const std::vector<Data>& data,
			       const std::vector<double>& param,
			       Kmat& kmat, Gmat& gmat)
  {
#ifdef HAVE_MPI
    uint n_data = data.size();
    uint n_cell = (n_data+1)*n_data/2;
    uint n_proc = MPI::COMM_WORLD.Get_size();
    uint rank = MPI::COMM_WORLD.Get_rank();

    std::vector<double> diag_k(data.size());
    boost::multi_array<double,2> diag_g(boost::extents[param.size()][data.size()]);

    // calculation (diag)
    std::vector<double> diag_k_l(n_data/n_proc+1);
    boost::multi_array<double,2> diag_g_l(boost::extents[param.size()][n_data/n_proc+1]);
    for (uint i=0, n=0; i!=data.size(); ++i) {
      if (i % n_proc == rank) {
	std::vector<double> g(param.size());
	diag_k_l[n] = BPLAKernel<double,MData>::
	  compute_gradients(data[i], data[i], param, g);
	for (uint l=0; l!=g.size(); ++l) diag_g_l[l][n] = g[l];
	++n;
      }
    }
    // synchronization (diag)
    std::vector<double> diag_k_g(n_data/n_proc+1);
    boost::multi_array<double,2> diag_g_g(boost::extents[param.size()][n_data/n_proc+1]);
    for (uint p=0; p!=n_proc; ++p) {
      if (p==rank) {
	std::copy(diag_k_l.begin(), diag_k_l.end(), diag_k_g.begin());
	for (uint l=0; l!=diag_g_g.size(); ++l)
	  std::copy(diag_g_l[l].begin(), diag_g_l[l].end(), diag_g_g[l].begin());
      }
      MPI::COMM_WORLD.Bcast(&diag_k_g[0], diag_k_g.size(), MPI::DOUBLE, p);
      for (uint l=0; l!=diag_g_g.size(); ++l)
	MPI::COMM_WORLD.Bcast(&diag_g_g[l][0], diag_g_g[l].size(), MPI::DOUBLE, p);

      for (uint i=0, n=0; i!=data.size(); ++i) {
	if (i % n_proc == p) {
	  diag_k[i] = diag_k_g[n];
	  for (uint l=0; l!=diag_g.size(); ++l) diag_g[l][i] = diag_g_g[l][n];
	  ++n;
	}
      }
    }

    // calculation
    std::vector<double> k_l(n_cell/n_proc+1);
    boost::multi_array<double,2> g_l(boost::extents[param.size()][n_cell/n_proc+1]);
    for (uint i=0, n=0, c=0; i!=data.size()-1; ++i) {
      for (uint j=i+1; j!=data.size(); ++j) {
	if (c++ % n_proc == rank) {
	  std::vector<double> g(param.size());
	  double k;
	  k = BPLAKernel<double,MData>::
	    compute_gradients(data[i], data[j], param, g);
	  double sq_ij = sqrt(diag_k[i]*diag_k[j]);
	  k_l[n] = k/sq_ij;
	  for (uint l=0; l!=g.size(); ++l) {
	    g_l[l][n] = g[l]/sq_ij
	      - k_l[n]/2*(diag_g[l][i]/diag_k[i]+diag_g[l][j]/diag_k[j]);
	  }
	  ++n;
	}
      }
    }
    // synchronization
    std::vector<double> k_g(n_cell/n_proc+1);
    boost::multi_array<double,2> g_g(boost::extents[param.size()][n_cell/n_proc+1]);
    for (uint i=0; i!=data.size(); ++i) {
      kmat[i][i] = 1.0; // diag_k[i]/sqrt(diag_k[i]*diag_k[i]);
      for (uint l=0; l!=gmat.size(); ++l) gmat[l][i][i] = 0.0;
    }
    for (uint p=0; p!=n_proc; ++p) {
      if (p==rank) {
	std::copy(k_l.begin(), k_l.end(), k_g.begin());
	for (uint l=0; l!=g_g.size(); ++l) 
	  std::copy(g_l[l].begin(), g_l[l].end(), g_g[l].begin());
      }
      MPI::COMM_WORLD.Bcast(&k_g[0], k_g.size(), MPI::DOUBLE, p);
      for (uint l=0; l!=g_l.size(); ++l)
	MPI::COMM_WORLD.Bcast(&g_g[l][0], g_g[l].size(), MPI::DOUBLE, p);

      for (uint i=0, n=0, c=0; i!=data.size()-1; ++i) {
	for (uint j=i+1; j!=data.size(); ++j) {
	  if (c++ % n_proc == p) {
	    kmat[i][j] = kmat[j][i] = k_g[n];
	    for (uint l=0; l!=g_g.size(); ++l)
	      gmat[l][i][j] = gmat[l][j][i] = g_g[l][n];
	    ++n;
	  }
	}
      }
    }
#else
    double k;
    std::vector<double> g(param.size());
    std::vector<double> diag_k(data.size());
    boost::multi_array<double,2> diag_g(boost::extents[param.size()][data.size()]);
    for (uint i=0; i!=data.size(); ++i) {
      diag_k[i] = BPLAKernel<double,MData>::
	compute_gradients(data[i], data[i], param, g);
      for (uint l=0; l!=g.size(); ++l) diag_g[l][i] = g[l];
    }

    for (uint i=0; i!=data.size(); ++i) {
      kmat[i][i] = 1.0; // diag_k[i]/sqrt(diag_k[i]*diag_k[i]);
      for (uint l=0; l!=g.size(); ++l) gmat[l][i][i] = 0.0;
    }
    for (uint i=0; i!=data.size()-1; ++i) {
      for (uint j=i+1; j!=data.size(); ++j) {
	k = BPLAKernel<double,MData>::
	  compute_gradients(data[i], data[j], param, g);
	double sq_ij = sqrt(diag_k[i]*diag_k[j]);
	kmat[i][j] = kmat[j][i] = k/sq_ij;
	for (uint l=0; l!=g.size(); ++l) {
	  gmat[l][i][j] = gmat[l][j][i] = g[l]/sq_ij
	    - kmat[i][j]/2*(diag_g[l][i]/diag_k[i]+diag_g[l][j]/diag_k[j]);
	}
      }
    }
#endif
  }
};

template < class LDF >
static
bool
read_data(const LDF& ldf,
	  const std::vector<std::string>& labels,
	  const std::vector<std::string>& files,
	  std::vector<int>& lb,
	  std::vector<typename LDF::Data>& data)
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
      typename LDF::Loader* loader=ldf.get_loader(p->c_str());
      if (loader==NULL) return false;
      if (os)
	*os << "loading " << *p
	    << " as label " << labels[i] << std::flush;
      while (true) {
	boost::timer tm;
	typename LDF::Data* d = loader->get();
	elapsed += tm.elapsed();
	if (d==NULL) break;
	lb.push_back(atoi(labels[i].c_str()));
	data.push_back(*d);
	delete d;
      }
      if (os) *os << " (" << elapsed << "s) done." << std::endl;
      delete loader;
    }
  }
  return true;
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
#ifdef HAVE_MPI
  MPIState mpi_state(argc, argv);
  if (MPI::COMM_WORLD.Get_rank()!=0) os=NULL;
#endif  
  
  BPMatrix::Options bp_opts;
  float gap;
  float alpha;
  float beta;
  float ext;
  double C;
  uint fold;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "show this message")
    ("cost,c",
     po::value<double>(&C)->default_value(1.0),
     "set the parameter C of C-SVC")
    ("nfold,v",
     po::value<uint>(&fold)->default_value(10),
     "n-fold cross validation mode");

  po::options_description f_desc("Folding Options");
  bp_opts.add_options(f_desc);

  po::options_description k_desc("Kernel Options");
  k_desc.add_options()
    ("normalize,n", "normalize the kernel matrix")
    //("noBP", "do not use base-pairing profiles, i.e. run as local alignment kernels.")
    ("gap,g", po::value<float>(&gap)->default_value(-8.0), "set gap weight")
    ("ext,e", po::value<float>(&ext)->default_value(-0.75), "set extension weight")
    ("alpha,a", po::value<float>(&alpha)->default_value(4.5), "set alpha")
    ("beta,b", po::value<float>(&beta)->default_value(0.11), "set beta");
  
  desc.add(k_desc).add(f_desc);
  po::variables_map vm;
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).
    options(desc).allow_unregistered().run();
  std::vector<std::string> extra_args =
    collect_unrecognized(parsed.options, po::include_positional);
  std::vector<po::option>::iterator new_end =
    std::remove_if(parsed.options.begin(), parsed.options.end(),
		   bind(&po::option::unregistered, _1) );
  parsed.options.erase(new_end, parsed.options.end());
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help") || extra_args.size()<2) {
    if (os) {
      *os << "Hyperparameter optimizer for BPLA Kernels" << std::endl
	  << "Usage:" << std::endl
	  << " " << argv[0]
	  << " [options] [label1 training-data1] ... \n\n"
	  << desc << std::endl;
    }
    return 1;
  }

  // parse extra args
  std::vector<std::string> label_names;
  std::vector<std::string> files;
  {
    uint n=extra_args.size()/2;
    label_names.resize(n);
    files.resize(n);
    n=0;
    for (uint i=0; i<extra_args.size(); i+=2) {
      label_names[n]=extra_args[i];
      files[n]=extra_args[i+1];
      ++n;
    }
  }

  bool res = false;
  try {
    typedef DataLoaderFactory<DataLoader<MData> > LDF;
    LDF ldf(bp_opts);
    std::vector<int> labels;
    std::vector<MData> data;
    res = read_data(ldf, label_names, files, labels, data);
    if (!res) return 1;

    std::vector<double> param(4);
    std::vector<double> lbd(4);
    std::vector<double> ubd(4);
    std::vector<long> nbd(4, LBFGSB::UNBOUND);
    param[0] = alpha; lbd[0] = 1e-3; nbd[0] = LBFGSB::LOWER_BOUND;
    param[1] = beta; lbd[1] = 1e-3; ubd[1]= 0.3; nbd[1] = LBFGSB::LOWER_UPPER_BOUND;
    param[2] = gap; ubd[2] = 0.0; nbd[2] = LBFGSB::UPPER_BOUND;
    param[3] = ext; ubd[3] = 0.0; nbd[3] = LBFGSB::UPPER_BOUND;
    if (vm.count("normalize")) {
      Optimizer<CalcMatrixN<MData>,GradientComputationAUC> opt;
      opt.optimize(labels, data, fold, param, C, lbd, ubd, nbd,
		   1e-3, 1e10, 1e-5, os);
    } else {
      Optimizer<CalcMatrix<MData>,GradientComputationAUC> opt;
      opt.optimize(labels, data, fold, param, C, lbd, ubd, nbd,
		   1e-3, 1e10, 1e-5, os);
    }
    if (os) {
      *os << std::endl
	  << "Optimized Parameters:" << std::endl
	  << "  C=" << C << ", "
	  << "alpha=" << param[0] << ", "
	  << "beta=" << param[1] << ", "
	  << "gap=" << param[2] << ", "
	  << "ext=" << param[3] << std::endl;
    }    
    
  } catch (const char* str) {
    res = false;
    if (os) *os << str << std::endl;
  }

  return res ? 0 : 1;
}
