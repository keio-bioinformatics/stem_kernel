// $Id: pf_wrapper.cpp 48 2006-11-21 10:06:09Z satoken $

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "pf_wrapper.h"
#include <algorithm>
#include <cmath>
#include "../common/cyktable.h"
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
#include <boost/thread.hpp>
#endif

float
PFWrapper::
fold(std::string& structure)
{
#if !defined(HAVE_MPI) && defined(HAVE_BOOST_THREAD)
  static boost::mutex mtx;
  boost::mutex::scoped_lock lock(mtx);
#endif
  ::noGU = noGU_ ? 1 : 0;
  ::no_closingGU = noCloseGU_ ? 1 : 0;
  ::noLonelyPairs = noLP_ ? 1 : 0;
  float ret=0.0;
  init_pf_fold(sz_-1);
  char *s = new char[sz_];
  ret=::pf_fold(const_cast<char*>(seq_.c_str()), s);
  structure=s;
  delete[] s;
  std::copy(pr, pr+sz_*(sz_+1)/2, pr_.begin());
  std::copy(iindx, iindx+sz_, iindx_.begin());
  free_pf_arrays();
  return ret;
}

