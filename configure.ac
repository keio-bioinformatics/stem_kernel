#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(stem_kernel, 0.1, no-define)
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([stem_kernel/main.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_CXX_HAVE_STL
AC_PROG_RANLIB

# Checks for libraries.
#ACX_MPI

dnl AC_ARG_WITH([mpi], AS_HELP_STRING([--with-mpi],[use MPICH]),
dnl    [ mpilib=$withval
dnl      if test "x$mpilib" = "xyes"; then
dnl        mpilib=mpi++
dnl      fi
dnl      AC_CHECK_LIB($mpilib, main)
dnl      AC_CHECK_HEADER(mpi.h, AC_DEFINE([HAVE_MPI], [], [have MPI headers])) 
dnl    ])
AC_ARG_WITH([mpi], AS_HELP_STRING([--with-mpi],[use MPICH]),
   [ AC_CHECK_HEADER(mpi.h, AC_DEFINE([HAVE_MPI], [], [have MPI headers])) ])

AX_BOOST_PROGRAM_OPTIONS
if test "$ax_cv_boost_program_options" = no; then
  AC_MSG_FAILURE("the Boost::Program_Options library is unavailable")  
else
  BOOST_LIBS="$BOOST_LIBS -l$BOOST_PROGRAM_OPTIONS_LIB"
fi
AX_BOOST_THREAD
if test "x$BOOST_THREAD_LIB" != "x"; then
  BOOST_LIBS="$BOOST_LIBS -l$BOOST_THREAD_LIB"
fi
AX_BOOST_PROGRAM_OPTIONS
if test "$ax_cv_boost_program_options" = no; then
  AC_MSG_FAILURE("the Boost::Program_Options library is unavailable")  
else
  BOOST_LIBS="$BOOST_LIBS -l$BOOST_PROGRAM_OPTIONS_LIB"
fi
AC_SUBST(BOOST_LIBS)

AC_CHECK_LIB(RNA, init_pf_fold)

# Checks for header files.
AC_CHECK_HEADER(boost/multi_array.hpp, ,
                AC_MSG_FAILURE("boost/multi_array.hpp is unavailable"))

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])
#ACX_LAPACK
AC_CONFIG_FILES([Makefile common/Makefile libsvm/Makefile string_kernel/Makefile stem_kernel/Makefile simpal/Makefile stem_kernel_lite/Makefile rla_kernel/Makefile])
AC_OUTPUT
