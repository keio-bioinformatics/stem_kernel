#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(stem_kernel, 296, no-define)
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([stem_kernel/main.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_CXX_HAVE_STL
AC_PROG_RANLIB
AM_PROG_AR

AC_ARG_ENABLE([ref-impl],
  AS_HELP_STRING([--enable-ref-impl],
    [build naive reference implimentation of stem kernels]))

AC_ARG_ENABLE([string-kernel],
  AS_HELP_STRING([--enable-string-kernel],
    [build naive string kernel kernels]))

AC_ARG_ENABLE([simpal],
  AS_HELP_STRING([--enable-simpal],
    [build simple palindrom kernels]))

# CONTRAfold
AC_ARG_WITH([contrafold], AS_HELP_STRING([--with-contrafold],[link the CONTRAfold library]),
    [
      AC_DEFINE([HAVE_LIBCONTRAFOLD], [], [link the CONTRAfold library])
      LIBCONTRAFOLD="-lcontrafold"
    ]
)
AC_SUBST(LIBCONTRAFOLD)
	    	   
# Checks for libraries.

AC_ARG_WITH([mpi], AS_HELP_STRING([--with-mpi],[use MPICH]),
   [ ACX_MPI([
       AC_DEFINE([HAVE_MPI], [], [have MPI headers])
       #CC=$MPICXX
       CXX=$MPICXX
   ]) 
])        
				    
AX_BOOST_PROGRAM_OPTIONS
if test "x$BOOST_PROGRAM_OPTIONS_LIB" != x; then
  BOOST_LIBS="$BOOST_LIBS $BOOST_PROGRAM_OPTIONS_LIB"
else
  AC_MSG_ERROR("the Boost::Program_Options library is unavailable")  
fi
AX_BOOST_SYSTEM
if test "x$BOOST_SYSTEM_LIB" != x; then
  BOOST_LIBS="$BOOST_LIBS $BOOST_SYSTEM_LIB"
else
  AC_MSG_ERROR("the Boost::System library is unavailable")  
fi
AX_BOOST_THREAD
if test "x$BOOST_THREAD_LIB" != "x"; then
  BOOST_LIBS="$BOOST_LIBS $BOOST_THREAD_LIB"
fi
AX_BOOST_IOSTREAMS
if test "x$BOOST_IOSTREAMS_LIB" != "x"; then
  BOOST_LIBS="$BOOST_LIBS $BOOST_IOSTREAMS_LIB"
fi
AC_SUBST(BOOST_LIBS)

AC_CHECK_LIB(RNA, init_pf_fold)

# Checks for header files.
AC_CHECK_HEADER([boost/multi_array.hpp], ,
                [AC_MSG_ERROR("boost/multi_array.hpp is unavailable")])

AC_CHECK_HEADER([ViennaRNA/ribo.h],
		[AC_DEFINE([HAVE_VIENNA18], [], [have the Vienna RNA package version 1.8 or above])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AM_CONDITIONAL([OPT_REF_IMPL], [test "$enable_ref_impl" = yes])
AM_CONDITIONAL([OPT_STRING_KERNEL], [test "$enable_string_kernel" = yes])
AM_CONDITIONAL([OPT_SIMPAL], [test "$enable_simpal" = yes])

AC_CONFIG_FILES([Makefile common/Makefile libsvm/Makefile optimizer/Makefile string_kernel/Makefile stem_kernel/Makefile simpal/Makefile stem_kernel_lite/Makefile bpla_kernel/Makefile])
AC_OUTPUT
