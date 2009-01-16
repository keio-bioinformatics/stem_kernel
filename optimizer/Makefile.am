bin_PROGRAMS = rbf_optimizer poly_optimizer sigmoid_optimizer
noinst_LIBRARIES = liboptimizer.a

liboptimizer_a_SOURCES = gradient.cpp gradient.h optimizer.cpp optimizer.h lbfgsb.c lbfgsb.h 

rbf_optimizer_SOURCES = rbf_optimizer.cpp rbf_kernel.cpp rbf_kernel.h
rbf_optimizer_LDADD = liboptimizer.a ../libsvm/libsvm.a ${BOOST_LIBS}

poly_optimizer_SOURCES = poly_optimizer.cpp poly_kernel.cpp poly_kernel.h
poly_optimizer_LDADD = liboptimizer.a ../libsvm/libsvm.a ${BOOST_LIBS}

sigmoid_optimizer_SOURCES = sigmoid_optimizer.cpp sigmoid_kernel.cpp sigmoid_kernel.h
sigmoid_optimizer_LDADD = liboptimizer.a ../libsvm/libsvm.a ${BOOST_LIBS}
