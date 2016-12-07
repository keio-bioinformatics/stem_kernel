bin_PROGRAMS = bpla_kernel la_kernel bpla_optimizer

bpla_kernel_SOURCES = main.cpp bpla_kernel.cpp bpla_kernel.h data.cpp data.h
bpla_kernel_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${LIBCONTRAFOLD} ${BOOST_LIBS}

la_kernel_SOURCES = la_main.cpp bpla_kernel.cpp bpla_kernel.h data.cpp data.h
la_kernel_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${LIBCONTRAFOLD} ${BOOST_LIBS}

bpla_optimizer_SOURCES = bpla_optimizer.cpp bpla_kernel.cpp bpla_kernel.h data.cpp data.h
bpla_optimizer_LDADD = ../common/libcommon.a ../optimizer/liboptimizer.a ../libsvm/libsvm.a \
	${LIBCONTRAFOLD} ${BOOST_LIBS}
