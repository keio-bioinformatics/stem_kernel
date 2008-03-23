bin_PROGRAMS = bpla_kernel rla_kernel

bpla_kernel_SOURCES = main.cpp bpla_kernel.cpp bpla_kernel.h \
	data.cpp data.h bpmatrix.cpp bpmatrix.h
bpla_kernel_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${BOOST_LIBS}

rla_kernel_SOURCES = rla-main.cpp string_kernel.cpp string_kernel.h
rla_kernel_LDADD = ../common/libcommon.a ${BOOST_LIBS}
