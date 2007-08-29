bin_PROGRAMS = stem_kernel_lite la_kernel

stem_kernel_lite_SOURCES = main.cpp \
	stem_kernel.cpp score_table.cpp data.cpp ribosum.cpp \
	string_kernel.cpp bpmatrix.cpp
stem_kernel_lite_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${BOOST_LIBS} 

la_kernel_SOURCES = la-main.cpp \
	data.cpp ribosum.cpp string_kernel.cpp bpmatrix.cpp
la_kernel_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${BOOST_LIBS} 
