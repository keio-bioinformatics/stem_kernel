bin_PROGRAMS = stem_kernel_lite la_kernel

stem_kernel_lite_SOURCES = main.cpp stem_kernel.cpp stem_kernel.h \
	score_table.cpp score_table.h data.cpp data.h \
	ribosum.cpp ribosum.h string_kernel.cpp string_kernel.h \
	def_kernel.h dptable.h dag.h
stem_kernel_lite_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${LIBCONTRAFOLD} ${BOOST_LIBS} 

la_kernel_SOURCES = la-main.cpp string_kernel.cpp string_kernel.h \
	data.cpp data.h ribosum.cpp ribosum.h dptable.h def_kernel.h
la_kernel_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${LIBCONTRAFOLD} ${BOOST_LIBS} 
