bin_PROGRAMS = stem_kernel_lite pp_stem_kernel_lite la_kernel pp_la_kernel

stem_kernel_lite_SOURCES = main.cpp \
	stem_kernel.cpp score_table.cpp data.cpp ribosum.cpp \
	string_kernel.cpp bpmatrix.cpp dag.cpp
stem_kernel_lite_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${BOOST_LIBS} 

pp_stem_kernel_lite_SOURCES = pp-main.cpp \
	stem_kernel.cpp score_table.cpp data.cpp ribosum.cpp \
	string_kernel.cpp bpmatrix.cpp dag.cpp
pp_stem_kernel_lite_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${BOOST_LIBS} 

la_kernel_SOURCES = la-main.cpp \
	data.cpp ribosum.cpp string_kernel.cpp bpmatrix.cpp dag.cpp
la_kernel_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${BOOST_LIBS} 

pp_la_kernel_SOURCES = pp-la-main.cpp \
	data.cpp ribosum.cpp string_kernel.cpp bpmatrix.cpp dag.cpp
pp_la_kernel_LDADD = ../common/libcommon.a ../libsvm/libsvm.a ${BOOST_LIBS} 
