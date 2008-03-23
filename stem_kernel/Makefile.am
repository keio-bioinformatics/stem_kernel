bin_PROGRAMS = stem_kernel #stem_train
stem_kernel_SOURCES = main.cpp stem_kernel.cpp stem_kernel.h phmm.cpp phmm.h \
	dptable.h log_value.h allocator.h
stem_kernel_LDADD = ${BOOST_LIBS} -L../common -lcommon
# stem_train_SOURCES = train.cpp stem_kernel.cpp example.cpp lbfgs.c
# stem_train_LDADD = -l${BOOST_PROGRAM_OPTIONS_LIB} -l${BOOST_THREAD_LIB} ${LAPACK_LIBS}
