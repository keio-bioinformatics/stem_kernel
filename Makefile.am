if OPT_REF_IMPL
  STEM_KERNEL = stem_kernel
endif
if OPT_STRING_KERNEL
  STRING_KERNEL = string_kernel
endif
if OPT_SIMPAL
  SIMPAL = simpal
endif
SUBDIRS =  common libsvm optimizer \
	$(STRING_KERNEL) $(STEM_KERNEL) $(SIMPAL) stem_kernel_lite bpla_kernel

EXTRA_DIST = README.rd README.html

README.html: README.rd
	rd2 README.rd > README.html
