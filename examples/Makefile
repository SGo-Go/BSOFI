# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Examples makefile.
# ------------------------------------------------------------------

TOP_DIR = ..
include ../config.mk

.PHONY  : clean install all lib examples ctest ftest_hybrid ftest_cpu ctest_gpu

all: ctest ftest_cpu ftest_hybrid ctest_gpu

ctest:
	($(CD) ./ctest; $(MAKE))

ctest_gpu:
	($(CD) ./ctest_gpu; $(MAKE))

ftest_hybrid:
	($(CD) ./ftest_hybrid; $(MAKE))

ftest_cpu:
	($(CD) ./ftest_cpu; $(MAKE))

clean:
	($(CD) ./ctest; $(MAKE) clean)
	($(CD) ./ctest_gpu; $(MAKE) clean)
	($(CD) ./ftest_hybrid; $(MAKE) clean)
	($(CD) ./ftest_cpu; $(MAKE) clean)
	$(RM) ./*~
