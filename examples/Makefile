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

SUBFOLDERS  = ctest ctest_gpu

ifneq "$(FC)" ""
SUBFOLDERS += ftest_cpu ftest_hybrid 
endif

.PHONY  : clean install all lib examples $(SUBFOLDERS)

all: $(SUBFOLDERS)

ctest:
	($(CD) ./ctest; $(MAKE))

ctest_gpu:
	($(CD) ./ctest_gpu; $(MAKE))

ifneq "$(FC)" ""
ftest_hybrid:
	($(CD) ./ftest_hybrid; $(MAKE))

ftest_cpu:
	($(CD) ./ftest_cpu; $(MAKE))
endif

clean:
	($(CD) ./ctest; $(MAKE) clean)
	($(CD) ./ctest_gpu; $(MAKE) clean)
ifneq "$(FC)" ""
	($(CD) ./ftest_hybrid; $(MAKE) clean)
	($(CD) ./ftest_cpu; $(MAKE) clean)
endif
	$(RM) ./*~
