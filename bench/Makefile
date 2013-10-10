# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Benchmarks makefile.
# ------------------------------------------------------------------

TOP_DIR = ..
include ../config.mk

.PHONY  : clean install all lib examples bench cpu hybrid

all: cpu hybrid

cpu:
	($(CD) ./cpu; $(MAKE))

hybrid:
	($(CD) ./hybrid; $(MAKE))

clean:
	($(CD) ./cpu; $(MAKE) clean)
	($(CD) ./hybrid; $(MAKE) clean)
	$(RM) ./*~

allclean: clean
	$(RM) -r ./results
