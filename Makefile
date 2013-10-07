# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Top makefile.
# ------------------------------------------------------------------

include ./config.mk

.PHONY  : clean install all lib examples

all: lib examples

lib:
	($(CD) src; $(MAKE))

examples:
	($(CD) examples; $(MAKE))

clean:
	$(RM) $(TOP_DIR)/*~
	($(CD) $(TOP_DIR)/src; $(MAKE) clean)
	($(CD) $(TOP_DIR)/examples; $(MAKE) clean)

allclean: clean
	$(RM) -R `find $(TOP_DIR) -name *~` $(TOP_DIR)/*.o
