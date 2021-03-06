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

.PHONY  : clean install all lib examples bench tune

all: lib examples bench tune

lib:
	($(CD) src; $(MAKE))

examples: lib
	($(CD) examples; $(MAKE))

bench: lib
	($(CD) bench; $(MAKE))

tune: lib
	($(CD) tune; $(MAKE))

clean:
	$(RM) $(TOP_DIR)/*~
	($(CD) $(TOP_DIR)/src; $(MAKE) clean)
	($(CD) $(TOP_DIR)/examples; $(MAKE) clean)
	($(CD) $(TOP_DIR)/bench; $(MAKE) clean)
	($(CD) $(TOP_DIR)/tune; $(MAKE) clean)

allclean: clean
	($(CD) $(TOP_DIR)/bench; $(MAKE) allclean)
	$(RM) -R `find $(TOP_DIR) -name *~` $(TOP_DIR)/*.o
	$(RM) -R `find $(TOP_DIR) -name *.pyc`
	#$(RM) $(TOP_DIR)/config.mk
