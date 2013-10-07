# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk/ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Auto-selection of configuration.

SYSTEM  = user
TOP_DIR = /home/sergiy/project/BSOFI

include $(TOP_DIR)/config/platforms/$(SYSTEM)-settings.mk
include $(TOP_DIR)/config/makefile-general-settings.mk
