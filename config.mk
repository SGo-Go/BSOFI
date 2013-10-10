# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenk/ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Template for config settings.

SYSTEM  = zwolf
TOP_DIR = /home/sergiy/project/BSOFI/trunk

include $(TOP_DIR)/config/platforms/$(SYSTEM)-settings.mk
include $(TOP_DIR)/config/makefile-general-settings.mk
