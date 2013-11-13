#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Module for running jobs on different platforms.
# ------------------------------------------------------------------

__all__ = ['job', 'platform']

from job import Job, PBSJobManager, BashJobManager
from platform import Platform

user        = Platform(BashJobManager())
zwolf       = Platform(BashJobManager())
dirac_fermi = Platform(PBSJobManager(r'config/jobs/dirac_fermi.sh.in'))
dirac_tesla = Platform(PBSJobManager(r'config/jobs/dirac_tesla.sh.in'))
