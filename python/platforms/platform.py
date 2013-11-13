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

from job import Job, PBSJobManager, BashJobManager


class Platform(object):
    def __init__(self, job_manager, **kw):
        super(Platform, self).__init__()
        self.job_manager = job_manager
        self.attrs = kw
        pass
