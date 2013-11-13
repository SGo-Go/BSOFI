#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko/ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Script for submitting job
# ------------------------------------------------------------------

import os, sys
os.environ['BSOFI_TOPDIR'] = r"/home/sergiy/project/BSOFI"

sys.path.append(os.path.join(os.getenv('BSOFI_TOPDIR', '.'), 'python'))

from bench import BSOFIPerformanceModel, BenchReportTable
from platforms import user as platform
from platforms import Job
#Job, Platform, BashJobManager, PBSJobManager #zwolf as platform

manager = platform.job_manager
manager.submit(Job(" ".join(sys.argv[1:])))