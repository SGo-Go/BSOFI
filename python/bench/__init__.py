#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Module for analysis results of performance tuning 
#   and benchmarking.
# ------------------------------------------------------------------

__all__ = ['report', 'tuning']

from report import BenchReportTable
from tuning import BSOFIPerformanceModel, \
    BSOFISteppedBenchParameter, BSOFIFittedBenchParameter

# from plot import *
# from flops import *
