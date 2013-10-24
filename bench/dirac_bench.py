#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Benchmarking script for Dirac (temporary solution).
# ------------------------------------------------------------------

import os

############################################################
# Parameters
############################################################

L    = 4
Lmax = 44
n_list = (1024, 512, 256)
nthreads = 12
implementations = ('hybrid',)# 'cpu')
platform = 'dirac_fermi'

ofolder = r"./results/%s" % platform

############################################################
def call_minimalistic():
    """
    Implementation of job running 
    with minimalistic requirements to Python interpreter
    """
    cmd_format          = r"./bench/%s/bench -t %s -n %s -L %s -f %s/%04i_%02i.blk"
    job_file_format     = r"%s/%s%04i.job"
    job_template_format = r"config/jobs/%s.sh.in"

    ntests = Lmax*(1024/n) - 3
    cmd = cmd_format % (implementation, ntests, n, L, 
                     ofolder_local, n, nthreads)
    job_template_handle = open(job_template_format % platform, "r")
    job_script = job_template_handle.read()
    job_template_handle.close()
    job_script = job_script.replace(r"@cmd@", cmd)
    # print job_script
    job_file = job_file_format % (ofolder_local, platform,n)
    job_script_handle = open(job_file, "w")
    job_script_handle.write(job_script)
    job_script_handle.close()
    os.system(r"qsub -S /bin/bash %s" % job_file)

# ############################################################
# def call_modern():
#     """
#     Implementation of job running 
#     for Python version >= 2.5
#     """
#     cmd_format = r"./bench/{implement}/bench -t {ntests} -n {blksize} -L {blknum}" \
#         " -f {ofolder}/{blksize:04d}_{nthreads:02d}.blk"
#     job_file_format = r"{ofolder}/{platform}{blksize:04d}.job"
#     job_template_format = r"config/jobs/{platform}.sh.in"

#     ntests = Lmax*(1024/n) - 3
#     cmd = cmd_format.format(implement = implementation, 
#                          ofolder = ofolder_local,
#                          nthreads=nthreads, ntests=ntests, 
#                          blknum=L, blksize=n)
#     with open(job_template_format.format(platform=platform), "r") as job_template_handle:
#         job_script = job_template_handle.read()
#     job_script = job_script.replace(r"@cmd@", cmd)
#     print job_script
#     job_file = job_file_format.format(ofolder = ofolder_local, 
#                                       platform = platform, blksize = n)
#     with open(job_file, "w") as job_script_handle:
#         job_script_handle.write(job_script)
#     os.system(r"qsub -S /bin/bash %s" % job_file)

############################################################ 
# Cycle over commands to run
############################################################ 
if not os.path.exists(ofolder): os.makedirs(ofolder)
for implementation in implementations:
    ofolder_local = os.path.join(ofolder, implementation)
    if not os.path.exists(ofolder_local): os.mkdir(ofolder_local)
    # for nthreads in (12,): # 6, 4, 2)
    # os.environ['OMP_NUM_THREADS'] = str(nthreads)
    # os.environ['MKL_NUM_THREADS'] = str(nthreads)

    for n in n_list:
        call_minimalistic() # call_modern() #
