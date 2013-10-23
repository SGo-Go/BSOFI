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

L = 4
n_list = (1024, 512, 256)
nthreads = 12
implementations = ('hybrid',)# 'cpu')
platform = 'dirac_fermi'

ofolder = r"./results/%s" % platform
cmd_all = r"./bench/{implement}/bench -t {ntests} -n {blksize} -L {blknum}"\
    " -f {ofolder}/{blksize:04d}_{nthreads:02d}.blk"
scmd_all = r"./bench/%s/bench -t %s -n %s -L %s -f %s/%04i_%02i.blk"

if not os.path.exists(ofolder): os.makedirs(ofolder)
for implementation in implementations:
    ofolder_local = os.path.join(ofolder, implementation)
    if not os.path.exists(ofolder_local): os.mkdir(ofolder_local)
    # for nthreads in (12,): # 6, 4, 2)
    # os.environ['OMP_NUM_THREADS'] = str(nthreads)
    # os.environ['MKL_NUM_THREADS'] = str(nthreads)

    for n in n_list:
        ntests = 44*(1024/n) - 3
        # cmd = cmd_all.format(implement = implementation, 
        #                      ofolder = ofolder_local,
        #                      nthreads=nthreads, ntests=ntests, 
        #                      blknum=L, blksize=n)
        cmd = scmd_all % (implementation, ntests, n, L, 
                          ofolder_local, n, nthreads)
        with open(r"config/jobs/%s.sh.in" % platform, "r") as job_script_handle:
            job_script = job_script_handle.read()
        job_script = job_script.replace(r"@cmd@", cmd)
        print job_script
        job_file = r"%s/%s%04i.job" % (ofolder_local, platform,n)
        with open(job_file, "w") as job_script_handle:
            job_script_handle.write(job_script)
        os.system(r"qsub -S /bin/bash %s" % job_file)
