#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Benchmarking script.
# ------------------------------------------------------------------

import os

L = 4
n_list = (1024, 512, 256)
nthreads = 12
implementations = ('hybrid', 'cpu')

ofolder = r"./results"
cmd_all = r"./{implement}/bench -t {ntests} -n {blksize} -L {blknum} > "\
    "{ofolder}/mkl_{blksize:04d}_{nthreads:02d}.blk"

if not os.path.exists(ofolder): os.mkdir(ofolder)
for implementation in implementations:
    ofolder_local = os.path.join(ofolder, implementation)
    if not os.path.exists(ofolder_local): os.mkdir(ofolder_local)
    # for nthreads in (12,): # 6, 4, 2)
    # os.environ['OMP_NUM_THREADS'] = str(nthreads)
    # os.environ['MKL_NUM_THREADS'] = str(nthreads)

    for n in n_list:
        ntests = 4*(1024/n) - 3
        print     cmd_all.format(implement = implementation, ofolder = ofolder_local,
                                 nthreads=nthreads, ntests=ntests, blknum=L, blksize=n)
        os.system(cmd_all.format(implement = implementation, ofolder = ofolder_local,
                                 nthreads=nthreads, ntests=ntests, blknum=L, blksize=n))
