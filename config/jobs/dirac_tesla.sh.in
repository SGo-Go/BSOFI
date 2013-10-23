#PBS -q dirac_small
#PBS -l nodes=1:ppn=8:tesla
#PBS -l walltime=02:00:00
#PBS -N my_job
#PBS -e job.$PBS_JOBID.err
#PBS -o job.$PBS_JOBID.out
#PBS -V

module load cuda
module unload pgi
module load intel
module load mkl
module load python

cd $PBS_O_WORKDIR
# ${PBS_O_WORKDIR}/tune/hybrid/tune -f 32 -l 2048 -s 32
@cmd@