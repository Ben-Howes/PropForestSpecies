#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=64gb

module load anaconda3/personal

cd $PBS_O_WORKDIR
Rscript spaMM_models.r
