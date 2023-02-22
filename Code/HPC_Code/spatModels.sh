#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=9:mem=920gb
#PBS -J 1-4

module load anaconda3/personal

cd $PBS_O_WORKDIR
Rscript spatModels.r $PBS_ARRAY_INDEX
