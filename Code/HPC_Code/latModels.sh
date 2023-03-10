#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -J 1-4000

module load anaconda3/personal

cd $PBS_O_WORKDIR
Rscript latModels.r $PBS_ARRAY_INDEX
