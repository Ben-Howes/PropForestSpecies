#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -J 1-4

module load anaconda3/personal
source activate Renv

cd $PBS_O_WORKDIR
Rscript spatModels.r $PBS_ARRAY_INDEX
