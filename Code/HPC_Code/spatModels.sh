#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -J 1-400

module load anaconda3/personal

cd $PBS_O_WORKDIR
Rscript spatModels.r $PBS_ARRAY_INDEX
