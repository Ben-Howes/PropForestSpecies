#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -J 1-1000

module load anaconda3/personal

cd $PBS_O_WORKDIR
Rscript bootstrap_Matrix_Model.R $PBS_ARRAY_INDEX