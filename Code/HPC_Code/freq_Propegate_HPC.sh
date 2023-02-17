#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=64:mem=128gb

module load anaconda3/personal

cd $PBS_O_WORKDIR
Rscript freq_Propegate_Matrix_Species.R
