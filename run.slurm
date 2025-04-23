#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --job-name="alt-sr"
#SBATCH --mem=30GB
# #SBATCH --nodes=4
# #SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --output=alt-sr.out

cd /home/s.defina/comorb-longit-project

module load R

Rscript 1.ALTSR.R '../alt_sr_results/' '../longi/raw_data.rds'

# Rscript 1.RICLPM.R './results/ri_pfree_transf'
# Rscript 1.RICLPM.R './results/ri_pstat_transf'
# Rscript 1.RICLPM.R './results/ri_pstat_std'
# Rscript 1.RICLPM.R './results/ri_pfree_std'

# Rscript 1.GCLPM.R './results/g_lfree_pstat'
# Rscript 1.GCLPM.R './results/g_lstat_pstat'
# Rscript 1.GCLPM.R './results/g_lstat_pfree'
# Rscript 1.GCLPM.R './results/g_lfree_pfree'
