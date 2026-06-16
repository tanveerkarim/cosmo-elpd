#!/bin/bash

#SBATCH --job-name=lgo-cv
#SBATCH --output=logs/elpd_%j.out
#SBATCH --error=logs/elpd_%j.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

#SBATCH --time=01:00:00

#SBATCH --account=def-rhlozek
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=tanveer.karim@utoronto.ca
#SBATCH --mail-type=ALL

# 1. Purge modules and load the exact environment stack
module purge
module load StdEnv/2023 gcc openmpi python/3.13 scipy-stack mpi4py

# 2. Activate your virtual environment
source /project/def-rhlozek/tanveerk/vEnvs/stats/bin/activate

# 3. Thread limiters (Critical for pure MPI)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 4. Critical environment overrides
export COBAYA_PACKAGES_PATH=/scratch/tanveerk/cobaya-external
export COBAYA_USE_FILE_LOCKING=false

source ~/vEnvs/stats/bin/activate

python compute_lgocv.py
