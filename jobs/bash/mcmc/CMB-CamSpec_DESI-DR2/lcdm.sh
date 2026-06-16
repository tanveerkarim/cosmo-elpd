#!/bin/bash
#SBATCH --account=def-rhlozek
#SBATCH --job-name=mc-d+c-l
#SBATCH --nodes=1
#SBATCH --ntasks=32              # Changed from ntasks-per-node to explicitly request 192 tasks
#SBATCH --cpus-per-task=6
#SBATCH --time=03:30:00            # Set to 30 mins for the debug trial
#SBATCH --output=/scratch/tanveerk/bayesian-model-workspace/slurm-log/20260616/%x-%j.out
#SBATCH --mail-user=tanveer.karim@utoronto.ca
#SBATCH --mail-type=ALL

# 1. Purge modules and load the exact environment stack
module purge
module load StdEnv/2023 gcc openmpi python/3.13 scipy-stack mpi4py

# 2. Activate your virtual environment
source /project/def-rhlozek/tanveerk/vEnvs/stats/bin/activate

# 3. Thread limiters (Critical for pure MPI)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 4. Critical environment overrides
export COBAYA_PACKAGES_PATH="/scratch/tanveerk/cobaya-external"
export COBAYA_USE_FILE_LOCKING=false

model="lcdm"
combo="CMB-CamSpec_DESI-DR2"

echo "=========================================="
echo "Job started at $(date)"
echo "Running on node(s): $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "Model: $model"
echo "Combo: $combo"
echo "=========================================="

# -------------------------------------------------
# Optional cleanup of stale lock files
# -------------------------------------------------

echo "Sweeping stale lock files..."

find "/scratch/tanveerk/bayesian-model-workspace/chains/mcmc/${combo}/${model}/" -name "*.locked" -type f -delete
find "/scratch/tanveerk/cobaya-external/data/planck_2018_CamSpec2021/" -name "*.locked" -type f -delete
echo "Cleanup complete."

# -------------------------------------------------
# Launch MPI production run
# -------------------------------------------------
# 6. Launch the job using native SLURM srun
echo "Starting PolyChord run at $(date)"
srun python -m cobaya run "/scratch/tanveerk/bayesian-model-workspace/yaml/${combo}/${model}.yaml"

RUN_EXIT=$?

echo "=========================================="
echo "Run finished at $(date)"
echo "Exit code: $RUN_EXIT"
echo "=========================================="

exit $RUN_EXIT
