#!/bin/bash
#SBATCH --account=def-rhlozek
#SBATCH --job-name=mpi-192-w0wa-highAcc-pc
#SBATCH --nodes=1
#SBATCH --ntasks=192              # Changed from ntasks-per-node to explicitly request 192 tasks
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00            # Set to 30 mins for the debug trial
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

echo "=========================================="
echo "Job started at $(date)"
echo "Running on node(s): $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="

# -------------------------------------------------
# Optional cleanup of stale lock files
# Only delete locks older than 30 minutes
# -------------------------------------------------

echo "Sweeping stale lock files..."

find /scratch/tanveerk/bayesian-model-workspace/ \
    -name "*.locked" -type f -delete
find /scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/w0waCDM/ -name "*.locked" -type f -delete
find /scratch/tanveerk/cobaya-external/data/planck_2018_CamSpec2021/ \
    -name "*.locked" -type f -delete

echo "Cleanup complete."

# -------------------------------------------------
# Optional serial warm-up test
# Helps populate caches before MPI launch
# -------------------------------------------------

#echo "Running serial warm-up test..."

#python -m cobaya run \
#    /scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-w0wacdm.yaml \
#    --test

#WARMUP_EXIT=$?

#if [ $WARMUP_EXIT -ne 0 ]; then
#    echo "Warm-up test failed with exit code $WARMUP_EXIT"
#    exit $WARMUP_EXIT
#fi

#echo "Warm-up test completed successfully."

##----------
# Clean warm-up lock files
# ----------
#echo "Cleaning temporary warm-up output..."

#rm -rf /scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/w0waCDM/polychord_raw
#rm -rf /scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/w0waCDM/clusters

#find /scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/w0waCDM/ -name "*.locked" -type f -delete
#find /scratch/tanveerk/bayesian-model-workspace/ -name "*.locked" -type f -delete
#echo "Warm-up cleanup complete."


# -------------------------------------------------
# Launch MPI production run
# -------------------------------------------------
# 6. Launch the job using native SLURM srun
echo "Starting PolyChord run at $(date)"
srun python -m cobaya run /scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-w0wacdm.yaml

RUN_EXIT=$?

echo "=========================================="
echo "Run finished at $(date)"
echo "Exit code: $RUN_EXIT"
echo "=========================================="

exit $RUN_EXIT
