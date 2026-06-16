#!/bin/bash
#SBATCH --account=def-rhlozek
#SBATCH --job-name=w0wa-pc-debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4       # Request all 192 cores on a Trillium node
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --partition=debug
#SBATCH --output=%x-%j.out          # Standard output log (JobName-JobID.out)
#SBATCH --mail-user=tanveer.karim@utoronto.ca
#SBATCH --mail-type=ALL

# 1. Purge modules and load the exact environment stack
module purge
module load StdEnv/2023 gcc openmpi python/3.13 scipy-stack mpi4py

# 2. Activate your virtual environment
source /project/def-rhlozek/tanveerk/vEnvs/stats/bin/activate

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 3. Critical environment overrides (pointing to the Scratch cache)
export COBAYA_PACKAGES_PATH=/scratch/tanveerk/cobaya-external
export COBAYA_USE_FILE_LOCKING=false

# 5. Launch the job using pure Python module execution
echo "Starting PolyChord run at $(date)"
#mpirun -n 4 python -m cobaya run /scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-w0wacdm.yaml
srun python -m cobaya run /scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-w0wacdm.yaml
echo "Finished at $(date)"
