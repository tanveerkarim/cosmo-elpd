#!/usr/bin/env python3

"""
Serial Cobaya ELPD evaluation over a subset of samples.

Designed for SLURM task-level parallelism (NOT multiprocessing).

Example usage:

python compute_lgocv_serial.py --model lcdm --start 0 --end 500
python compute_lgocv_serial.py --model lcdm --start 500 --end 1000
"""

# ============================================================
# Imports
# ============================================================

import argparse
import os
import numpy as np

from anesthetic import read_chains
from cobaya.model import get_model
from tqdm import tqdm

# ============================================================
# CLI args
# ============================================================

parser = argparse.ArgumentParser()

parser.add_argument("--model", type=str, required=True, choices=["lcdm", "w0wacdm"])
parser.add_argument("--start", type=int, required=True)
parser.add_argument("--end", type=int, required=True)

args = parser.parse_args()

MODEL_NAME = args.model
START = args.start
END = args.end

# ============================================================
# Paths
# ============================================================

PRECOMPUTED_PATH = "/project/def-rhlozek/tanveerk/cosmo-bayesian-model/data/lgo_cv_ingredients.npz"

YAML_PATHS = {
    "lcdm": "/scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-lcdm-highAc.yaml",
    "w0wacdm": "/scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-w0wacdm-db.yaml"
}

CHAIN_PATHS = {
    "lcdm": "/scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/lcdm-highacc/polychord_raw/",
    "w0wacdm": "/scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/w0waCDM-db/polychord_raw/"
}

OUTPUT_PATH = f"outputs/loglik_{MODEL_NAME}_{START}_{END}.npz"

# ============================================================
# Parameter names
# ============================================================

cosmo_params = [
    "logA", "ns", "omch2", "ombh2", "tau", "H0", "amp_143x217"
]

nuisance_params = [
    "amp_143x217", "A_planck", "n_217", "n_143",
    "amp_143", "calEE", "calTE", "n_143x217", "amp_217"
]

param_names = np.concatenate((cosmo_params, nuisance_params))

# ============================================================
# Gaussian helper
# ============================================================

def gaussian_logpdf(delta, Sigma_inv, norm):
    return -0.5 * (delta @ Sigma_inv @ delta + norm)

# ============================================================
# Load precomputed data
# ============================================================

print("Loading precomputed ingredients...")

precomputed = np.load(PRECOMPUTED_PATH, allow_pickle=True)

y_joint = precomputed["y_joint"]
q_joint = precomputed["q_joint"]
chunk_indices = precomputed["chunk_indices"][()]
sigma_cond_list = precomputed["sigma_cond_list"]

K = len(chunk_indices)

print(f"Chunks: {K}")

# ============================================================
# Precompute Gaussian matrices
# ============================================================

print("Precomputing Gaussian matrices...")

gaussian_cache = []

for k, Sigma_cond in enumerate(sigma_cond_list):

    sign, logdet = np.linalg.slogdet(Sigma_cond)

    if sign <= 0:
        raise ValueError(f"Chunk {k} covariance not positive definite")

    Sigma_inv = np.linalg.inv(Sigma_cond)

    N = Sigma_cond.shape[0]
    norm = logdet + N * np.log(2.0 * np.pi)

    gaussian_cache.append({
        "Sigma_inv": Sigma_inv,
        "norm": norm
    })

print("Gaussian precomputation complete")

# ============================================================
# Load chain
# ============================================================

print(f"Loading chain: {MODEL_NAME}")

samples_df = read_chains(CHAIN_PATHS[MODEL_NAME]).posterior_points()

TOTAL = len(samples_df)

END = min(END, TOTAL)
S = END - START

print(f"Running samples [{START}, {END}) => {S} samples")

# ============================================================
# Initialize Cobaya
# ============================================================

print("Initializing Cobaya model...")

model = get_model(YAML_PATHS[MODEL_NAME])

like_sn = model.likelihood["sn.pantheonplus"]
like_bao = model.likelihood["bao.desi_dr2"]

print("Cobaya ready")

# ============================================================
# Output
# ============================================================

log_lik_matrix = np.zeros((1, S, K))

# ============================================================
# Main loop (SERIAL)
# ============================================================

print("Starting evaluation...")

for i, s in enumerate(tqdm(range(START, END))):

    row = samples_df.iloc[s]

    theta = row.droplevel(1)[param_names].to_dict()

    model.logposterior(theta)

    provider = model.provider

    # SN
    da_sn = provider.get_angular_diameter_distance(like_sn.zcmb)

    mu_sn = 5.0 * np.log10(
        (1.0 + like_sn.zhel)
        * (1.0 + like_sn.zcmb)
        * da_sn
    )

    # BAO
    mu_bao = np.array([
        like_bao.theory_fun(z, obs)
        for z, obs in zip(
            like_bao.data["z"],
            like_bao.data["observable"]
        )
    ]).T[0]

    mu_joint = np.concatenate([mu_sn, mu_bao])

    # residual
    g = q_joint @ (y_joint - mu_joint)

    logL_chunks = np.zeros(K)

    for k, idx in enumerate(chunk_indices):

        Sigma_cond = sigma_cond_list[k]
        cache = gaussian_cache[k]

        Sigma_inv = cache["Sigma_inv"]
        norm = cache["norm"]

        y_I = y_joint[idx]
        g_I = g[idx]

        mu_cond = y_I - Sigma_cond @ g_I

        delta = y_I - mu_cond

        logL_chunks[k] = gaussian_logpdf(delta, Sigma_inv, norm)

    log_lik_matrix[0, i, :] = logL_chunks

# ============================================================
# Save
# ============================================================

os.makedirs("outputs", exist_ok=True)

out_file = OUTPUT_PATH

np.savez_compressed(
    out_file,
    log_lik_matrix=log_lik_matrix,
    start=START,
    end=END,
    model=MODEL_NAME
)

print(f"Saved: {out_file}")
print("DONE")
