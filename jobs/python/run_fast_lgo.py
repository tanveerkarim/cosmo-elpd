import numpy as np
import pandas as pd
import multiprocessing as mp
from cobaya.model import get_model
from scipy.stats import multivariate_normal
from tqdm import tqdm
from anesthetic import read_chains


# --- GLOBALS: Loaded once in the parent process ---
# These will be inherited by workers via copy-on-write.
MODEL_PATH = "/scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-lcdm-highAcc.yaml"
INGREDIENTS_PATH = "/project/def-rhlozek/tanveerk/cosmo-bayesian-model/data/lgo_cv_ingredients.npz"

def init_worker():
    global model, ingredients, like_sn, like_bao
    # Reload ingredients per-worker (small, fast array)
    ingredients = np.load(INGREDIENTS_PATH, allow_pickle=True)
    # The 'model' is already inherited from parent process memory
    like_sn = model.likelihood['sn.pantheonplus']
    like_bao = model.likelihood['bao.desi_dr2']

def worker_eval(args):
    s, theta_dict = args
    # Update cosmology and compute
    model.logposterior(theta_dict)
    provider = model.provider
    
    # Compute theory
    da_sn = provider.get_angular_diameter_distance(like_sn.zcmb)
    mu_sn = 5.0 * np.log10((1.0 + like_sn.zhel) * (1.0 + like_sn.zcmb) * da_sn)
    mu_bao = np.array([like_bao.theory_fun(z, obs) for z, obs in zip(like_bao.data["z"], like_bao.data["observable"])]).T[0]
    mu_joint = np.concatenate([mu_sn, mu_bao])
    
    # Calculate LogL
    y_joint = ingredients['y_joint']
    q_joint = ingredients['q_joint']
    chunk_indices = ingredients['chunk_indices']
    sigma_cond_list = ingredients['sigma_cond_list']
    
    logL_chunks = np.zeros(len(chunk_indices))
    g = q_joint @ (y_joint - mu_joint)
    for k, idx in enumerate(chunk_indices):
        logL_chunks[k] = multivariate_normal.logpdf(x=y_joint[idx], mean=y_joint[idx] - sigma_cond_list[k] @ g[idx], cov=sigma_cond_list[k])
    return s, logL_chunks

if __name__ == '__main__':
    # Initialize the heavy model ONCE here
    model = get_model(MODEL_PATH)
    
    # Load samples
    df = read_chains("/scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/lcdm-highacc/polychord_raw/").posterior_points()  # Update path
    params = ['logA', 'ns', 'H0', 'ombh2', 'omch2', 'tau']
    
    tasks = [(s, df.iloc[s].to_dict()) for s in range(len(df))]
    
    # Run Pool
    with mp.Pool(processes=8, initializer=init_worker) as pool:
        results = list(tqdm(pool.imap(worker_eval, tasks), total=len(tasks)))
    
    # Aggregate and save results
    # 1. Initialize empty matrix: (chains=1, draws=S, observations=K)
    S = len(df)
    K = len(ingredients['chunk_indices'])
    final_logL = np.zeros((1, S, K))

    # 2. Sort and fill the results
    # (imap does not guarantee order, so we use the index 's' to place them correctly)
    for s, logL_K in results:
        final_logL[0, s, :] = logL_K

    # 3. Create ArviZ InferenceData object
    print("\nAggregation complete. Converting to NetCDF...")
    idata = az.from_dict(
        log_likelihood={"lgo_chunks": final_logL}
    )

    # 4. Save to disk
    out_filename = "lgocv_idata.nc"
    idata.to_netcdf(out_filename)
    print(f"Successfully saved to {out_filename}")
