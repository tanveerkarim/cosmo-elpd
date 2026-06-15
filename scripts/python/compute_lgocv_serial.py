import numpy as np
from anesthetic import read_chains
from cobaya.model import get_model

import arviz as az

import yaml 

# Load pre computed files 
precomputed = np.load("../data/lgo_cv_ingredients.npz", allow_pickle=True)
y_joint = precomputed['y_joint']
q_joint = precomputed['q_joint']
chunk_indices = precomputed['chunk_indices'][()]
sigma_cond_list = precomputed['sigma_cond_list']

#Load nested sample chains
nested_samples = {}

nested_samples['lcdm'] = read_chains("/scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/lcdm-highacc/polychord_raw/").posterior_points()
nested_samples['w0wacdm'] = read_chains("/scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/w0waCDM-db/polychord_raw/").posterior_points()

#Setup ELPD calculations
# 2. Load your MCMC chain (e.g., from Cobaya/GetDist outputs)
# Let's assume 'chain' is an array of shape (S, num_params)

log_lik_matrix = {}

K = len(chunk_indices)

for model in nested_samples:
    S = len(nested_samples[model])    

    # Initialize the exact array ArviZ wants
    # Shape is (chains, draws, observations/chunks) -> (1, S, K)
    log_lik_matrix[model] = np.zeros((1, S, K))

    print(log_lik_matrix[model].shape)

# Setup cosmology pipeline 
model = {}

# Load YAML configurations
yaml_path = "/scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-lcdm-highAc.yaml"
model['lcdm'] = get_model(yaml_path)

yaml_path = "/scratch/tanveerk/bayesian-model-workspace/DESI-DR2_Pantheon+_CMB-CamSpec-w0wacdm-db.yaml"
model['w0wacdm'] = get_model(yaml_path)

# Setup nonfactorizable logL
def get_mu_joint(theta_sample, param_names, model, like_sn, like_bao):
    """
    Extracts the uncalibrated SNIa distance moduli and the BAO theory vector.
    
    Parameters:
    -----------
    theta_sample : pd.Series
        A single row from your Pandas DataFrame (retains the MultiIndex).
    param_names : list of str
        The exact parameter strings Cobaya needs (e.g., ['logA', 'ns', ...]).
    """
    # 1. droplevel(1) removes the LaTeX formatting level on the fly.
    # 2. [param_names] slices out only the cosmological parameters we need.
    # 3. to_dict() packages it for Cobaya.
    # The original DataFrame remains completely untouched.
    theta_dict = theta_sample.droplevel(1)[param_names].to_dict()
    
    # Force Cobaya to compute the cosmology for this draw.
    model.logposterior(theta_dict)
    provider = model.provider
    
    # ---------------------------------------------------------
    # PART A: Pantheon+ Theory (Uncalibrated)
    # ---------------------------------------------------------
    da_sn = provider.get_angular_diameter_distance(like_sn.zcmb)
    mu_sn = 5.0 * np.log10((1.0 + like_sn.zhel) * (1.0 + like_sn.zcmb) * da_sn)
    
    # ---------------------------------------------------------
    # PART B: DESI DR2 BAO Theory
    # ---------------------------------------------------------
    # The BAO base class iterates through its data dataframe and evaluates 
    # the specific observable (e.g., 'DM_over_rs', 'Hz_rs') at each redshift.
    # Because `like_bao` shares the updated provider, theory_fun() works instantly.
    
    mu_bao_list = [
        like_bao.theory_fun(z, obs)
        for z, obs in zip(like_bao.data["z"], like_bao.data["observable"])
    ]
    
    # Replicate the exact array casting used in the BAO logp() source code
    # The .T[0] handles extracting the scalars if the provider returned 1D arrays
    mu_bao = np.array(mu_bao_list).T[0]
            
    # ---------------------------------------------------------
    # PART C: Concatenate
    # ---------------------------------------------------------
    mu_joint = np.concatenate([mu_sn, mu_bao])
    
    return mu_joint

# Calculate logL per block 
from scipy.stats import multivariate_normal
def cluster_chunk_loglikes(mu_joint, y_joint, q_joint, chunk_indices, sigma_cond_list):
    K = len(chunk_indices)
    logL_chunks = np.zeros(K)
    g = q_joint @ (y_joint - mu_joint) # Global residual
    
    for k, idx in enumerate(chunk_indices):
        Sigma_cond = sigma_cond_list[k]
        y_I = y_joint[idx]
        g_I = g[idx]
        
        mu_cond = y_I - Sigma_cond @ g_I
        logL_chunks[k] = multivariate_normal.logpdf(x=y_I, mean=mu_cond, cov=Sigma_cond)
        
    return logL_chunks

#Main execution
cosmo_params = ['logA', 'ns', 'omch2', 'ombh2', 'tau', 'H0', 'amp_143x217', ]
nuisance_params = ['amp_143x217', 'A_planck', 'n_217', 'n_143', 'amp_143', 'calEE', 'calTE', 'n_143x217', 'amp_217']
param_names = np.concatenate((cosmo_params, nuisance_params))

from tqdm import tqdm
for cModel in log_lik_matrix:
    print(cModel)
    # load the instances of the likelihoods
    like_sn = model[cModel].likelihood['sn.pantheonplus']
    like_bao = model[cModel].likelihood['bao.desi_dr2']

    S = log_lik_matrix[cModel].shape[1]
    K = log_lik_matrix[cModel].shape[2]
    print(f"Evaluating {S} samples across {K} chunks...")

    for s in tqdm(range(S)):
        
        # per equal weighted point
        row_series = nested_samples[cModel].iloc[s]

        # compute model mean 
        mu_joint_s = get_mu_joint(row_series, param_names, model[cModel], like_sn, like_bao)

        # compute logL
        logL_K = cluster_chunk_loglikes(mu_joint_s, y_joint, q_joint, chunk_indices, sigma_cond_list)
        
        # store value 
        log_lik_matrix[cModel][0, s, :] = logL_K
