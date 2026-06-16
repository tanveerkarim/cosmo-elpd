import os
import numpy as np
from anesthetic import read_chains as anesthetic_read
import xarray as xr
from scipy.special import logsumexp
import arviz as az 



def read_chains(path, SNIaNames, cosmoNames, type='mcmc', is_compress=True, NSamples=3000):
    """
    Reads MCMC or nested sample chains in this analysis.
    
    Args
    ----
        path: str
            Path to the directory containing the chains.
        type: str
            Type of chains to read. Options are 'mcmc' or 'nested'. Default is 'mcmc'.
        is_compress: bool
            Whether to subsample unweighted MCMC or use equal weighted for nested. Default is True.
        NSamples: int
            Number of samples to subsample if is_compress is True. Default is 3000.
    """
    # SNIaNames = ['pantheon+', 'desy5', 'desdovekie']
    # cosmoNames = ['lcdm', 'w0wacdm']
    samples = {}
    IS_MCMC = type == 'mcmc'
    
    for SNIa in SNIaNames:
        for cModel in cosmoNames:
            model = f'{SNIa}_{cModel}'
            
            # Construct folder path for the specific model with a trailing slash
            if IS_MCMC:
                model_path = os.path.join(path, model) + "/"
            else:
                model_path = os.path.join(path, model, "polychord_raw") + "/"
            
            # Check if directory exists
            if os.path.exists(model_path):
                try:
                    if IS_MCMC:
                        samples[model] = anesthetic_read(model_path).remove_burn_in(burn_in=0.3)
                    else:
                        samples[model] = anesthetic_read(model_path)
                except Exception as e:
                    print(f"Warning: Failed to load chain for {model}: {e}")
                    
    if is_compress: 
        if IS_MCMC:
            for model in list(samples.keys()): 
                try:
                    samples[model] = samples[model].compress(NSamples)
                except Exception as e:
                    print(f"Warning: Failed to compress MCMC samples for {model}: {e}")
        else: 
            for model in list(samples.keys()): 
                try:
                    samples[model] = samples[model].posterior_points()
                except Exception as e:
                    print(f"Warning: Failed to get posterior points for nested samples of {model}: {e}")
                    
    return samples

def setup_idata(samples, logL_matrix, is_w0wa=False, is_CMB=False):
    """
    Sets up the ArviZ InferenceData object for each model based on the provided log-likelihood matrix and the samples from the chains.
    
    Args:
    samples: anesthetic.Samples 
        Anesthetic Samples object containing the posterior samples of shape (Nsamples, NParameters). 
    logL_matrix: np.ndarray
        A 2D array of shape (Nchunk, NSamples) containing the log-likelihood values for each chunk and each sample.
    is_w0wa: bool
        Whether the model includes w0 and wa parameters. Default is False.
    is_CMB: bool
        Whether the model includes CMB data. Default is False.
    """
    
    # Add the 'chains' dimension -> Shape becomes (1, NSamples, NChunks)
    logL_arviz_ready = logL_matrix[np.newaxis, :, :]
    logL_arviz_ready = np.swapaxes(logL_arviz_ready, 1, 2) # (1, NChunks, NSamples)


    # Define the specific parameters to track
    if is_w0wa:
        par_list = ['logA','ns','H0','omch2','ombh2','tau', 'w', 'wa']
    else:
        par_list = ['logA','ns','H0','omch2','ombh2','tau']

    if not is_CMB:
        remove = ['logA', 'ns', 'tau']
        par_list = [x for x in par_list if x not in remove]

    # Dynamically build the dictionary in one line
    posterior_dict = {
        p: samples[p].values[np.newaxis, :] 
        for p in par_list
    }

    # Build the InferenceData Object
    idata = az.from_dict({
        "posterior":posterior_dict,
        "log_likelihood": {"chunk": logL_arviz_ready}}
        )
    
    return idata 

def get_chunk_shifted_stats(idata, chunk_idx, param_names):
    """
    Computes the PSIS-reweighted posterior mean and covariance using ArviZ loo outputs.
    
    Parameters:
    - idata: InferenceData object containing the posterior.
    - chunk_idx: The integer index of the chunk being removed.
    - param_names: List of cosmological parameters to track (e.g., ['w', 'wa']).
    """

    loo_result = az.loo(idata, pointwise=True) 
    # 1. Extract the pre-smoothed log weights for this specific chunk
    # loo_result.log_weights is an xarray; we extract the values to a numpy array.
    # Assuming the first dimension corresponds to the chunks (obs dimension).
    chunk_log_weights = loo_result.log_weights.values[chunk_idx, ...]
    
    # 2. Normalize and convert to linear weights safely
    log_norm_constant = logsumexp(chunk_log_weights)
    linear_weights = np.exp(chunk_log_weights - log_norm_constant)
    
    # Flatten the weights to 1D to match the flattened samples
    linear_weights_flat = linear_weights.flatten()
    
    # 3. Extract and stack the posterior samples
    samples = np.stack([
        idata.posterior[p].values.flatten() for p in param_names
    ], axis=-1)
    # 4. Calculate the re-weighted statistics
    mean_shifted = np.average(samples, weights=linear_weights_flat, axis=0)
    cov_shifted = np.cov(samples.T, aweights=linear_weights_flat)
    
    return mean_shifted, cov_shifted