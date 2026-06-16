# AI Assistant Context: Cosmo-ELPD

## Developer Profile & Rules
* **User:** Postdoctoral researcher in observational cosmology.
* **Context:** Operating in a high-performance computing environment (Slurm) dealing with heavy MCMC/Nested Sampling chains.
* **Rules for AI:**
  1. Do not hallucinate physics parameters or simulation names. 
  2. Provide production-ready, modular Python code (prefer `src/` module imports).
  3. Keep explanations rigorous but concise. I am familiar with Bayesian statistics and large-scale structure.
  4. If confused, ask for detailed clarification than assuming.
  5. Ask one question at a time to keep your context clean. 

## Data Access
- You have access to the local version of the Github to check for code, maths and physics-based bugs. The main chains live in the cluster system which you do not have access to. 
- - You do not have access to the `scratch/` directory as a result, so you cannot execute those.
- - If you need to test those, ask for a prototype version. 
- Two representative examples of nested and mcmc chains are in @chains/mcmc and @chains/nested. 
- - The MCMC chains are produced by `cobaya` and the nested chain is produced by `polychord` via `cobaya`. 
- The data is accessed using the `anesthetic` package format for ease of use. 

## Leave-Redshift-Out Cross-Validation 
- We implement chunking of BAO and SNIa data by redshift chunks. 
- Both DESI BAO and SNIa likelihoods are Gaussian, so we can use non-factorizable log likelihoods to compute ELPD even if the covariance matrices are dense. 
- Non-Factorized Protocol (Bürkner et al. 2021)
 - For correlated data (e.g., SN, CMB), implement optimized LOO-CV approximation to maintain $O(SN^k)$ complexity from Bürkner et al. (2021). 
 - - Compute conditional mean: $\tilde{\mu}_i = y_i - \frac{g_i}{\bar{\sigma}_{ii}}$ 
 - - Compute conditional variance: $\tilde{\sigma}_i^2 = \frac{1}{\bar{\sigma}_{ii}}$ 
 - - Final Point-wise Log-Likelihood: $$\log p(y_i|y_{-i}, \theta) = -\frac{1}{2}\log(2\pi\tilde{\sigma}_i^2) - \frac{1}{2}\frac{(y_i - \tilde{\mu}_i)^2}{\tilde{\sigma}_i^2}$$
- These data are stored in @data/postprocess-chains/logL 

## Current State
* **Active Branch:** `analysis/elpd-vs-zp`
* **Current Objective:** Calculate the ELPD versus pivot redshift ($z_p$) metric to observe how $z_p$ shifts when specific data chunks are removed using PSIS inferred posteriors.

## Task Checklist
- [ ] Write PSIS posterior reader in `src/data_io.py` to parse chunk-dropped weights.
- [ ] Implement the mathematical $z_p$ shifting optimization in `src/metrics.py`.
- [ ] Write the executable plotting script in `jobs/generate_zp_plots.py`.
- [ ] Validate YAML configs in `config/` match the base cosmology.

## Active Blockers / Scratchpad
* *[Leave this blank for now. When you hit a bug, paste the traceback or error message here before asking for help.]*