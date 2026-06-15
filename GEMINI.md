# AI Assistant Context: Cosmo-ELPD

## Developer Profile & Rules
* **User:** Postdoctoral researcher in observational cosmology.
* **Context:** Operating in a high-performance computing environment (Slurm) dealing with heavy MCMC/Nested Sampling chains.
* **Rules for AI:**
  1. Do not hallucinate physics parameters or simulation names. 
  2. Provide production-ready, modular Python code (prefer `src/` module imports).
  3. Keep explanations rigorous but concise. I am familiar with Bayesian statistics and large-scale structure.

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