"""
Generates pivot redshift related plots. 
"""

import numpy as np 
import arviz as az
import sys 
sys.path.insert(0, "../src/")
from data_io import read_chains, setup_idata, get_chunk_shifted_stats
from metrics import compute_zp_from_cov

import matplotlib.pyplot as plt

import os
os.makedirs("../figures", exist_ok=True)

### Load chains ### 
SNIaNames = ['pantheon+', 'desy5', 'desdovekie']
cosmoNames = ['lcdm', 'w0wacdm']

NBAO_chunks = 7 # Number of DESI BAO chunks 
isDESI = True 
isCamSpec = True 

if isCamSpec:
    combo = "DESI-DR2_SNIa_CMB-CamSpec" 

project_folder = "/home/tanveerk/links/projects/def-rhlozek/tanveerk/cosmo-bayesian-model"

samples = {}
for SNIa in SNIaNames:
    for cModel in cosmoNames:
        samples[f'{SNIa}_{cModel}'] = read_chains(project_folder, SNIaNames, cosmoNames, type='nested')

### Load log-likelihood matrix ###
logL_matrix_BAO = {} 
for cModel in samples:
    logL_matrix_BAO[cModel] = np.load(f"{project_folder}/data/postprocess-chains/logL/BAO_{cModel}_desi-{isDESI}+camspec-{isCamSpec}.npz")['logL_NXS']


### Setup iData for arviz ### 
idata = {}
for cModel in samples:
    if 'w0wa' in cModel:
        is_w0wa = True
    else:
        is_w0wa = False
    idata[cModel] = setup_idata(samples[cModel], logL_matrix_BAO[cModel], is_w0wa=is_w0wa, is_CMB=isCamSpec)

### Compute ELPD_i for each chunk and each model ###
ELPD_BAO = {}
ELPD_i_BAO = {}

for cModel in samples:
    #print(cModel)
    loo_result = az.loo(idata[cModel], pointwise=True)
    #print(loo_result)
    #print(loo_result.pareto_k)
    #print("----")

    ELPD_BAO[cModel] = loo_result.elpd
    ELPD_i_BAO[cModel] = loo_result.elpd_i
    
### Compute pivot redshift from covariance matrix of (w0, wa) ###
zp_dict = {}
for cModel in idata:
    if 'lcdm' in cModel:
        continue
    
    zp_dict[cModel] = np.zeros(NBAO_chunks)
    
    for i in range(NBAO_chunks):
        _, tmpcov = get_chunk_shifted_stats(idata[cModel], i, ['w', 'wa'])
        zp_dict[cModel][i] = compute_zp_from_cov(tmpcov)
        
### Plotting helpers ### 
# Order is Pantheon+, DES Y5, DES Dovekie
markers = ['o', 's', '^']
labelName = ['Pantheon+', 'DES Y5', 'DES Dovekie']
colors = ['#cb6a49', '#a46cb7', '#7aa457']
ls = ['solid', 'dashed', 'dotted']

zBAO = np.array([0.295, 0.51 , 0.706, 0.934, 1.321, 1.484, 2.33])

### Plot 1: zp versus zBAO ### 
for i, k in enumerate(zp_dict):
    plt.plot(zBAO, zp_dict[k], marker = markers[i], c=colors[i], ls='--', label = labelName[i])
plt.xlabel(r"$z_{\rm BAO}$", fontsize=15)
plt.ylabel(r"$z_p$", fontsize=15)
plt.legend()

### Plot 2: Delta ELPD versus zBAO and cumulative Delta ELPD versus zBAO ###
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 4.5),
                        layout="constrained")
plt.savefig("../figures/zp_vs_zBAO.pdf", dpi=200, bbox_inches='tight')

### --- delta ELPD versus z --- ###
for ii in range(len(SNIaNames)):
    delta_elpd_BAO = ELPD_i_BAO[f'{SNIaNames[ii]}_lcdm'] - ELPD_i_BAO[f'{SNIaNames[ii]}_w0wacdm']
    axs[0].plot(zBAO, delta_elpd_BAO, c=colors[ii], marker=markers[ii], label = labelName[ii])

    # Cumulative sum
    cum_delta_elpd = np.cumsum(delta_elpd_BAO[::-1]) # in reverse so that high-z is on the left
    axs[1].plot(zBAO[::-1], cum_delta_elpd, c=colors[ii], marker=markers[ii], label = labelName[ii])
    
axs[0].legend()
axs[0].set_xlabel(r"Redshift $z$", fontsize=15)
axs[0].set_ylabel(r"$\Delta_{\rm ELPD} = {\rm ELPD}_{\Lambda} - {\rm ELPD}_{w_0 w_a}$", fontsize=15)
axs[0].axhline(0, c = 'k', ls ='-.')
axs[0].invert_xaxis()
axs[0].text(0.37, 0.08, "BGS", fontsize=12, color='b')
axs[0].text(0.6, 0.35, "LRG1", fontsize=12, color='b')
axs[0].text(0.80, -1.1, "LRG2", fontsize=12, color='b')
axs[0].text(0.9, -0.5, "LRG3 + ELG1", fontsize=12, color='b')
axs[0].text(1.3, -0.1, "ELG2", fontsize=12, color='b')
axs[0].text(1.5, 0.1, "QSO", fontsize=12, color='b')
axs[0].text(2.4, 0.1, r"Ly-$\alpha$", fontsize=12, color='b')
axs[0].tick_params(axis='both', labelsize=14)

### --- Plot 2: cumulative delta ELPD versus z --- ###
axs[1].invert_xaxis()
axs[1].set_xlabel(r"Redshift $z$", fontsize=15)
axs[1].set_ylabel(r"Cumulative $\Delta_{\mathrm{ELPD}} = {\rm ELPD}_{\Lambda} - {\rm ELPD}_{w_0 w_a}$")
axs[1].axhline(0, linestyle='--', c ='k')
axs[1].text(0.37, -.9, "BGS", fontsize=12, color='b')
axs[1].text(0.6, -1.05, "LRG1", fontsize=12, color='b')
axs[1].text(0.80, -1.5, "LRG2", fontsize=12, color='b')
axs[1].text(0.9, -0.5, "LRG3 + ELG1", fontsize=12, color='b')
axs[1].text(1.3, 0.1, "ELG2", fontsize=12, color='b')
axs[1].text(1.55, 0.1, "QSO", fontsize=12, color='b')
axs[1].text(2.4, 0.1, r"Ly-$\alpha$", fontsize=12, color='b')
axs[1].tick_params(axis='both', labelsize=14)
plt.savefig("../figures/delta_elpd_vs_zBAO.pdf", dpi=200, bbox_inches='tight')

### Plot 3: delta ELPD versus zp ###
for ii in range(len(SNIaNames)):
    delta_elpd_BAO = ELPD_i_BAO[f'{SNIaNames[ii]}_lcdm'] - ELPD_i_BAO[f'{SNIaNames[ii]}_w0wacdm']
    plt.scatter(zp_dict[f'{SNIaNames[ii]}_w0wacdm'], delta_elpd_BAO, c=zBAO, marker=markers[ii], cmap='viridis',  label = labelName[ii])
    
plt.legend()
plt.xlabel(r"$z_p$", fontsize=15)
plt.ylabel(r"$\Delta_{\rm ELPD} = {\rm ELPD}_{\Lambda} - {\rm ELPD}_{w_0 w_a}$", fontsize=15)
plt.axhline(0, c = 'k', ls ='-.')
plt.colorbar(label=r'$z_{\rm BAO}$')
plt.savefig("../figures/delta_elpd_vs_zp.pdf", dpi=200, bbox_inches='tight')

### Plot 4: Impact of removing each chunk on the pivot redshift on w0-wa plane ### 
from matplotlib.patches import Ellipse

def draw_ellipse(mean, cov, ax, color, label, linestyle='-'):
    """Helper function to mathematically draw 1-sigma and 2-sigma error ellipses."""
    vals, vecs = np.linalg.eigh(cov)
    
    # Sort eigenvalues in descending order
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    
    # Calculate rotation angle of the ellipse
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    # Draw 1-sigma (68%) and 2-sigma (95%) contours (scale factors 1.515, 2.447)
    for n_sig, alpha in zip([1.515, 2.447], [0.8, 0.4]):
        width, height = 2 * n_sig * np.sqrt(vals)
        ellip = Ellipse(xy=mean, width=width, height=height, angle=theta, 
                        facecolor='none', edgecolor=color, linestyle=linestyle, 
                        linewidth=2, alpha=alpha)
        ax.add_patch(ellip)
    
    # Dummy line for the legend
    ax.plot([], [], color=color, linestyle=linestyle, label=label)


# ---------------------------------------------------------
# Grid Setup
# ---------------------------------------------------------
chunk_names = ['BGS', 'LRG1', 'LRG2', 'LRG3+ELG1', 'ELG2', 'QSO', 'Lya']
w0wa_models = [m for m in idata.keys() if 'lcdm' not in m]

# Setup the 7x3 figure. Sharing axes by column keeps the SN baseline fixed for visual comparison
fig, axes = plt.subplots(len(chunk_names), len(w0wa_models), 
                         figsize=(18, 28), sharex='col', sharey='col')

# ---------------------------------------------------------
# Plotting Loop
# ---------------------------------------------------------
for row_idx, chunk_name in enumerate(chunk_names):
    for col_idx, cModel in enumerate(w0wa_models):
        
        # Isolate the specific subplot
        ax = axes[row_idx, col_idx]
        
        # 1. Baseline Statistics (All Data)
        samples_base = np.stack([
            idata[cModel].posterior['w'].values.flatten(),
            idata[cModel].posterior['wa'].values.flatten()
        ], axis=-1)
        
        mean_base = np.mean(samples_base, axis=0)
        cov_base = np.cov(samples_base.T)
        
        # 2. Shifted Statistics (Iterative Chunk Removed)
        mean_shifted, cov_shifted = get_chunk_shifted_stats(idata[cModel], row_idx, ['w', 'wa'])
        
        # 3. Draw Ellipses
        draw_ellipse(mean_base, cov_base, ax, color='black', label='Baseline')
        draw_ellipse(mean_shifted, cov_shifted, ax, color='indigo', 
                     label=f'Shifted (No {chunk_name})', linestyle='--')
        
        # 4. Crosshairs
        ax.axvline(-1, color='gray', linestyle=':', alpha=0.5)
        ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
        
        # ---------------------------------------------------------
        # Formatting & Labeling
        # ---------------------------------------------------------
        ax.margins(0.2)
        
        # Only label the X-axis on the very bottom row
        if row_idx == len(chunk_names) - 1:
            ax.set_xlabel(r'$w_0$', fontsize=16)
            
        # Only label the Y-axis on the far left column
        if col_idx == 0:
            ax.set_ylabel(r'$w_a$', fontsize=16)
            
        # Only place SN Anchor titles on the very top row
        if row_idx == 0:
            title_str = cModel.split('_')[0].capitalize()
            ax.set_title(f'SN Anchor: {title_str}', fontsize=16)
            
        # Add legend to every box to identify the chunk removed
        ax.legend(loc='best', fontsize=10)

# Reduce whitespace between the shared axes
plt.tight_layout()
plt.subplots_adjust(hspace=0.05, wspace=0.05)
plt.savefig("../figures/zBAO_shift_w0wa.pdf", dpi=200, bbox_inches='tight')