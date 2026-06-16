from anesthetic import read_chains
import time

start = time.time()
print("Starting load...")
df = read_chains("/scratch/tanveerk/DESI-DR2_SN-Pantheon+_CMB-CamSpecLensing/lcdm-highacc/polychord_raw/").posterior_points()
print(f"Load finished in {time.time() - start:.2f} seconds. Rows: {len(df)}")
