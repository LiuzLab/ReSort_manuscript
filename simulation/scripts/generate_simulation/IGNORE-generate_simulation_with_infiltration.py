import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from matplotlib import pyplot as plt
import seaborn as sns
from tqdm import trange
from joblib import load, dump
from time import time
import os

def sample_by_cell_count(meta, ref0, frac=0.2):
    ref = ref0.raw.to_adata()
    np.random.seed(2021)
    cts = meta.columns.tolist()
    spots = meta.index.tolist()
    genes = ref.var.index.tolist()
    spot_reads_dict = {}
    
    for i in trange(len(spots)):
        spot = spots[i]
        spot_reads = []
        for ct in cts:
            # sample cells from the same cell type
            nct = meta.loc[spot,ct]
            cells_ct = ref.obs.index[ref.obs.CT == ct].tolist()
            cells = np.random.choice(cells_ct, nct, replace=False)
            # sample reads from each cell
            for cell in cells:
                vals = np.ravel(ref[cell, :].X).astype("int64")
                reads = np.repeat(genes, vals) # make gene counts as reads
                sampled_reads = list(np.random.choice(reads, int(frac * len(reads)), replace=False)) # sample reads
                spot_reads += sampled_reads
        spot_reads_dict[spot] = spot_reads
    
    # Reshape the reads data to count data
    sampled_data = np.zeros((len(spots), len(genes)))
    for i in trange(len(spots)):
        spot = spots[i]
        for j in range(len(genes)):
            gene = genes[j]
            sampled_data[i,j] = spot_reads_dict[spot].count(gene)
    sampled_data = ad.AnnData(X=sampled_data, var=ref.var, obs=meta)
    return sampled_data

def add_infiltraion(mode='cancer_inf', prop=0.01, width = 15, length=10):
    np.random.seed(2021)
    simulated_pure_regions = pd.read_csv("../simulation/simulated_spot_counts_ground_truth.csv", index_col=0)
    simulated_with_infilt = simulated_pure_regions.copy()
    
    if mode == 'cancer_inf':
        pure_cancer_spots = simulated_pure_regions.loc[(simulated_pure_regions.Ductal == 0) &  (simulated_pure_regions.Stroma==0)].index.to_numpy()
        n = len(pure_cancer_spots)
        infiltration_spots = np.random.choice(pure_cancer_spots,int(prop * n), replace=False)

        ps = np.random.normal(loc=0.5, scale=0.2, size=len(infiltration_spots))
        ps[ps > 1] = 1
        ps[ps < 0] = 0

        for i, spot in enumerate(infiltration_spots):
            cancer_cell_count = simulated_pure_regions.loc[spot, 'Cancer']
            inf_num = int(ps[i] * cancer_cell_count)
            simulated_with_infilt.loc[spot, 'Stroma'] = inf_num
            simulated_with_infilt.loc[spot, 'Cancer'] = cancer_cell_count - inf_num

    elif mode == 'random':
        spots = simulated_with_infilt.index.to_numpy()
        regions = ['Cancer', 'Stroma', 'Ductal']
        for i, region in enumerate(regions):
            np.random.seed(2021-i)
            n = simulated_with_infilt.loc[:, region].sum()
            infiltration = pd.Series(np.random.choice(spots,int(prop * n), replace=True)).value_counts()
            inf_dict = dict(infiltration)
            for spot in inf_dict.keys():
                inf_num = inf_dict[spot]
                simulated_with_infilt.loc[spot, region] += inf_num
    elif mode == 'cancer_subregion':
        pure_cancer_spots = simulated_pure_regions.loc[(simulated_pure_regions.Ductal == 0) &  (simulated_pure_regions.Stroma==0)].index.tolist()
        xs = [int(i.split("x")[0]) for i in pure_cancer_spots]
        ys = [int(i.split("x")[1]) for i in pure_cancer_spots]

        ps = np.random.normal(loc=0.5, scale=0.2, size = len(pure_cancer_spots))
        ps[ps > 1] = 1
        ps[ps < 0] = 0
        for x in range(int(np.median(xs) - width/2), int(np.median(xs) + width/2) + 1):
            for y in range(int(np.median(ys) - length/2), int(np.median(ys) + length/2) + 1):
                cancer_cell_count = simulated_with_infilt.loc[f'{x}x{y}', 'Cancer']
                stroma_count = int(cancer_cell_count * np.random.choice(ps, 1)[0])
                simulated_with_infilt.loc[f'{x}x{y}', 'Stroma'] = stroma_count
                simulated_with_infilt.loc[f'{x}x{y}', 'Cancer'] = cancer_cell_count - stroma_count
    return simulated_with_infilt

import sys
mode = sys.argv[1]
assert mode in ['cancer_inf', 'random', 'cancer_subregion']

t0 = time()
## ADD INFILTRATION
if mode in ['cancer_inf', 'random']:
    prop = float(sys.argv[2])
    folder = f"{mode}_{prop}"
    simulated_with_infilt = add_infiltraion(mode=mode, prop=prop)
else:
    if len(sys.argv) > 3:
        width = int(sys.argv[2])
        length = int(sys.argv[3])
        folder = f"{mode}_len_{length}_width_{width}"
        simulated_with_infilt = add_infiltraion(mode=mode, width=width, length=length)
    else:
        folder = f"{mode}_len_10_width_15"
        simulated_with_infilt = add_infiltraion(mode=mode)
        
if not os.path.exists(folder):
    os.mkdir(folder)
    
simulated_with_infilt.to_csv(f'{folder}/simulated_cell_counts_ground_truth.csv')
simulated_with_infilt_props = simulated_with_infilt.divide(simulated_with_infilt.sum(axis=1), axis=0)
simulated_with_infilt_props.to_csv(f"{folder}/simulated_cell_props_ground_truth.csv")

## PLOT PATTERN
pat = np.zeros((100,100,3))
inds = simulated_with_infilt_props.index.tolist()
xs = [int(i.split("x")[0]) for i in inds]
ys = [int(i.split("x")[1]) for i in inds]
for i in range(len(inds)):
    pat[xs[i], ys[i], :] = simulated_with_infilt_props.values[i]
plt.imshow(pat)
plt.axis("off")
plt.savefig(f"{folder}/simulated_cell_props_pattern_ground_truth.png", dpi=300, bbox_inches='tight')

## SAMPLE READS
scRefA = sc.read_h5ad('../simulation/PDACA_single_cell_filtered_meta_cell_types.h5')
simulated_mixture = sample_by_cell_count(simulated_with_infilt, scRefA)
dump(simulated_mixture, f'{folder}/simulated_mixture.job')
t1 = time()
print(f"Pipeline finished in {(t1-t0):.2f} seconds.")
