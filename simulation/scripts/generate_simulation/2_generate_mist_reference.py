import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from matplotlib import pyplot as plt
import seaborn as sns
from tqdm import trange
from time import time
from joblib import load
from scipy.stats import ttest_ind
import sys
sys.path.append("../MIST_v0/")
import region_detection_utils
from Data import Data
import neighbors
from sklearn.metrics import adjusted_rand_score, rand_score

def region_assign(data, epi, min_region=5):
    ccs = neighbors.spatialCCs(data.nodes, data.cormat, epi)
    ccs = [[c.name for c in cc] for cc in ccs]
    
    nb_df = assign_membership(ccs, min_region)
    region_df = nb_df.loc[nb_df['region_size'] > min_region,:]
    
    regions = [cc for cc in ccs if len(cc) == min_region]
    
    for cc1 in ccs:
        if len(cc1) > min_region:
            for cc2 in ccs:
                if len(cc2) > min_region:
                    test_res =  test_identity(cc1, cc2, data.cormat, epi)
                    if test_res:
                        cc1_ind = region_df.loc[cc1, "cluster_ind"].tolist()[0]
                        cc2_ind = region_df.loc[cc2, "cluster_ind"].tolist()[0]
                        region_df.loc[region_df.cluster_ind == cc2_ind, "cluster_ind"] = cc1_ind
                        
    region_df['spot'] = region_df.index.tolist()
    region_grps  = region_df.groupby("cluster_ind")
    
    regions = []
    for name, rg in region_grps:
        regions.append(rg.spot.tolist())
        
    region_df2 = assign_membership(regions, min_region)
    region_df2 = region_df2.loc[region_df.index,:]
    region_df2['region_ind'] = region_df['cluster_ind'].tolist()
    return region_df2, regions

def test_identity(cc1, cc2, cormat, ep, min_region=5):
    #zero out diag of upper traingle of the cormat
    s1 = np.ravel(np.tril(cormat.loc[cc1, cc1], -1))
    #get lower triangle elements of cormat
    s1 = s1[np.where(s1)]
    s2 = np.ravel(np.tril(cormat.loc[cc2, cc2], -1))
    s2 = s2[np.where(s2)]
    s12 = np.ravel(cormat.loc[cc1, cc2])
    
    mean_s1 = np.mean(s1)
    mean_s2 = np.mean(s2)
    mean_s12 = np.mean(s12)
    
    diff1 = mean_s12 - mean_s1
    if len(cc1) > min_region:
        pval1 = ttest_ind(s1, s12)[1]
    else:
        pval1 = 0
        
    diff2 = mean_s12 - mean_s2
    if len(cc2) > min_region:
        pval2 = ttest_ind(s2, s12)[1]
    else:
        pval2 = 0
    mean_pass = (np.absolute(mean_s1 - mean_s12) <= 0.05) and (np.absolute(mean_s2 - mean_s12) <= 0.05)
    p_pass = (pval1 > 0.01) or (pval2 > 0.01)
    max_pass = (np.max(s12) > ep)
    result = (mean_pass or p_pass) and max_pass
    return result

def assign_membership(ccs, min_region=5):
    dfs = []
    k = 0
    for cc in ccs:
        xs = [int(c.split("x")[0]) for c in cc]
        ys = [int(c.split("x")[1]) for c in cc]
    
        df = pd.DataFrame({'x':xs, 'y':ys}, index = cc)
        df['region_size'] = len(cc)
        if len(cc) <  min_region:
            df['cluster_ind'] = -1
        else:
            df['cluster_ind'] = k
            k += 1
        dfs.append(df)
    if len(dfs) > 0:
        dfs = pd.concat(dfs)
    else:
        dfs = pd.DataFrame(columns=['x','y', 'cluster_ind'])
    return dfs

def obj_scores(ccs, cormat, min_region=5, sigma=0.02):
    #zero out diag of upper traingle of the cormat  
    term1 = 0
    n1 = 0
    for cc1 in ccs:
        if len(cc1) >= min_region: 
            s1 = np.ravel(np.tril(cormat.loc[cc1, cc1], -1))
            #get lower triangle elements of cormat
            s1 = s1[np.where(s1)]
            n1 += len(s1)
            term1 += np.sum(s1)
#             mean_s1 = np.mean(s1)
#             term1 += len(cc1) * mean_s1
    if n1 > 0:
        term1 = term1 / n1
    
    term2 = 0
    n2 = 0
    for cc1 in ccs:
        if len(cc1) >= min_region: 
            for cc2 in ccs:
                if len(cc2) >= min_region and (cc1[0] not in cc2):
                    s12 = np.ravel(cormat.loc[cc1, cc2])
                    n2 += len(s12)
                    term2 += np.sum(s12)
#                     mean_s12 = np.mean(s12)
#                     term2 += (mean_s12 * (len(cc1) + len(cc2))/2)
#    term2 = np.sqrt(np.absolute(term2))
    if n2 > 0:
        term2 = term2 / n2
    
    term3 =  sigma * len(ccs)
    scores = -1 * term1 + term2 + term3
    return scores

def select_epsilon(data, min_sim=0.4, max_sim=0.91, gap=0.02, min_region=5):
    st = time()
    eps = np.arange(min_sim, max_sim, gap)
    scores = []
    for i in trange(len(eps)):
        ep = eps[i]
        _, ccs = region_assign(data, ep,min_region)
        scores.append(obj_scores(ccs, data.cormat, min_region))
    
    ind = np.argmin(scores)
    ep, score = eps[ind], scores[ind]
    f, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(16, 7))
    ax1.plot(eps, scores)
    ax1.vlines(ep, ymin=score, ymax=np.max(scores), color='red', ls='--')
    ax1.set_ylim(score-0.1, np.max(scores) + 0.1)
    region_df, ccs = region_assign(data, ep, min_region)
    
    if len(region_df.cluster_ind.drop_duplicates().tolist()) <= 10:
        palette = 'tab10'
    else:
        palette = 'tab20'
    
    sns.scatterplot(data=region_df, x='y', y='x', hue='cluster_ind', 
                     palette=palette, legend=False, ax=ax2, size=5)
    ax2.invert_yaxis()

    ax2.set_xticks([])
    ax2.set_yticks([])
    end = time()
    print("Epsilon %.3f is selected in %.2f seconds." %(ep, end-st))
    return f, region_df

if __name__ == "__main__":
    t0 = time()
    folder = sys.argv[1]
    mixture_fn = f"{folder}/simulated_mixture.h5"
    simulated_mixture = sc.read_h5ad(mixture_fn)
    simulated_mixture.raw = simulated_mixture

    # Normalize
    sc.pp.normalize_total(simulated_mixture)
    sc.pp.log1p(simulated_mixture)
    sc.settings.set_figure_params(dpi=80, facecolor='white', figsize=(4,4))
    sc.tl.pca(simulated_mixture)
    sc.pl.pca(simulated_mixture, color=['Cancer', 'Ductal', 'Stroma'])
    t1 = time()
    print(f"Data preprocessed in {(t1-t0):.2f} seconds.")

    inds = simulated_mixture.obs.index
    xs = [int(ind.split('x')[0]) for ind in inds]
    ys = [int(ind.split('x')[1]) for ind in inds]

    mixture_meta = pd.DataFrame({'x':xs, 'y':ys}, index=inds)
    #print(simulated_mixture.var, simulated_mixture.obs, simulated_mixture.X.todense())

    mixture_df = pd.DataFrame(data=simulated_mixture.X.toarray(), 
                              index=simulated_mixture.obs.index.tolist(), 
                             columns=simulated_mixture.var.index.tolist())
    t11 = time()
    mixture_data = Data(count=mixture_df, meta=mixture_meta)
    t2 = time()
    print(f"MIST Data created in {(t2-t11):.2f} seconds.")
    f0, mixture_regions = select_epsilon(mixture_data,min_sim=0.5, 
                                         max_sim=0.99, gap=0.05 ,min_region=20)

    f0.savefig(f"{folder}/epsilon_path.png", dpi=200, bbox_inches='tight')

    t3 = time()
    print(f"Region detection finished in {(t3-t11):.2f} seconds.")

    # Plot regions
    f = plt.figure(figsize=(6,6))
    sns.set_style("dark")
    sns.scatterplot(data=mixture_regions, x="y", y="x", hue='cluster_ind', s=15, 
                    palette={0:'red', 1:'blue', 2:'green'})

    # mixture_regions.cluster_ind = mixture_regions.cluster_ind.astype(str)
    # sns.scatterplot(data=mixture_regions, x="y", y="x", hue='cluster_ind', s=15)

    plt.gca().invert_xaxis()
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('')
    plt.ylabel('')
    plt.legend(bbox_to_anchor=(1.2, 0.7))

    f.savefig(f'{folder}/MIST_detected_regions.png', bbox_inches='tight', dpi=200)
    t4 = time()
    print(f"Region plot generated in {(t4-t3):.2f} seconds.")

    simulated_meta = pd.read_csv(f"{folder}/simulated_spot_counts_ground_truth.csv", index_col=0)
    simulated_meta = simulated_meta.divide(simulated_meta.sum(axis=1), axis=0)
    simulated_meta[simulated_meta < 1] = 0
    simulated_meta.sum(axis=0)
    simulated_meta['region'] = 'Others'
    simulated_meta.loc[simulated_meta.Cancer == 1, "region"] = 'Cancer'
    simulated_meta.loc[simulated_meta.Ductal == 1, "region"] = 'Ductal'
    simulated_meta.loc[simulated_meta.Stroma == 1, "region"] = 'Stroma'

    mist_regions = []
    for ind in simulated_meta.index.tolist():
        if ind not in mixture_regions.index:
            mist_regions.append(-1)
        else:
            mist_regions.append(mixture_regions.loc[ind, 'cluster_ind'])

    ari = adjusted_rand_score(simulated_meta.region.to_numpy(), mist_regions)
    ri = rand_score(simulated_meta.region.to_numpy(), mist_regions)
    print(f"ARI = {ari}, RI = {ri}")
    mixture_raw = simulated_mixture.raw.to_adata()
    sc.pp.calculate_qc_metrics(mixture_raw, inplace=True)

    #### Plot UMIs of mixture ###
    f2 = plt.figure(figsize=(4,4))
    plt.hist(mixture_raw.obs.total_counts, bins=30)
    f2.savefig(f"{folder}/mixture_counts_hist.png", dpi=200, bbox_inches='tight')

    ### SAVE mixture counts ###
    mixture_raw_df = pd.DataFrame(data=mixture_raw.X.toarray().astype(int), index=mixture_raw.obs.index, columns=mixture_raw.var.index)
    mixture_raw_df.columns.name =None
    mixture_raw_df.to_csv(f"{folder}/simulated_mixture_raw_counts.csv")
    mixture_raw_df.to_csv(f"{folder}/simulated_mixture_raw_counts.tsv", sep='\t')
    print(f"Mixture count matrix saved.")

    ### SAVE mixture coordinates ###
    xs, ys = [i.split("x")[0] for i in mixture_raw_df.index], [i.split("x")[1] for i in mixture_raw_df.index]
    mixture_mixture_meta = pd.DataFrame({'x':xs, 'y':ys}, index=mixture_raw_df.index)
    mixture_mixture_meta.to_csv(f"{folder}/simulated_mixture_coordinates.csv")
    print(f"Mixture coordinate matrix saved.")

    ### SAVE MIST detected region reference ###
    mixture_regions['region'] = mixture_regions.cluster_ind.map({0:'Cancer', 1:'Stroma', 2:'Ductal'})
    mist_ref_meta = mixture_regions[['region']]

    mist_ref_meta.columns=['bio_celltype']
    mist_ref_meta.to_csv(f"{folder}/MIST_reference_meta.csv")
    mist_ref_meta.to_csv(f'{folder}/MIST_reference_meta.tsv', sep="\t")

    mist_ref_count = mixture_raw_df.loc[mist_ref_meta.index,:]
    mist_ref_count.to_csv(f"{folder}/MIST_reference_count.csv")
    mist_ref_count.to_csv(f"{folder}/MIST_reference_count.tsv", sep="\t")
    print(f"Regions reference saved.")
