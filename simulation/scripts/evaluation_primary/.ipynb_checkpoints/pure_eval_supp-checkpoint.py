import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np
import sys
import warnings
from scipy.stats import entropy

def read_process_res(fp, method, cell_types):
    
    if method == 'stereoscope':
        sep = "\t"
    else:
        sep = ","
        
    df = pd.read_csv(fp, sep=sep, index_col=0)

    if method == 'spatialDWLS':
        df.set_index('cell_ID', drop=True, inplace=True)
        df.index.name = None
        
    df = df[cell_types]
    
    for i in range(df.shape[0]):
        if (df.iloc[i, :] == 0).all():
            df.iloc[i, :] = 0.1
            
    df = df.divide(df.sum(axis=1), axis=0)
    return df

def evaluate_performance(est, tru, cell_types):
    """calculate cell type level Mean Absolute Percentage Error"""
    est = est.loc[tru.index, cell_types]
    tru = tru.loc[:, cell_types]
    mae = np.nanmean(np.absolute(est.values - tru.values))
    rs, kls = [], []
    for i in range(est.shape[0]):
        try:
            rs.append(pearsonr(est.values[i], tru.values[i])[0])
        except:
            rs.append(np.nan)
        kl = entropy(est.values[i]+ 0.001, tru.values[i]+ 0.001)
        
        kls.append(kl)
    ks_df = pd.DataFrame({'KL': kls}, index = tru.index)
    rs_df = pd.DataFrame({'r': rs}, index = tru.index)
    
    return rs_df, ks_df

def get_results(method, cell_types, folder, truth):
    
    if method == 'stereoscope':
        fmt = "tsv"
    else:
        fmt = "csv"
        
    f = f"{folder}{method}_results"
    
    if method != 'BayesSpace':
        # internal
        df_int = read_process_res(f"{f}/estimated_proportions_ref_internal.{fmt}", method, cell_types)  
        rs_int, kl_int = evaluate_performance(df_int, truth, cell_types)
        rs_int['reference'] = 'Internal'
        kl_int['reference'] = 'Internal'

        # external
        df_ext = read_process_res(f"{f}/estimated_proportions_ref_external.{fmt}", method, cell_types)
        rs_ext, kl_ext = evaluate_performance(df_ext, truth,cell_types)
        rs_ext['reference'] = 'External'
        kl_ext['reference'] = 'External'

        # mist
        df_mist = read_process_res(f"{f}/estimated_proportions_ref_MIST.{fmt}", method, cell_types)
        rs_mist, kl_mist = evaluate_performance(df_mist, truth,cell_types)
        rs_mist['reference'] = 'ReSort'
        kl_mist['reference'] = 'ReSort'


        rs_dfs = pd.concat([rs_int, rs_ext, rs_mist])
        kl_dfs = pd.concat([kl_int, kl_ext, kl_mist])
    else:
        df = read_process_res(f"{f}/estimated_proportions_BayesSpace.{fmt}", method, cell_types)  
        rs_dfs , kl_dfs= evaluate_performance(df, truth, cell_types)
        rs_dfs['reference'] = 'BayesSpace'
        kl_dfs['reference'] = 'BayesSpace'
        
    rs_dfs['method'] = method
    kl_dfs['method'] = method
    return rs_dfs, kl_dfs