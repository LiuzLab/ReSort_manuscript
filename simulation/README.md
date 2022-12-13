# ReSort simulation results and figures reproducing

This folder contains subfolders:

## scripts:
    * generate_simulation: The folder that generates the ST simulate data, and ReSort (MIST) reference.
    * MIST_v0: the MIST code used for the simulation study.
    * run_deconv_models: code to run other deconvolution methods.
    * evaluation_primary: code to evaluate primary deconvolution results and reproduce figures Fig1, 2 in the paper.
    * evaluation_secondary: code to evaluate secondary deconvolution results and reproduce Fig3 in the paper.
        
## model results:
    * primary:
        * pure:
            * [METHOD]_results: (result using specific deconvolution method)
                * estimated_proportions_ref_[ref].csv: results using reference internal, external, or ReSort (MIST).
        * infiltrated:
            * [METHOD]_results: (result using specific deconvolution method)
                * estimated_proportions_ref_[ref].csv: results using reference internal, external, or ReSort (MIST).
        
    * secondary_results:
        * [METHOD]: (for each deconvolution method)
            * props_internal.csv: results using internal reference
            * props_external.csv: results using external reference
            * props_immune.csv: results using immune only external reference (for the second step of ReSort)
            * props_immune.csv: results using MIST for primary cell type deconvolution (for the second step of ReSort)
    
    
## source_data: 
    h5 files to generate the internal and external references, and simulate ST samples using the internal reference. Data is too big to be put here. It is deposited at Zenodo (doi:).

## reference:
    Reference data for simulation. Data is big and deposited at Zenodo (doi:).
    
    * primary: internal, external and ReSort references to run deconvolution models for primary cell types.
        * Internal
            * count: count matrix in both .csv and .txt formats due to different algorithms' needs
            * meta: meta file specify the primary cell types in both .csv and .txt formats due to different algorithms' needs
        * External
            * count: count matrix in both .csv and .txt formats due to different algorithms' needs
            * meta: meta file specify the primary cell types in both .csv and .txt formats due to different algorithms' needs
        * ReSort
            * count: count matrix in both .csv and .txt formats due to different algorithms' needs
            * meta: meta file specify the primary cell types in both .csv and .txt formats due to different algorithms' needs
    * secondary: internal, external references to run deconvolution models for secondary cell types.
        * Internal
            * count: count matrix in both .csv and .txt formats due to different algorithms' needs
            * meta: meta file specify the primary cell types in both .csv and .txt formats due to different algorithms' needs
        * External
            * count: count matrix in both .csv and .txt formats due to different algorithms' needs
            * meta: meta file specify the primary cell types in both .csv and .txt formats due to different algorithms' needs

## simulated data:
    Simulated ST data. Data is big and deposited at Zenodo (doi:).
    
    * pure
        * H5 format of simulated data with pure regions: simulated_mixture.h5
        * CSV/TSV format of simulated count matrix with pure regions: simulated_mixture_raw_counts.c(t)sv
        * CSV/TSV format of simulated spatial coordinates for each spot: simulated_mixture_mixture_coordinates.csv
        * Ground truths files: simulated_spot_props_ground_truth_primary.csv
        * Simulated spots' UMI distributions: mixture_counts_hist.png
    * infiltrated
        * Exact the same files but with some spots from cancer regions with infiltrating immune cells.
        * Ground truth file for secondary deconvolution: simulated_spot_props_ground_truth_primary.csv