# Using ReSort to analyze epithelial-mesenchymal transition in mouse breast cancer

This folder contains subfolders:

## Raw data:
    * ./data/ (Too big, put on zenodo: doi:)

## Detect major regions for primary cell types:

    * code: 1-EMT-MIST-regions.ipynb
    * results: ./data/MIST (Too big, put on zenodo: doi:)
        
## Fig4a-e&g, Fig5:

    * 2-ReSort-EMS-proportions.ipynb (Fig4a-e&g, Fig5)

## Figure 4f, h & i:

    * ./TCGA_TNBC_validation/validate_tcga_tnbc_Fig4fhi.R
    
## ReSort primary cell type deconvolution

    * code: ./run_RCTD_primary/run_RCTD.R
    * results: ./data/ReSort_primary/primary_props.csv

## CIBERSORTx (secondary cell type)

    * results: ./data/CIBERSORTx_LM22/CIBERSORTx_LM22_result.csv
    
## ReSort secondary cell type deconvolution

    * code: 2-ReSort-EMS-proportions.ipynb
    * results: ./data/ReSort_secondary/ReSort_props.csv
