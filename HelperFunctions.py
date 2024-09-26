
### Helper functions

import scanpy as sc
import gc
import scipy
import scrublet as scr
import numpy as np

def qc_and_filter(adata, doublet_threshold=0.3, gene_threshold=200, cell_threshold=10):
    """
    Thresholds genes, cells, and doublets. Calculates QC metrics.
    
    Args:
        adata (anndata.AnnData): Input AnnData
        sharepoint_id (str): Dataset ID
        doublet_method (str): Doublet detection method ('scrublet' or 'solo').
        doublet_threshold (float): Doublet detection threshold
    """

    # General filters
    sc.pp.filter_cells(adata, min_genes=gene_threshold)
    gc.collect()
    sc.pp.filter_genes(adata, min_cells=cell_threshold)
    gc.collect()

    # Get percent mt
    adata.var['mito'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)
    gc.collect()
    adata.obs = adata.obs.astype({'n_genes_by_counts': 'int32','total_counts': 'int32'})
    adata.obs = adata.obs.rename(columns={'pct_counts_mt': 'percent.mt', 'n_genes_by_counts': 'nFeature_RNA', 'total_counts': 'nCount_RNA'})
    adata.var = adata.var.drop('mito', axis=1)
    try:
        adata=adata[adata.obs["percent.mt"] < 20,].copy()
    except KeyError:
        print('Error with "percent.mt"')
    gc.collect()

    # Set up layers
    adata.layers['counts'] = scipy.sparse.csr_matrix(adata.X.copy())
    gc.collect()
    adata.layers['log2norm'] = sc.pp.log1p(adata.layers['counts'].copy(), base=2)
    gc.collect()
    
    # Scrublet
    scrub = scr.Scrublet(adata.X)
    doublet_scores = scrub.scrub_doublets()[0]
    adata.obs['DoubletScores'] = list(doublet_scores)
    total_cells = len(adata.obs.index)
    num_doublets = len([i for i in list(doublet_scores) if i >= doublet_threshold])
    adata = adata[adata.obs['DoubletScores'] < doublet_threshold,].copy()
    print(f'Percent doublets: {num_doublets/float(total_cells)}')

    # Plot
    sc.pl.violin(adata, ['nFeature_RNA', 'nCount_RNA', 'percent.mt'], jitter=0.4, multi_panel=True, show=False)

    return(adata)
