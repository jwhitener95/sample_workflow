
### Helper functions

import scanpy as sc
import gc
import scipy
import scrublet as scr
import numpy as np

def qc_and_filter(adata, doublet_threshold=0.3, gene_threshold=200, cell_threshold=10):
    """
    Thresholds genes, cells, and doublets. Calculates QC metrics.
    """

    # General filters
    sc.pp.filter_cells(adata, min_genes=gene_threshold)
    gc.collect()
    sc.pp.filter_genes(adata, min_cells=cell_threshold)
    gc.collect()

    # Get percent mt
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    gc.collect()
    adata.obs = adata.obs.astype({'n_genes_by_counts': 'int32','total_counts': 'int32'})
    adata.obs = adata.obs.rename(columns={'pct_counts_mt': 'percent.mt', 'n_genes_by_counts': 'nFeature_RNA', 'total_counts': 'nCount_RNA'})
    adata.var = adata.var.drop('mt', axis=1)
    adata=adata[adata.obs["percent.mt"] < 20,].copy()
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

def label_cells(adata):
    """
    Label cells manually in lieu of available atlas
    """

    # Marker genes from scanpy documentation
    marker_genes = {
        "CD14+ Mono": ["FCN1", "CD14"],
        "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
        # Note: DMXL2 should be negative
        "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
        "Erythroblast": ["MKI67", "HBA1", "HBB"],
        # Note HBM and GYPA are negative markers
        "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
        "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
        "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
        "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
        # Note IGHD and IGHM are negative markers
        "B cells": [
            "MS4A1",
            "ITGB1",
            "COL4A4",
            "PRDM1",
            "IRF4",
            "PAX5",
            "BCL11A",
            "BLK",
            "IGHD",
            "IGHM",
        ],
        "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
        # Note PAX5 is a negative marker
        "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
        "CD4+ T": ["CD4", "IL7R", "TRBC2"],
        "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
        "T naive": ["LEF1", "CCR7", "TCF7"],
        "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
    }

    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    # Get most variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    # PCA
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    # Nearest neighbors and UMAP
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color="sample",
        # Setting a smaller point size to get prevent overlap
        size=2,
    )
    # # Get leiden clusters
    # sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    # for res in [0.02, 0.5, 2.0]:
    #     sc.tl.leiden(
    #         adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    #     )

    # sc.pl.umap(
    #     adata,
    #     color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    #     legend_loc="on data",
    # )

    # sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

    return adata