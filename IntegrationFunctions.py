
### Integration functions

import gc
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
from collections import Counter
import scvi

def downsample(adata):
    """
    Separate rare and common celltypes. Minimize common cell types and then concatenate.
    """

    celltype_threshold = 5000
    max_per_study = 2000
    category = "cell_type"

    counts = adata.obs[category].value_counts()
    rare_cells = [i for i in counts.index if counts[i] < celltype_threshold]
    common_cells = [i for i in list(set(adata.obs[category])) if i not in rare_cells]
    # Separate rare and common celltypes
    rare = adata[adata.obs[category].isin(rare_cells),].copy()
    common = adata[adata.obs[category].isin(common_cells),].copy()
    for j in list(pd.unique(common.obs[category])):
        minimized = common[common.obs[category] == j].copy()
        if len(minimized.obs.index) > max_per_study:
            sc.pp.subsample(minimized, n_obs=max_per_study)
    adata = anndata.concat([rare, minimized])
    return adata

def load_datasets(path):
    """
    Read files and downsample common cell types.
    """

    adata = sc.read_h5ad(path)
    adata.obs_names_make_unique()
    adata.obs_names = (path.split('/')[-1].split('_')[0] + "_" + adata.obs_names)
    adata = downsample(adata)
    gc.collect()
    return adata

def concat_datasets(datasets):
    """
    Load datasets (read and downsample) and then concatenate. 
    """

    adatas = []
    for ds in datasets:
        adata1 = load_datasets(ds)
        adata1.var_names_make_unique()
        adatas += [adata1]
    adata = anndata.concat(adatas)
    gc.collect()
    adata.obs_names_make_unique()
    return adata

def process_datasets(dataset):
    """
    1. Filter genes. 
    2. Filter cells. 
    3. Remove mitochondrial genes. 
    4. Filter batches based on threshold number of cells.
    5. Filter doublets. 
    6. Filter highly variable genes.
    """

    batch_id = "ds_sample"
    min_genes = 200
    min_cells = 20
    mitoch_max = 20
    batch_count = 50
    n_gen = 4000
    doubl_lim = 0.3

    # 1
    sc.pp.filter_genes(dataset, min_cells=min_cells)
    gc.collect()
    # 2
    sc.pp.filter_cells(dataset, min_genes=min_genes)
    gc.collect()
    # 3
    try:
        dataset = dataset[dataset.obs["percent.mt"] < mitoch_max].copy()
        gc.collect()
    except:
        pass # percent.mt not in dataset.obs
    # 4
    if ("sample" in dataset.obs.columns) and ("ds_id" in dataset.obs.columns):
        dataset.obs["ds_sample"] = (dataset.obs["ds_id"].astype(str) + "_" + dataset.obs["sample"].astype(str))
    t = Counter(dataset.obs[batch_id])
    df = pd.DataFrame.from_dict(t, orient="index").reset_index()
    df = df.rename(columns={"index": batch_id, 0: "count"})
    to_keep = [df.loc[i, batch_id] for i in df.index if df.loc[i, "count"] >= batch_count]
    dataset = dataset[dataset.obs[batch_id].isin(to_keep)].copy()
    gc.collect()
    # 5
    dataset = dataset[dataset.obs["DoubletScores"] < doubl_lim].copy()
    gc.collect()
    # 6
    sc.pp.highly_variable_genes(
        dataset,
        layer="log2norm",
        n_top_genes=n_gen,
        batch_key=batch_id,
        flavor="cell_ranger",
        inplace=True
    )
    hvg = [i for i in dataset.var.index if dataset.var.loc[i, "highly_variable"] == True]
    dataset = dataset[:, hvg].copy()
    gc.collect()
    return dataset

def integrate_datasets(adata):
    """
    1. Fill nans.
    2. Setup model and train.
    3. Compute nearest neighbrs, get leiden clusters, and plot
    """

    batch = "ds_sample"
    labels = "cell_type"
    layer = "counts"
    method = 'umap'
    hyperparameters = {"n_layers": 5,"n_hidden": 240,"n_latent": 10,"dropout_rate": 0.15}
    scvi.settings.seed = 0
    train_args = {"max_epochs": 100, "early_stopping": True, "check_val_every_n_epoch": 1, "early_stopping_patience": 4,
                  "early_stopping_min_delta": 1, "early_stopping_monitor": "reconstruction_loss_validation"}

    # 1
    nan_columns = adata.obs.isna().sum().index[adata.obs.isna().sum() > 0].tolist()
    for col in nan_columns:
        if adata.obs[col].dtype in ["float64", "float32", "int32", "int64"]:
            adata.obs.loc[:, col] = adata.obs.loc[:, col].fillna(0)
        else:
            adata.obs.loc[:, col] = adata.obs.loc[:, col].astype(str).fillna("Unknown")
    # 2
    scvi.model.SCANVI.setup_anndata(adata, layer=layer, batch_key=batch, labels_key=labels, unlabeled_category='Undetermined')
    model = scvi.model.SCANVI(adata, **hyperparameters)
    model.train(**train_args, n_samples_per_label=100)
    adata.obsm["X_scANVI"] = model.get_latent_representation()
    # 3
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=15, method=method)
    sc.tl.leiden(adata)
    sc.tl.umap(adata, method=method)
    return (adata, model)
