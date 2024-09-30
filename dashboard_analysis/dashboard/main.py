
import streamlit as st
# import scanpy as sc
# import matplotlib.pyplot as plt

tabs = st.tabs(['Processing'])

with tabs[0]:
    # Dataset id
    ds_ids = ['Miller2023']
    cols = st.columns(4)
    with cols[0]:
        ds_id = st.sidebar.selectbox('Select dataset to view', ds_ids)

    # Header
    st.subheader(f'Processing output: {ds_id}')
    st.write('---')

    # Display processing plots
    base_path = f'plots/{ds_id}/'
    plots_to_display = [('mito_plot1.png', 'Computed QC metrics'),
                        ('mito_plot2.png', 'Consider QC metrics jointly'),
                        ('variance_ratio.png', 'Contribution of single PCs to the total variance in the data'),
                        ('umap_sample.png', 'UMAP according to sample'),
                        ('leiden_plots.png', 'Leiden clusters at different resolutions')]
    for plot, cap in plots_to_display:
        st.write('**'+cap+'**')
        st.image(base_path + plot)


# # Base file path
# base = f'D:/GitHub/Data/NSCLC/{ds_id}/'

# # Path to AnnData
# adata_path = f'{base}{ds_id}_anndata_processed.h5ad'

# @st.cache_data
# def load_data():
#     adata = sc.read_h5ad(adata_path)
#     return adata
# adata = load_data()

# st.write(adata)

# fig, ax = plt.subplots()
# sc.tl.umap(adata)
# sc.pl.umap(
#     adata,
#     color="sample",
#     size=2,
# )
# cols = st.columns(3)
# with cols[0]:
#     fig, ax = plt.subplots()
#     sc.pl.umap(
#         adata,
#         color="leiden_res_0.02",
#         legend_loc="on data",
#         ax=ax
#     )
#     st.pyplot(fig)
# with cols[1]:
#     fig, ax = plt.subplots()
#     sc.pl.umap(
#         adata,
#         color="leiden_res_0.50",
#         legend_loc="on data",
#         ax=ax
#     )
#     st.pyplot(fig)
# with cols[2]:
#     fig, ax = plt.subplots()
#     sc.pl.umap(
#         adata,
#         color="leiden_res_2.00",
#         legend_loc="on data",
#         ax=ax
#     )
#     st.pyplot(fig)