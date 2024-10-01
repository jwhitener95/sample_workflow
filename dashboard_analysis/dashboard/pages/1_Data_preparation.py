
import streamlit as st
from PIL import Image

# Dataset id
ds_ids = ['Miller2023']
ds_id = st.sidebar.selectbox('Select dataset to view', ds_ids)

tabs = st.tabs(['Processing', 'Labelling', 'Integrating'])

with tabs[0]:
    # Header
    st.subheader(f'Processing output: {ds_id}')
    st.write('---')

    # Display processing plots
    plots_to_display = [(f'violin_QC_{ds_id}.png', 'Computed QC metrics'),
                        (f'scatter_jointQC_{ds_id}.png', 'Consider QC metrics jointly'),
                        (f'pca_variance_ratio_{ds_id}.png', 'Contribution of single PCs to the total variance in the data'),
                        (f'umap_sample_{ds_id}.png', 'UMAP according to sample')]
    for plot, cap in plots_to_display:
        st.write('**'+cap+'**')
        image = Image.open(f'process_datasets/figures/{plot}')
        st.image(image)

with tabs[1]:
    # Header
    st.subheader(f'Labelling output: {ds_id}')
    st.write('---')

    # Display processing plots
    res = 'leiden_res_0.50'
    plots_to_display = [(f'umap_leiden_{ds_id}.png', 'Leiden clusters at different resolutions'),
                        (f'dotplot_{res}_{ds_id}.png', 'Leiden clusters plotted with marker genes'),
                        (f'dotplot_{res}_ranked_{ds_id}.png', 'Leiden clusters plotted with top differentially expressed genes'),
                        (f'umap_celltypes_{ds_id}.png', 'UMAP of manual labelling results')]
    for plot, cap in plots_to_display:
        st.write('**'+cap+'**')
        image = Image.open(f'/workspaces/sample_workflow/process_datasets/figures/{plot}')
        st.image(image)

with tabs[2]:
    # Header
    st.subheader(f'Integration results: {ds_ids}')
    st.write('---')







    # im = cv2.imread(f'dashboard_analysis/dashboard/static/{ds_id}/{plots_to_display[0][0]}')
    # # im_resize = cv2.resize(im, (500, 500))
    # is_success, im_buf_arr = cv2.imencode(".jpg", im)
    # byte_im = im_buf_arr.tobytes()
    # st.image(byte_im, caption=['Original Image'])


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