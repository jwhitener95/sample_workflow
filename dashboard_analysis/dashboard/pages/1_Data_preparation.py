
import streamlit as st
from PIL import Image
import os

# Dataset id
ds_ids = ['DeZuani2024', 'Miller2023']
ds_id = st.sidebar.selectbox('Select dataset to view', ds_ids)
st.sidebar.caption('*Selected dataset is used for the "Processing" and "Labeling" tabs only')

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
        image = Image.open(f'process_datasets/figures/{plot}')
        st.image(image)

with tabs[2]:
    # Header
    st.subheader(f'Integration results: {ds_ids}')
    st.write('---')

    cols = st.columns(2)
    plot_names = [i for i in os.listdir('integrate_datasets/figures/') if len(i.split('_')) > 2]
    celltypes = [i.split('_')[2] for i in plot_names]
    with cols[0]:
        cellS = st.selectbox('Choose cell type to narrow down marker selection', sorted(list(set(celltypes))))
    with cols[1]:
        genes = [i.split('_')[3].split('.')[0] for i in plot_names if i.split('_')[2] == cellS]
        geneS = st.selectbox('Choose gene to view expression', genes)

    # Display processing plots
    plots_to_display = [(f'umap_integrated.png', 'UMAP of integration results')]
    plots_to_display += [(f'umap_integrated_{cellS}_{geneS}.png', f'Single-cell expression of {geneS}')]
    # col = st.columns(2)
    # with col[0]:
    plot, cap = plots_to_display[0]
    st.write('**'+cap+'**')
    image = Image.open(f'integrate_datasets/figures/{plot}')
    st.image(image)
    # with col[1]:
    plot, cap = plots_to_display[1]
    st.write('**'+cap+'**')
    image = Image.open(f'integrate_datasets/figures/{plot}')
    st.image(image)
