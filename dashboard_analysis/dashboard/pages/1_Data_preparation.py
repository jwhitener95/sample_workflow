
import streamlit as st
from PIL import Image
import os

st.set_page_config(layout="wide")

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
    plots_to_display = [(f'violin_QC_{ds_id}.png', 'Computed QC metrics', 'Violin plots are used to display the computed QC metrics -- \
                         the number of genes expressed in the count matrix, \
                         the total counts per cell, the percentage of counts in mitochondrial genes.'),
                        (f'scatter_jointQC_{ds_id}.png', 'Consider QC metrics jointly', 'A scatter plot shows the correlation between the QC \
                         metrics observed in the violin plots -- the number of genes expressed in the count matrix (y-axis), \
                         the total counts per cell (x-axis), the percentage of counts in mitochondrial genes (color).'),
                        (f'pca_variance_ratio_{ds_id}.png', 'Contribution of single PCs to the total variance in the data', 'This gives us a sense of \
                         how many PCs to use when clustering.'),
                        (f'umap_sample_{ds_id}.png', 'UMAP according to sample', 'UMAP created by embedding the neighborhood graph of cells \
                         into two dimensions for the purpose of visualization. This accounts for the PCA representation of the data matrix. \
                         The samples, which were taken as batches, should theoretically overlap. \
                         The lack of overlap can be attributed to batch effects in the clustering algorithm.')]
    for plot, cap, note in plots_to_display:
        st.write('**'+cap+'**')
        st.write(note)
        image = Image.open(f'process_datasets/figures/{plot}')
        st.image(image)

with tabs[1]:
    # Header
    st.subheader(f'Labelling output: {ds_id}')
    st.write('---')

    # Display processing plots
    res = 'leiden_res_0.50'
    plots_to_display = [(f'umap_leiden_{ds_id}.png', 'Leiden clusters at different resolutions', 'By looking at different resolutions, \
                         one can determine which resolution would be most helpful when assigning cell type labels. For the purpose of this \
                         exercise, "leiden_res_0.50" was used. Due to the small amount of samples (because of current computational power), \
                         the granularity and breadth of the resulting cell type labels may be less than satisfactory. However, the workflow \
                         shows how the analysis could be done given the proper resources. As a note, cloud services such as AWS could be used \
                         to expand the workflow if needed.'),
                        (f'dotplot_{res}_{ds_id}.png', 'Leiden clusters plotted with marker genes', 'This and the following plot were used \
                         to manually label the cells. Notably, only cell types that were fairly clear were marked in order to demonstrate the \
                         workflow. In this plot, pre-determined marker genes were highlighted to bring out the correlation between the leiden \
                         clusters and the cell types.'),
                        (f'dotplot_{res}_ranked_{ds_id}.png', 'Leiden clusters plotted with top differentially expressed genes', 'This is an alternate \
                         way to gain insight into the cell types. This requires familiarity with the genes that are most differentially expressed.'),
                        (f'umap_celltypes_{ds_id}.png', 'UMAP of manual labelling results', 'The UMAP helps to visualize the resulting cluster \
                         annotations in two dimensions.')]
    for plot, cap, note in plots_to_display:
        st.write('**'+cap+'**')
        st.write(note)
        image = Image.open(f'process_datasets/figures/{plot}')
        st.image(image)

with tabs[2]:
    # Header
    st.subheader(f'Integration results: {ds_ids}')
    st.write('---')
    st.write('This highlights the result from the integration of the two datasets. Only three random patients from each dataset were used \
             due to limited resources. Additionally, only a few cell types were labelled. This was due to the lack of confidence in assigning \
             cell types when so few cells were present. Interestingly, the marker genes do, for the most part, align with the correlated cell types. \
             The selection boxes help to distinguish the marker genes that were used. First, one selects the cell type of interest. Then, one selects \
             from the available markers for that cell type. The result can then be viewed.')

    cols = st.columns(4)
    plot_names = [i for i in os.listdir('integrate_datasets/figures/') if len(i.split('_')) > 2]
    celltypes = [i.split('_')[2] for i in plot_names]
    with cols[0]:
        cellS = st.selectbox('Choose cell type to narrow down marker selection', sorted(list(set(celltypes))))
    with cols[1]:
        genes = [i.split('_')[3].split('.')[0] for i in plot_names if i.split('_')[2] == cellS]
        geneS = st.selectbox('Choose gene to view expression', sorted(genes))

    # Display processing plots
    plots_to_display = [(f'umap_integrated.png', 'UMAP of integration results')]
    plots_to_display += [(f'umap_integrated_{cellS}_{geneS}.png', f'Single-cell expression of {geneS}')]
    col = st.columns(2)
    with col[0]:
        plot, cap = plots_to_display[0]
        st.write('**'+cap+'**')
        image = Image.open(f'integrate_datasets/figures/{plot}')
        st.image(image)
    with col[1]:
        plot, cap = plots_to_display[1]
        st.write('**'+cap+'**')
        image = Image.open(f'integrate_datasets/figures/{plot}')
        st.image(image)
