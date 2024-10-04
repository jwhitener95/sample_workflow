
import streamlit as st

# Header
st.subheader(f'Sample scRNA-seq workflow')
st.write('---')

st.write('⚠️ The current dashboard is in construction ⚠️')
st.write('---')
url_dezuani = 'https://doi.org/10.1038/s41467-024-48700-8'
url_Miller = 'https://pubmed.ncbi.nlm.nih.gov/37711198/'
# url_Wang = 'https://pubmed.ncbi.nlm.nih.gov/39113235/'
st.write('**The following datasets will be used to demonstrate the workflow**')
st.write('De Zuani, M., Xue, H., Park, J.S. et al. Single-cell and spatial transcriptomics analysis of non-small cell lung cancer. Nat Commun 15, 4388 (2024). [https://doi.org/10.1038/s41467-024-48700-8](%s)' % url_dezuani)
st.write('Miller YE, Ghosh M, Merrick DT, Kubala B et al. Phase Ib trial of inhaled iloprost for the prevention of lung cancer with predictive and response biomarker assessment. Front Oncol 2023;13:1204726. PMID: [37711198](%s)' % url_Miller)
# st.write('Wang D, Li S, Yang Z, Yu C et al. Single-cell transcriptome analysis deciphers the CD74-mediated immune evasion and tumour growth in lung squamous cell carcinoma with chronic obstructive pulmonary disease. Clin Transl Med 2024 Aug;14(8):e1786. PMID: [39113235](%s)' % url_Wang)

st.write('---')
st.subheader('⇦ Select a page in the sidebar to begin')
st.write('---')

git_link = 'https://github.com/jwhitener95/sample_workflow'
st.write('Visit the [GitHub page](%s) for further details' % git_link)
