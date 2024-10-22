# Sample scRNA-seq workflow

#### ⚠️ *The current repository is in construction. Pardon the incomplete and unpolished nature of the files.* ⚠️

---

### The following datasets will be used to demonstrate the workflow

De Zuani, M., Xue, H., Park, J.S. et al. Single-cell and spatial transcriptomics analysis of non-small cell lung cancer. Nat Commun 15, 4388 (2024). [https://doi.org/10.1038/s41467-024-48700-8](https://doi.org/10.1038/s41467-024-48700-8)

Miller YE, Ghosh M, Merrick DT, Kubala B et al. Phase Ib trial of inhaled iloprost for the prevention of lung cancer with predictive and response biomarker assessment. Front Oncol 2023;13:1204726. PMID: [37711198](https://pubmed.ncbi.nlm.nih.gov/37711198/)

---

### The Streamlit dashboard can be viewed [here](https://sampleworkflow-mdgj7vhdjc5vzbj4qpjw43.streamlit.app/).

---

### collect &rarr; format &rarr; process &rarr; integrate &rarr; analyze and dashboard

---

### Collect:
Two non-small cell lung cancer (NSCLC) datasets were chosen. From these, three random patients were selected from each study due to limited computational resources (could be expanded with cloud computing if needed). Datasets with individual files per patient were chosen to lower the computational burden in processing. This resulted in features.tsv, barcodes.tsv, and matrix.mtx files for each sample (the standard 10xGenomics output). Composite CSV files and the like were not used because of the RAM needed to load the files in the initial step.

### Format:
The formatting files were used to create an AnnData file from the TSV and MTX inputs. This step would also be the step in which manual metadata could be added (a cell in each of the formatting notebooks is designated for this). The resulting AnnData object is then fed into the processing step.

### Process:
Processing primarily involves manually labelling the cell types. Prior to labelling the cells, the AnnData object is filtered and goes through a QC step. Afterwards, leiden clusters are used in conjunction with pre-determined [cell markers](https://www.sc-best-practices.org/cellular_structure/annotation.html) to label the cell types. 

### Integrate:
After the two datasets were processed and had manually labelled cell types, they were integrated using [scANVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html). The results can be viewed on the [Streamlit dashboard](https://sampleworkflow-mdgj7vhdjc5vzbj4qpjw43.streamlit.app/). The benefits of scANVI and the motivation for its implementation are its comprehensive features and scalability. scANVI builds off of scVI and additionaly allows for prediction of unlabelled cells which is particularly useful in single-cell atlas generation. A future step may involve building this into an atlas that can be used to label novel datasets via transfer learning. To accomplish this, scANVI can be used with [scArches](https://docs.scarches.org/en/latest/scanvi_surgery_pipeline.html). In this approach, the cell type information gathered from the integrated atlas are "transferred" onto the new dataset, labelling the cells. 

#### scANVI deep learning model: [source](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html)
![image](https://github.com/user-attachments/assets/7fb42a8f-a139-4e22-8060-31c102ae3fee)

#### scArches atlas generation model: [source](https://docs.scarches.org/en/latest/about.html)
![image](https://github.com/user-attachments/assets/eb800500-d90c-4553-8b34-2e2801f29f4f)

### Analyze and dashboard: 
To be elaborated
