# HSC single cell analysis code

These notebooks are part of the revisions of the resubmitted publication in August 2021. 

Different parts of the analysis are:

**Notebooks 01: Preprocessing - further tests for batch correction:**  

01a: Preprocessing of all samples (treatments) together without batch correction + clustering.  

01b: Preprocessig of all samples (treatments) separately + clustering. 

01c: Overlap of gene enriched in clusters to determine which clusters are similar. Turns out GCSF, ct and indo have similar clusters (by enriched genes), therefore cells should largely overlap. DmPGE2 and poly(I:C) signatures are actually different.  

01d: Tried harmony batch correction method as an alternative to Combat. 

01e: Tried scanorama batch correction method as an alternative to Combat. 

**Notebook 02: Transcriptional scores**

Calculation of transcriptional scores in HSCs and LSKs. Calculation of Entropy and Gini score per cell (wasn't included in the final manuscript). Evaluation of mean scores in other clusters. Also random UMAP plots to update in the main manuscript. 

**Notebooks 03: Overlap of enriched genes in clusters and treatments**  

03a: Overlap of genes enriched in HSC and LSK clusters. Pathway enrichment of these enriched genes. 

03b: Overlap of genes enriched in clusters and by treatments for both HSCs and LSKs. Pathway enrichment of genes induced by treatments.

03c: Various heatmaps with genes identified in Notebook 03b.

**Notebook 04: Replotting of scATAC-Seq data**

**Notebook 05: Reexport to cellbrowser**

With updated cluster names for LSKs

**Notebook 06: Export source table for Figure 4**


**helper_functions.py**

added a number of functions