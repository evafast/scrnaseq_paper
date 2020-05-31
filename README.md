# HSC single cell analysis code

This repo contains jupyter notebooks containing all the analysis for the following manuscript XXX. All analysis was run in Docker containers that are available here XXX. LT and HSC is used interchangably to denote long-term (LT) HSCs - this still needs to be fixed. Poly (I:C) is denoted by pIC.

Different parts of the analysis are:

**Notebooks 00: Preprocessing steps required to initate analysis.**

00a_preprocess_CITEseq_demuxEM: This takes a raw Hashtag-derived oligo (HTO) count by cell matrix and assigns MPP type identies to it. MPP identies are used to demultiplex pooled MPP samples and repool with correct MPP numbers according to FACS

00b_preprocess_MPP_LT_overlay_correct_cell_proportions: Assigns cells in MPP samples with MPP subtype (MPP, MPP1, MPP2, MPP3/4) identity and then repools samples according to FACS frequencies

00c_preprocess_MPP_split_treatments: splits MPP by MPP subtype (MPP, MPP1, MPP2, MPP3/4) for each treatment and does filtering and processing. This is needed for MAST differential expression analysis

**Notebooks 01: Filtering, QC, batch correction and clustering of LT and MPP samples.**

Silhouette score and Davis-Boulding index are used to infer optimal clustering hyperparameters


**Notebooks 02: Validation of Clustering methodology, comparison of biological replicates**

02a and 02b: Preprocessing two independent biological replicates of control HSCs and infer optimal cluster hyperparmeters by Silhouette score and Davis-Boulding index. 

02c: Comparison of two independent biological replicates. Overlay of top enriched genes in each cluster. Comparison of cluster proportion similarity by Differential proportion analysis (DPA)

**Notebooks 03: Plotting of clusters in MPPs and LTs (Figure 1 and 3)**

03a LT: Cluster identiy assignment, gene enrichment analysis, stacked barcharts of proportions in different clusters, plotting of graphs for Figure 1

03a MPP: Cluster identiy assignment, gene enrichment analysis, stacked barcharts of proportions in different clusters, plotting of graphs for Figure 3

03b: overlap of top genes in MPP and LT clusters, DPA analysis between MPP clusters and LT clusters

03c: DPA analysis of various subdivisions in MPPs (by cluster, treatment, surface phenotype)

03d: export LT and MPP objects to cell browser

**Notebooks 04: Differential gene expression analysis using MAST**

04a: Differential expression analysis of LT using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control overall and within each cluster, controlling for number of counts and sex. 

04b: Differential expression analysis of MPP using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control overall and within each cluster, controlling for number of counts and sex. 

04c: Differential expression analysis of surfac MPPs using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control within each surface phenotype (MPP, MPP1, MPP2, MPP3/4), controlling for number of counts and sex. 


**Notebooks 05: Analysis and plotting (Figure 2 and 3) of differential gene expression**

05a: extract and overlay differential expression analysis within each treatment and for each cluster and celltype

05b: For each treatment compare gene expression between clusters and identify uniquely expressed genes

05c: Hierarchical clustering of average gene expression for differentially expressed genes for each treatment in LT

05d: Hierarchical clustering of average gene expression for differentially expressed genes for each treatment in MPP

05e: Overlap analysis of differentially expressed genes in LT and MPP. Export of average gene expression .xls

05f: Diffusion pseudotime calculation for indomethacin vs control (including GCSF)

**Notebooks 06: Sexual dimorphism in HSCs**

06a: run MAST comparing male vs female on two biological replicates of HSCs controling for cluster and number of genes

06b: split MPP and LT objects by sex

06c: Differential expression analysis of LT and MPP split by sex using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control controlling for number of counts and cluster. 

06d: Extract MAST results and analyze overlap of gene sets. Plot correlations of differential expression coefficient for male and female cells

06e: Plot results from differential gene expression analysis. Plot cell within cluster distribution by sex. Plot heatmaps for indomethacin sexual dimorphic genes

06f: test differences in cell proportions by DPA

**Notebooks 07: scATAC analysis in HSCs and MPPs**

07a: preprocess old control LT scRNAseq dataset with old scanpy version - this is for scATAC - scRNAseq integration. Needs to be processed with old scanpy version otherwise integration fails

07b: Assigns cells in MPP samples with MPP subtype (MPP, MPP1, MPP2, MPP3/4) identity and then repools samples according to FACS frequencies

07c: preprocess old control MPP scRNAseq dataset with old scanpy version - this is for scATAC - scRNAseq integration. Needs to be processed with old scanpy version otherwise integration fails

07d: preprocess and cluster LT scATAC data, run ChromVar and plot

07e: testload LT scRNAseq data into Seurat

07f: integrate scRNAseq and scATAC data for LT-HSCs (unsucessful)

07d: preprocess and cluster MPP scATAC data, run ChromVar and plot

07e: testload MPP scRNAseq data into Seurat

07f: integrate scRNAseq and scATAC data for MPPs (unsucessful)