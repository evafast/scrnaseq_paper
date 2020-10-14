# HSC single cell analysis code

This repo contains jupyter notebooks containing all the analysis for the following manuscript XXX. All analysis was run in Docker containers that are available here https://hub.docker.com/u/evafast1. LT and HSC is used interchangably to denote long-term (LT) HSCs - this still needs to be fixed. MPP nomenclature was at some point clarified to LSK = HSC + MPP, MPP is only progenitors, but there might still be inconsistent naming in parts of the notebooks. Poly (I:C) is denoted by pIC.

Different parts of the analysis are:

**Notebooks 00: Preprocessing steps required to initate analysis.**

00a: This takes a raw Hashtag-derived oligo (HTO) count by cell matrix and assigns MPP type identies to it. MPP identies are used to demultiplex pooled MPP samples and repool with correct MPP numbers according to FACS

00b: Assigns cells in MPP samples with MPP subtype (MPP, MPP1, MPP2, MPP3/4) identity and then repools samples according to FACS frequencies

00c: splits MPP by MPP subtype (MPP, MPP1, MPP2, MPP3/4) for each treatment and does filtering and processing. This is needed for MAST differential expression analysis

**Notebooks 01: Filtering, QC, batch correction and clustering of LT and MPP samples.**

Silhouette score and Davis-Boulding index are used to infer optimal clustering hyperparameters. Notebook 01c is to investigate dmPGE2 induced shift in surface marker poroportions. 

**Notebooks 02: Validation of Clustering methodology, comparison of biological replicates**

02a and 02b: Preprocessing two independent biological replicates of control HSCs and infer optimal cluster hyperparmeters by Silhouette score and Davis-Boulding index. 

02c: Comparison of two independent biological replicates. Overlay of top enriched genes in each cluster. Comparison of cluster proportion similarity by Differential proportion analysis (DPA)

**Notebooks 03: Plotting of clusters in MPPs and LTs (Figure 1 and 3)**

03a: HSC cluster identiy assignment, pathway enrichment analysis, stacked barcharts of proportions in different clusters, plotting for Figure 1

03b: LSK cluster identiy assignment, pathway enrichment analysis, stacked barcharts of proportions in different clusters, plotting of graphs for Figure 2

03c: overlap of top genes in MPP and LT clusters, DPA analysis between MPP clusters and LT clusters

03d: DPA analysis of various subdivisions in MPPs (by cluster, treatment, surface phenotype)

03e: export LT and MPP objects to cell browser

**Notebooks 04: Differential gene expression analysis using MAST**

04a: Differential expression analysis of HSCs using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control overall and within each cluster, controlling for number of counts and sex. 

04b: Differential expression analysis of LSKs using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control overall and within each cluster, controlling for number of counts and sex. 

04c: Differential expression analysis of surface MPPs using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control within each surface phenotype (MPP, MPP1, MPP2, MPP3/4), controlling for number of counts and sex. 

**Notebooks 05: Analysis and plotting (Figure 3 and 4) of differential gene expression**

05a: extract and overlay differential expression analysis within each treatment and for each cluster and celltype

05b: For each treatment compare gene expression between clusters and identify uniquely expressed genes

05c: Hierarchical clustering of average gene expression for differentially expressed genes for each treatment in HSC

05d: Hierarchical clustering of average gene expression for differentially expressed genes for each treatment in LSK

05e: Overlap analysis of differentially expressed genes in LT and LSK. Export of average gene expression .xls

05f: Diffusion pseudotime calculation for indomethacin vs control (including GCSF)

**Notebooks 06: Sexual dimorphism in HSCs**

06a: run MAST comparing male vs female on two biological replicates of HSCs controling for cluster and number of genes

06b: split LSK and HSC objects by sex

06c: Differential expression analysis of HSC and LSK split by sex using MAST. Differential expression of raw counts from each drug treatment (dmPGE2, poly (I:C), G-CSF, indomethacin) is compared vs control controlling for number of counts and cluster. 

06d: Extract MAST results and analyze overlap of gene sets. Plot correlations of differential expression coefficient for male and female cells

06e: Plot results from differential gene expression analysis. Plot cell within cluster distribution by sex. Plot heatmaps for indomethacin sexual dimorphic genes

06f: test differences in cell proportions by DPA

**Notebooks 07: scATAC analysis in HSCs and LSKs**

07a: preprocess and cluster HSC scATAC data and run ChromVar

07b: preprocess MPP scATAC dataset in R, extract counts from MACS2 peaks

07c: filter MPP scATAC dataset using Signac

07c: merge MPP and HSC scATAC dataset, this is now the LSK dataset

07e: cluster LSK scATAC dataset using Signac

07f: run ChromVar on LSK scATAC dataset

07g: plotting of scATAC data