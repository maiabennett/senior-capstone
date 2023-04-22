# An analysis of circulating, colon- and intestine-residing innate lymphoid cell (ILC) phenotypes and functions

This is the repository for Maia Bennett's Senior Capstone in Bioinformatics (BS) at the University of Nebraska at Omaha. This project, "Investigating the roles of innate immune cells in intestinal immunity and gut microbiome signaling," will be completed in May 2023.

Literature strongly suggests connections and interactions between ILCs and microbiota, specifically through metabolic activity and ILC plasticity. Given the important supporting role ILCs play in innate immunity, there is an urgent need to identify the mechanisms by which ILCs are influenced by microbiota. This is furthered by the heterogenicity of gut microbiota between individuals and the responsiveness of microbiota to various medications, disorders, and environments. Investigation of these interactions may lead to therapeutic targets in cases of dysregulated immunity or microbiome. These factors drive the question, how are innate lymphoid cell (ILC) plasticity and metabolic activity influenced by gut microbiota presence and signaling? Previous data supporting and informing these relationships have been gathered using mouse models; although they are a useful model system, there is a need for similar investigations using primary human data. Utilization of human data, especially scRNA-seq data, is likely to provide even greater insight into the biological and molecular interplay between ILCs and gut microbiota. For this project, I aim to utilize primary single-cell sequencing data to investigate metabolic activity and plasticity of ILCs relative to gut microbiota.

## The main goal of this project is to identify potential ILC-microbiome cross-regulation.
- The first aim of this project is to quantify the phenotypic composition of ILCs across three public human scRNA-seq datasets using Seurat (v4).
- The second aim of this project is to investigate phenotypic differences between ILCs in humans with normal and dysregulated (ulcerative colitis; UC) gut conditions using Seurat (v4).
- The third aim of this project is to establish potential involvement of ILCs in microbial processes using pathfindR to query KEGG gene sets.

## General overview of the analysis

![analysis_overview](https://user-images.githubusercontent.com/123126475/225156104-99060923-9e8e-4812-88e2-ed81447f32b2.png)

## Implementation
This GitHub includes all code used to generate the final results of the indicated senior capstone project. The order of all indicidual code components, as well as their contents and general descriptions of code utilization, can be found in [this R markdown file](https://github.com/maiabennett/senior-capstone/blob/main/senior-capstone.Rmd). Seurat objects were saved at each point of analysis (pre-processing and each iteration of clustering/annotation) to increase the reproducibility of this project without necessitating computationally heavy processes. These RDS files may be obtained by contacting Maia Bennett (maiabennett@unomaha.edu) and can be imported into R for analysis using the steps detailed [in this R script](https://github.com/maiabennett/senior-capstone/blob/main/readRDS.R). Any extraneous questions on implementation and use can also be directed to Maia Bennett at the above email. 
