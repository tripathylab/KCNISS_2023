# KCNI Summer School 2023
## Project 1 - Single Cell Transcriptomics

**Integrative analysis of human SMART-seq and 10x single-cell gene expression data.** How do the methodologies compare in defining cell types, and how can we use single-cell data for case-control analyses?

### What’s this project about? 

**Main idea:** perform integrative analysis of two human neocortical single-cell datasets, which will lead into an analysis of cell type-specific differential gene expression between individuals with dementia and controls.

**Key questions:**

1. How does gene expression compare between cell types defined by SMART-seq vs. those defined by 10x sequencing?
    1. Are certain cell types better defined by one method relative to the other?
2. How does cell type-specific gene expression change in the context of dementia?
    1. What types of genes are most affected by the condition within a given cell type?

**What datasets are available to help answer these questions?**
[Allen Institute for Brain Sciences Cell Types database](https://celltypes.brain-map.org/)
* [Human Multiple Cortical Areas SMART-seq](https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq)
* [Human MTG 10x SEA-AD](https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad)

**Transcriptomics TAs:** Mel Davie, Derek Howard

## Resources

Seurat tutorials:
* [Dataset integration workflow](https://satijalab.org/seurat/articles/integration_introduction.html)
* [Differential expression analysis](https://satijalab.org/seurat/articles/de_vignette.html)
* [Data visualization methods](https://satijalab.org/seurat/articles/visualization_vignette.html)
* [Cell type annotation](https://satijalab.org/seurat/articles/integration_mapping.html)

[Pseudobulk differential expression tutorial tutorial](https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html)

## Schedule
**Day 1:** [Intro to single-cell RNAseq analysis, R, and Seurat](https://github.com/tripathylab/KCNISS_2023/blob/main/code/Day_1.Rmd)

**Day 2:** [Intro to differential expression, cell type identification & visualizations](https://github.com/tripathylab/KCNISS_2023/blob/main/code/Day_2.Rmd)

**Day 3:** [Dataset integration & automated cell type annotation](https://github.com/tripathylab/KCNISS_2023/blob/main/code/Day_3.Rmd)

**Day 4:** [Case-control differential expression with pseudobulks](https://github.com/tripathylab/KCNISS_2023/blob/main/code/Day_4.Rmd)

**Day 5:** [Spatial biology talks @ SickKids & final presentation!](https://github.com/tripathylab/KCNISS_2023/blob/main/code/Day_5.Rmd)
