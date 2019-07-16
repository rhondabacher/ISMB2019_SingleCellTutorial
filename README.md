## ISMB 2019 Tutorial on "Recent Advances in Statistical Methods and Computational Algorithms for Single-Cell Omics Analysis"

This tutorial is focused on advanced statistical and computational methods that are recently developed for 
single-cell omics data. It is intended for an audience with genomics/computational background, who are
interested in cutting-edge developments of single-cell research, including both method development and application. 

All tutorial materials and extra links are provided here.

### Schedule

Official schedule for this tutorial can be found [here](https://www.iscb.org/ismbeccb2019-program/tutorials#am2).

### Course Instructors (equal contributors; alphabetical order)

Rhonda Bacher, University of Florida, rbacher@ufl.edu

Yuchao Jiang, University of North Carolina-Chapel Hill, yuchaoj@email.unc.edu

Jingshu Wang, University of Chicago, wangjingshususan@gmail.com

Please contact any of us regarding comments or questions.

### Tutorial Feedback

If you attend our tutorial at ISMB 2019, please provide feedback via this survery:
https://goo.gl/forms/0sR1kfVO6nj4X8bO2


### Slides

* [Introduction to Single-Cell RNA-seq and Pre-processing](https://github.com/rhondabacher/ISMB2019_SingleCellTutorial/blob/master/slides/1_introduction.pdf)

* [Single-Cell Visualization, Alignment, and Denoising](https://github.com/rhondabacher/ISMB2019_SingleCellTutorial/blob/master/slides/2_visualization_alignment_denoising.pdf)

* [Single-Cell Trajectory Inference](https://github.com/rhondabacher/ISMB2019_SingleCellTutorial/blob/master/slides/3_pseudotime.pdf)

* [Single-Cell Immune Profiling](https://github.com/rhondabacher/ISMB2019_SingleCellTutorial/blob/master/slides/4_immunology.pdf)

* [Single-Cell ATAC-Seq and Multimodal Alignment](https://github.com/rhondabacher/ISMB2019_SingleCellTutorial/blob/master/slides/5_multimodal_alignment.pdf)

* [Single-Cell Cancer Genomics](https://github.com/rhondabacher/ISMB2019_SingleCellTutorial/blob/master/slides/6_cancer_genomics.pdf)

### List of Methods


#### Single-cell quality control

* scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R ([paper](https://www.ncbi.nlm.nih.gov/pubmed/28088763), [software](http://bioconductor.org/packages/release/bioc/html/scater.html))


#### Single-cell normalization

* SCnorm: robust normalization of single-cell RNA-seq data ([paper](https://www.ncbi.nlm.nih.gov/pubmed/28418000), [software](https://bioconductor.org/packages/release/bioc/html/SCnorm.html))

* scran: Pooling across cells to normalize single-cell RNA sequencing data with many zero counts ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7), [software](http://bioconductor.org/packages/release/bioc/html/scran.html))

* scTransform: Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression ([paper](https://www.biorxiv.org/content/10.1101/576827v1), [software](https://github.com/ChristophH/sctransform))

#### Single-cell visualization

* t-SNE: t-Distributed Stochastic Neighbor Embedding ([paper](https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf))

* UMAP: Uniform Manifold Approximation and Projection ([paper](https://arxiv.org/pdf/1802.03426.pdf))
[software](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)


#### Single-cell batch correction

* mnnCorrect: Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors ([paper](https://www.ncbi.nlm.nih.gov/pubmed/29608177), [software](http://bioconductor.org/packages/release/bioc/html/scran.html))

* multiCCA: ([paper1](https://www.ncbi.nlm.nih.gov/pubmed/29608179), [paper2](https://www.cell.com/cell/pdf/S0092-8674(19)30559-8.pdf), [software](https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html))


#### Denoising
* SAVER: gene expression recovery for single-cell RNA sequencing ([paper](https://www.nature.com/articles/s41592-018-0033-z), [software](https://github.com/mohuangx/SAVER))

* DCA: Single-cell RNA-seq denoising using a deep count autoencoder ([paper](https://www.nature.com/articles/s41467-018-07931-2), [software](https://github.com/theislab/dca))

* scVI: Deep generative modeling for single-cell transcriptomics ([paper](https://www.nature.com/articles/s41592-018-0229-2), [software](https://github.com/YosefLab/scVI))

* SAVER-X: Transfer learning in single-cell transcriptomics improves data denoising and pattern discovery ([paper](https://www.biorxiv.org/content/10.1101/457879v2), [software](https://github.com/jingshuw/SAVERX))


#### Transfer learning
* SAVER-X: see above

* scGen: Generative modeling and latent space arithmetics predict single-cell perturbation response across cell types, studies and species ([paper](https://www.biorxiv.org/content/10.1101/478503v1), [software](https://github.com/theislab/scgen))

* cTP-net: Surface protein imputation from single cell transcriptomes by deep neural networks ([paper](https://www.biorxiv.org/content/10.1101/671180v1.full), [software](https://github.com/zhouzilu/cTPnet))


#### Single-cell pseudotime

* Full list at: [https://github.com/agitter/single-cell-pseudotime](https://github.com/agitter/single-cell-pseudotime)

* TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis ([paper](https://www.ncbi.nlm.nih.gov/pubmed/27179027), [software](https://github.com/zji90/TSCAN))

* Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics ([paper](https://www.ncbi.nlm.nih.gov/pubmed/29914354), [software](https://bioconductor.org/packages/release/bioc/html/slingshot.html))

* Monocle2/3: Reversed graph embedding resolves complex single-cell trajectories ([paper](https://www.ncbi.nlm.nih.gov/pubmed/28825705), [software](https://github.com/cole-trapnell-lab/monocle-release))

* Benchmarking: [Saelens et al., “A comparison of single-cell trajectory inference methods”. Nature Biotechnology. 2019.](https://www.ncbi.nlm.nih.gov/pubmed/30936559)


#### Single-cell clustering

* SINCERA: A Pipeline for Single-Cell RNA-Seq Profiling Analysis ([paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004575), [software](https://github.com/xu-lab/SINCERA))

* pcaReduce: hierarchical clustering of single cell transcriptional profiles ([paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0984-y), [software](https://github.com/JustinaZ/pcaReduce))

* CIDR: Ultrafast and accurate clustering through imputation for single-cell RNA-seq data ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1188-0), [software](https://github.com/VCCRI/CIDR))

* SNN-Cliq: Identification of cell types from single-cell transcriptomes using a novel clustering method. ([paper](https://www.ncbi.nlm.nih.gov/pubmed/25805722), [software](https://github.com/BIOINSu/SNN-Cliq))

* SOUP: Semisoft clustering of single-cell data ([paper](https://www.pnas.org/content/116/2/466), [software](https://github.com/lingxuez/SOUPR))

* SC3: consensus clustering of single-cell RNA-seq data ([paper](https://www.nature.com/articles/nmeth.4236), [software](https://bioconductor.org/packages/release/bioc/html/SC3.html))



#### Single-cell differential features

* SCDE: Bayesian approach to single-cell differential expression analysis ([paper](https://www.nature.com/articles/nmeth.2967), [software](http://bioconductor.org/packages/release/bioc/html/scde.html))

* MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data ([paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4676162/), [software](http://bioconductor.org/packages/release/bioc/html/MAST.html))

* BASiCS: Bayesian Analysis of Single-Cell Sequencing Data ([paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004333), [software](https://bioconductor.org/packages/release/bioc/html/BASiCS.html))

* DECENT: differential expression with capture efficiency adjustmeNT for single-cell RNA-seq data ([paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz453/5514046), [software](https://github.com/cz-ye/DECENT))

* scDD: A statistical approach for identifying differential distributions in single-cell RNA-seq experiments ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y), [software](https://bioconductor.org/packages/release/bioc/html/scDD.html))

* DESCEND: Gene expression distribution deconvolution in single-cell RNA sequencing ([paper](https://www.pnas.org/content/115/28/E6437), [software](https://github.com/jingshuw/descend))


#### Single-cell classification

* SingleCellNet: a computational tool to classify single cell RNA-Seq data across platforms and across species ([paper](https://www.biorxiv.org/content/10.1101/508085v1), [software](https://github.com/pcahan1/singleCellNet/))

* ACTINN: Automated identification of Cell Types in Single Cell RNA Sequencing ([paper](https://www.biorxiv.org/content/10.1101/532093v1), [software](https://github.com/mafeiyang/ACTINN))


#### Single-cell immune profiling

* TraCeR: T cell fate and clonality inference from single-cell transcriptomes ([paper](https://www.nature.com/articles/nmeth.3800),[software](https://github.com/Teichlab/tracer))

* VDJPuzzle: B-cell receptor reconstruction from single-cell RNA-seq ([paper](https://www.ncbi.nlm.nih.gov/pubmed/29659703),[software](https://bitbucket.org/kirbyvisp/vdjpuzzle2))



#### Single-cell epigenomics

* SCRAT: Single-cell regulome data analysis ([paper](https://academic.oup.com/bioinformatics/article/33/18/2930/3823309), [software](https://github.com/zji90/SCRAT))

* scABC: Unsupervised clustering and epigenetic classification of single cells ([paper](https://www.nature.com/articles/s41467-018-04629-3),[software](https://github.com/SUwonglab/scABC))

* Destin: Toolkit for single-cell analysis of chromatin accessibility ([paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz141/5367832),[software](https://github.com/urrutiag/destin))

#### Single-cell multimodal alignment

* PECA: Modeling gene regulation from paired expression and chromatin accessibility data ([paper](https://www.pnas.org/content/114/25/E4914),[software](http://web.stanford.edu/~zduren/PECA/))

* MATCHER: Manifold alignment reveals correspondence between single cell transcriptome and epigenome dynamics
([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1269-0),[software](https://github.com/jw156605/MATCHER))

* CCA (Seurat): Comprehensive Integration of Single-Cell Data ([paper](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub),[software](https://satijalab.org/seurat/))

#### Single-cell cancer genomics

* Canopy: Assessing intratumor heterogeneity and tracking longitudinal and spatial clonal evolutionary history by next-generation sequencing ([paper](https://www.pnas.org/content/113/37/E5528),[software](https://github.com/yuchaojiang/Canopy))

* MARATHON: Integrative pipeline for profiling DNA copy number and inferring tumor phylogeny ([paper](https://academic.oup.com/bioinformatics/article/34/12/2126/4838234), [software](https://github.com/yuchaojiang/MARATHON))

* InferCNV: Inferring CNV from Single-Cell RNA-Seq ([paper](https://science.sciencemag.org/content/344/6190/1396.long), [software](https://github.com/broadinstitute/infercnv))

* HoneyBADGER: Linking transcriptional and genetic tumor heterogeneity through allele analysis of single-cell RNA-seq data ([paper](https://genome.cshlp.org/content/early/2018/06/13/gr.228080.117), [software](https://github.com/JEFworks/HoneyBADGER))

* Cardelino: Integrating whole exomes and single-cell transcriptomes to reveal phenotypic impact of somatic variants ([paper](https://www.biorxiv.org/content/10.1101/413047v2), [software](https://github.com/PMBio/cardelino))

* SCOPE: A normalization and copy number estimation method for single-cell DNA sequencing ([paper](https://www.biorxiv.org/content/10.1101/594267v1), [software](https://github.com/rujinwang/SCOPE))



### Resources

#### Other tutorials and workflows
* Analysis of single cell RNA-seq data: https://hemberg-lab.github.io/scRNA.seq.course

* Awesome single cell: https://github.com/seandavi/awesome-single-cell

* simpleSingleCell: http://bioconductor.org/packages/simpleSingleCell


#### Datasets
* 10X Genomics: https://support.10xgenomics.com

* Single Cell Portal: https://portals.broadinstitute.org/single_cell

* Single Cell Expression Atlas: https://www.ebi.ac.uk/gxa/sc

* conquer: http://imlspenticton.uzh.ch:3838/conquer

* Human Cell Atlas: https://www.humancellatlas.org/

* Mouse Cell Atlas: http://bis.zju.edu.cn/MCA 

* JingleBells: http://jinglebells.bgu.ac.il

* scQuery: https://scquery.cs.cmu.edu/processed_data/

#### Specific reviews

* Cancer: [Baslan, T., & Hicks, J. (2017). Unravelling biology and shifting paradigms in cancer with single-cell sequencing. Nature Reviews Cancer, 17(9), 557.](https://www.ncbi.nlm.nih.gov/pubmed/28835719)
* Immunology: [Papalexi, E., & Satija, R. (2018). Single-cell RNA sequencing to explore immune cell heterogeneity. Nature Reviews Immunology, 18(1), 35.](https://www.ncbi.nlm.nih.gov/pubmed/28787399)
* Technology: [Kolodziejczyk, A. A., Kim, J. K., Svensson, V., Marioni, J. C., & Teichmann, S. A. (2015). The technology and biology of single-cell RNA sequencing. Molecular cell, 58(4), 610-620.](https://www.ncbi.nlm.nih.gov/pubmed/26000846)
* Design and Methods Overview: [Bacher, R., & Kendziorski, C. (2016). Design and computational analysis of single-cell RNA-sequencing experiments. Genome biology, 17(1), 63.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4823857/)
