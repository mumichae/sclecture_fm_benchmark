# Benchmarking Single Cell Foundation models beyond cell type

## Background

Large-scale single-cell studies across multiple labs have increased the available single-cell data across tissues, organs, and organisms. This wealth of data enables us to investigate single-cell biology at a population-scale level by building single-cell atlases, which form the basis for large-scale single-cell foundation models (scFMs).
As more single-study atlases and scFMs emerge, the challenge of batch integration persists. Datasets from different institutions contain strong technical biases that mask the biological signal of interest. Various modeling approaches address this by disentangling technical from biological variation, ranging from linear models to cVAE models to large-scale transformer-based scFMs. However, while the task is inherently unsupervised, model performance evaluations often prioritize the preservation of cell type annotations at the expense of other biological signals. Expanding the quantitative assessment of batch integration is especially important for scFMs, where we expect the models to learn more complex biological signals.
Further research is needed to characterize biological signals beyond cell types in single-cell representations for more nuanced evaluation of scFMs. Using informative genes as ground truth or quantifying local structure differences e.g. through topological dissimilarity across different embeddings can add a purely data-driven perspective to the quantitative evaluation of models.

## Objective

The aim is to set up a large-scale benchmark of scFMs and compare the global and local structure of the different representations. It also provides the opportunity of getting involved in an actively developed software project.
Impact: You will evaluate the performance of scFMs on batch integration beyond cell type annotations using novel quantitative approaches. The results of this project will contribute to a paper that is currently in progress. Your contribution will be acknowledged in the paper, with the chance of co-authorship if desired.

## Datasets

* [Single-cell atlas of the human retina v1.0](https://data.humancellatlas.org/hca-bio-networks/eye/atlases/retina-v1-0)

## Research questions

1. How are single cell foundation models currently evaluated and what are the limitations?
2. How do scFMs compare to each other (and task-specific baselines) in batch correction and biological conservation, when considering local structure?

## Methods

You will be using and developing the scAtlasTb codebase as well as investigating different scFM representations globally and locally using existing metrics. The scAtlasTb abstracts the coding overhead of complex benchmark setups, providing a framework for scalable computation of state-of-the art batch integration methods and model performance evaluation with label-based and label-free evaluation metrics.
You will be working with the following tools:

* https://github.com/pytorch/pytorch
* https://github.com/helicalAI/helical
* https://github.com/scverse/scanpy
* https://github.com/scverse/anndata
* https://github.com/snakemake/snakemake
* https://github.com/HCA-integration/scAtlasTb


## Learning outcomes

* large-scale benchmarking of foundation models
* quantitative and qualitative comparison of single cell representations
* quantify local structure of real-world single cell data beyond data annotations
* principles of good collaborative software design on a real-world software project

## References

1. Luecken, M.D., Büttner, M., Chaichoompu, K. *et al. “*Benchmarking atlas-level data integration in single-cell genomics”. Nat Methods 19, 41–50 (2022). https://doi.org/10.1038/s41592-021-01336-8
2. Baek, S., Song, K. & Lee, I. “Single-cell foundation models: bringing artificial intelligence into cell biology.” Exp Mol Med 57, 2169–2181 (2025). https://doi.org/10.1038/s12276-025-01547-5
3. Li, J., Wang, J., Ibarra, I.L., et al. “Integrated multi-omics single cell atlas of the human retina.” bioRxiv 2023.11.07.566105; doi: https://doi.org/10.1101/2023.11.07.566105
4. Kedzierska, K.Z., Crawford, L., Amini, A.P. et al. Zero-shot evaluation reveals limitations of single-cell foundation models. Genome Biol 26, 101 (2025). https://doi.org/10.1186/s13059-025-03574-x
5. Wang, H., Leskovec, J. & Regev, A. “Limitations of cell embedding metrics assessed using drifting islands.” Nat Biotechnol (2025). https://doi.org/10.1038/s41587-025-02702-z