## ViScoreR

ViScoreR is an R package which evaluates dimensionality reduction (DR) of single-cell data by comparing neighbourhood relations of labelled cell populations in the high-dimensional (HD) input data and the low-dimensional (LD) embedding, and between different LD embeddings.
Each DR tool introduces artifacts.
ViScoreR is a diagnostic tool that identifies them easily, so as to prevent faulty reasoning about data based on a misleading embedding.
(If you're interested in unsupervised evaluation, we recommend the use of R<sup>NX</sup> curves instead.)

Our metrics are computed at a cell population level, rather than having a single score for an embedding.

The first tool we implement is the xNPE (Extended Neighbourhood-Proportion-Error), which is, roughly speaking, a shape distortion measure.
It is based on the [Neighbourhood Proportion Error](https://github.com/akonstodata/NPE).
For each population (reference), the xNPE measures how the distribution of same-vs-differently labelled cells in the reference’s neighbourhood change, as we go from HD to LD.
**xNPE is intended for comparing shape distortions between different LD embeddings of the same HD data.**

The second tool are NCPs (neighbourhood composition plots).
They reveal which cell populations are represented in the near neighbourhood of a given population (reference), in HD data or in any given LD embedding of it.
**NCPs serve to identify sources of positional distortion.**

See the **Algorithms** section below for more information.

### Install

To install ViScoreR, run the following code in your R session.

```
devtools::install_github('saeyslab/ViScoreR')
```

ViScoreR depends on `tidyverse`, available in CRAN, and `emdist` and `BiocNeighbors`, available on Bioconductor.

### Usage

To view a 2-dimensional embedding of data, use `ViScoreR::PlotEmbedding`.

To compute *K*-nearest-neighbour graph for HD or LD data, use `ViScoreR::ComputeKNN`.

To compute the xNPE, use `ViScoreR::xNPE`.
Then use `ViScoreR::PlotLD` to inspect the likeness distributions for specific populations and `ViScoreR::PlotxNPE` to contrast the xNPE values of multiple embeddings of the same HD data.

To compute neighbourhood composition for a chosen population in the HD data or any LD embedding, use `ViScoreR::NeighbourhoodComposition`.
To plot the composition, use `ViScoreR::PlotNC`.

All plots are made with `ggplot2`, and as such can be modified (eg. to change font sizes or theme), combined with `gridExtra` or saved using the `ggplot2::ggsave` function.

### Documentation

Each ViScoreR exported function is documented.
Use the `?ViScoreR::function_name` syntax to view the help file for a function of interest.

### Shekhar retina dataset case study

A short case study is shown in the `poster.pdf` file, using the [Shekhar retina](https://rdrr.io/github/LTLA/scRNAseq/man/ShekharRetinaData.html) scRNA-seq dataset.

To download RDS files with the pre-processed HD data (50 principal components of the re-scaled count data), cell labels, a UMAP embedding and a t-SNE embedding, use [Git LFS](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage).

<hr>

### Algorithms

#### xNPE

For a matrix of coordinates *X*, the *K*-ary neighbourhood of each point in *X* is found.
The label for each respective identified neighbour are retrieved.

These neighbourhoods-per-point are then aggregated for each cell population.
Thus, each cell population *p* of size *N<sup>p</sup>* has *N<sup>p</sup>* vectors of length *K*.

For each cell in each population *p*, we record the number of neighbours that come from the same population, out of all the *K* neighbours.
Each population now has *N<sup>p</sup>* values (each between 0 and *K*).
We record these values as a distribution (the *likeness distribution*, referring to the number of 'like' cells).

These likeness distributions are calculated for each population, using a matrix of HD coordinates *X<sup>HD</sup>* and a matrix of LD coordinates *X<sup>LD</sup>*.
The xNPE score for a gives population is then the (Earth mover's distance)[https://en.wikipedia.org/wiki/Earth_mover%27s_distance] between its likeness distribution in HD and in LD.

The lower this value, the less distortion occurs in the process of the embedding.

#### NCP

For a select reference population *p* in a datasets (HD or LD), we look at the *K*-ary neighbourhood of each cell belonging to it.
We constrain the search by ignoring points belonging to *p* or to any chosen 'unassigned'/debris population.

The relative share of cells from each population in the data in different neighbourhood ranges (by rank) in the *K*-ary neighbourhood of the reference is then computed.

We produce a stacked area plot, where the y-axis represents the share of cells from a population (stacked by population, so that it adds up to 1) and the x-axis represents a segment of the neighbourhood, aggregated by a step size (of 10, by default).
(The plot is not cumulative along the x-axis.)

