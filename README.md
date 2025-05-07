
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fetchNOMe

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8402785.svg)](https://doi.org/10.5281/zenodo.8402785)

[![R-CMD-check](https://github.com/fmi-basel/gpeters-fetchNOMe/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fmi-basel/gpeters-fetchNOMe/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Overview

**fetchNOMe** is an R package for fast and efficient retrieval of
NOMe-seq data from BAM files. NOMe-seq (Nucleosome Occupancy and
Methylome Sequencing) simultaneously captures endogenous CpG methylation
and M.CviPI-induced GpC methylation, providing a high-resolution map of
chromatin accessibility and DNA methylation.

This package is specifically designed to process BAM files generated
from pipelines aligning bisulfite-converted DNA, such as Bismark, QuasR
and BISCUIT aligners.

## Features

- Extract GpC (M.CviPI-induced) and CpG (endogenous) methylation
  directly from BAM files (`get_data_matrix_from_bams`,
  `get_molecule_data_list_from_bams`)
- Extract co-occurrence statistics from BAM files for footprint spectral
  analysis by the **nomeR** package (`get_ctable_from_bams`)
- Extract statistics of protected states (`get_protect_stats_from_bams`)

## Requirements

- Indexed BAM files from NOMe-seq experiments
- Index (.bai) file in the same location as the BAM

## Installation

You can install the development version of fetchNOMe from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fmi-basel/gpeters-fetchNOMe")
```

## Usage

``` r
library(fetchNOMe)

## Define BAM file path
bam <- system.file("extdata","QuasR_test.bam",package = "fetchNOMe")

## Define genome FASTA file path
genome <- system.file("extdata","random_genome_700bp.fa",package = "fetchNOMe")

## extract GCH protection states for a region of interest (ROI)
GCH_tbl_mat <- get_data_matrix_from_bams(bamfiles = bam,
                                                                                         samplenames = "quasr_pe",
                                                                                         regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                                                                                         strand = "+",
                                                                                                                                                         IRanges::IRanges(start = 200,end = 500)),
                                                                                         genome = genome,
                                                                                         whichContext ="GCH",
                                                                                         remove_nonunique = F,
                                                                                         clip_until_nbg = 0L,
                                                                                         max_bisC_meth = 1)
                                                                                         
## extract WCG endogenous methylation states for a region of interest (ROI)
WCG_tbl_mat <- get_data_matrix_from_bams(bamfiles = bam,
                                                                                         samplenames = "quasr_pe",
                                                                                         regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                                                                                         strand = "+",
                                                                                                                                                         IRanges::IRanges(start = 200,end = 500)),
                                                                                         genome = genome,
                                                                                         whichContext ="WCG",
                                                                                         remove_nonunique = F,
                                                                                         clip_until_nbg = 0L,
                                                                                         max_bisC_meth = 1)

## extract GCH protection for all molecules within ROI
GCH_tbl_methlist <- get_molecule_data_list_from_bams(bamfiles = bam,
                                                     samplenames = "quasr_pe",
                                                     regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                      strand = "+",
                                                                                      IRanges::IRanges(start = 200,end = 500)),
                                                     genome = genome,
                                                     whichContext ="GCH",
                                                     remove_nonunique = F,
                                                     clip_until_nbg = 0L,
                                                     max_bisC_meth = 1,
                                                     min_frag_data_len = 50L,
                                                     min_frag_data_dens = 0.05)
                                                     
## extract co-occurrence count tables for footprint spectral analysis
ctable <- get_ctable_from_bams(bamfiles = bam,
                                             samplenames = "quasr_pe",
                                             regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                              strand = "+",
                                                                              IRanges::IRanges(start = 200,end = 500)),
                                             genome = genome,
                                             alignerUsed = "QuasR",
                                             min_frag_data_len = 0L,
                                             min_frag_data_dens = 0.0,
                                             max_spacing = 300,
                                             remove_nonunique = F,
                                             clip_until_nbg = 0L,
                                             max_bisC_meth = 1)
                                                                            
```

## Citation

If you use `fetchNOMe` in your work, please cite:

> Evgeniy A. Ozonov. **fetchNOMe**: fast retrieval of NOMe-seq data from
> BAM files (<https://doi.org/10.5281/zenodo.8402785>). Available at:
> <https://github.com/fmi-basel/gpeters-fetchNOMe>

## License

This package is licensed under the MIT License. Copyright (c) 2023
Friedrich Miescher Institute for Biomedical Research
