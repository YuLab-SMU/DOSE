# DOSE: Disease Ontology Semantic and Enrichment analysis

[![](https://img.shields.io/badge/release%20version-4.6.0-green.svg)](https://www.bioconductor.org/packages/DOSE)
[![](https://img.shields.io/badge/devel%20version-4.7.3-green.svg)](https://github.com/YuLab-SMU/DOSE)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/DOSE.svg)](https://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#since)
[![codecov](https://codecov.io/gh/GuangchuangYu/DOSE/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/DOSE/)

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![platform](http://www.bioconductor.org/shields/availability/devel/DOSE.svg)](https://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#archives)
[![Build
Status](http://www.bioconductor.org/shields/build/devel/bioc/DOSE.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/DOSE/)
[![](https://img.shields.io/badge/download-1463455/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/DOSE)
[![](https://img.shields.io/badge/download-28438/month-blue.svg)](https://bioconductor.org/packages/stats/bioc/DOSE)

Provides ontology-aware methods for disease and phenotype knowledge
mining. DOSE supports semantic similarity analysis of disease and
phenotype ontology terms, genes, and gene clusters using methods
including Resnik, Schlicker, Jiang, Lin, and Wang. It also provides
over-representation analysis and gene set enrichment analysis for
interpreting gene vectors and ranked gene lists in disease, phenotype,
and cancer contexts.

## :writing_hand: Authors

Guangchuang YU <https://yulab-smu.top>

School of Basic Medical Sciences, Southern Medical University

Learn more at <https://yulab-smu.top/contribution-knowledge-mining/>.

## :arrow_double_down: Installation

Get the released version from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DOSE")
```

Or install the development version from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("YuLab-SMU/DOSE")
```

Please cite the following article when using `DOSE`:

Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
R/Bioconductor package for Disease Ontology Semantic and Enrichment
analysis. Bioinformatics. 2015, 31(4):608-609.

## :sparkling\_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
