<!-- AddToAny BEGIN -->
<div class="a2a_kit a2a_kit_size_32 a2a_default_style">
<a class="a2a_dd" href="//www.addtoany.com/share"></a>
<a class="a2a_button_facebook"></a>
<a class="a2a_button_twitter"></a>
<a class="a2a_button_google_plus"></a>
<a class="a2a_button_pinterest"></a>
<a class="a2a_button_reddit"></a>
<a class="a2a_button_sina_weibo"></a>
<a class="a2a_button_wechat"></a>
<a class="a2a_button_douban"></a>
</div>
<script async src="//static.addtoany.com/menu/page.js"></script>
<!-- AddToAny END -->

<link rel="stylesheet" href="https://guangchuangyu.github.io/css/font-awesome.min.css">

[![releaseVersion](https://img.shields.io/badge/release%20version-2.10.7-blue.svg?style=flat)](https://bioconductor.org/packages/DOSE)
[![develVersion](https://img.shields.io/badge/devel%20version-2.11.11-blue.svg?style=flat)](https://github.com/GuangchuangYu/DOSE)
[![total](https://img.shields.io/badge/downloads-45877/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/DOSE)
[![month](https://img.shields.io/badge/downloads-2240/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/DOSE)

This package implements five methods proposed by _Resnik_, _Schlicker_, _Jiang_, _Lin_ and _Wang_ respectively for measuring semantic similarities among DO terms and gene products. Enrichment analyses including hypergeometric model and gene set enrichment analysis are also implemented for discovering disease associations of high-throughput biological data. 


`DOSE` is released within the [Bioconductor](https://bioconductor.org/packages/DOSE) project and the source code is hosted in <a href="https://github.com/GuangchuangYu/DOSE"><i class="fa fa-github fa-lg"></i> GitHub</a>.


## <i class="fa fa-user"></i> Author

Guangchuang Yu, School of Public Health, The University of Hong Kong.

## <i class="fa fa-book"></i> Citation

Please cite the following article when using `DOSE`:

[![doi](https://img.shields.io/badge/doi-10.1093/bioinformatics/btu684-blue.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btu684)
[![citation](https://img.shields.io/badge/cited%20by-19-blue.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=16627502277303919270)
[![Altmetric](https://img.shields.io/badge/Altmetric-27-blue.svg?style=flat)](https://www.altmetric.com/details/2788597)

__Yu G__, Wang L, Yan G and He QY<sup>*</sup>. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. **Bioinformatics**, 2015, 31(4):608-609.


## <i class="fa fa-pencil"></i> Featured Articles

<img src="featured_img/c5mb00663e-f1_hi-res.gif" width="500">

<i class="fa fa-hand-o-right"></i> Find out more on <i class="fa fa-pencil"></i> [Featured Articles](https://guangchuangyu.github.io/DOSE/featuredArticles/).


## <i class="fa fa-download"></i> Installation

Install `DOSE` is easy, follow the guide in the [Bioconductor page](https://bioconductor.org/packages/DOSE/):

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("DOSE")
```

## <i class="fa fa-cogs"></i> Overview

#### <i class="fa fa-angle-double-right"></i> Semantic similarity measurement

+ DO term semantic similarity
+ Gene semantic similarity

#### <i class="fa fa-angle-double-right"></i> Enrichment Analysis

+ DO (Disease Ontology)
+ NCG (Network of Cancer Genes)
+ DisGeNet (gene-disease and SNP-disease associations)

#### <i class="fa fa-angle-double-right"></i> Visualization

+ barplot
+ cnetplot
+ dotplot
+ enrichMap
+ gseaplot
+ simplot
+ upsetplot


<i class="fa fa-hand-o-right"></i> Find out details and examples on <i class="fa fa-book"></i> [Documentation](https://guangchuangyu.github.io/DOSE/documentation/).

## <i class="fa fa-code-fork"></i> Projects that depend on _DOSE_


#### <i class="fa fa-angle-double-right"></i> Bioconductor packages
+ [bioCancer](https://www.bioconductor.org/packages/bioCancer): Interactive Multi-Omics Cancers Data Visualization and Analysis
+ [ChIPseeker](https://www.bioconductor.org/packages/ChIPseeker): ChIPseeker for ChIP peak Annotation, Comparison, and Visualization
+ [clusterProfiler](https://www.bioconductor.org/packages/clusterProfiler): statistical analysis and visualization of functional profiles for genes and gene clusters
+ [debrowser](https://www.bioconductor.org/packages/debrowser): Interactive Differential Expresion Analysis Browser
+ [eegc](https://www.bioconductor.org/packages/eegc): Engineering Evaluation by Gene Categorization (eegc)
+ [facopy](https://www.bioconductor.org/packages/facopy): Feature-based association and gene-set enrichment for copy number alteration analysis in cancer
+ [LINC](https://www.bioconductor.org/packages/LINC): annotation of noncoding RNAs and coexpressed genes
+ [ReactomePA](https://www.bioconductor.org/packages/ReactomePA): Reactome Pathway Analysis


## <i class="fa fa-comment"></i> Feedback
<ul class="fa-ul">
	<li><i class="fa-li fa fa-hand-o-right"></i> Please make sure you [follow the guide](https://guangchuangyu.github.io/2016/07/how-to-bug-author/) before posting any issue/question</li>
	<li><i class="fa-li fa fa-bug"></i> For bugs or feature requests, please post to <i class="fa fa-github-alt"></i> [github issue](https://github.com/GuangchuangYu/DOSE/issues)</li>
	<li><i class="fa-li fa fa-question"></i>  For user questions, please post to [Bioconductor support site](https://support.bioconductor.org/) and [Biostars](https://www.biostars.org/). We are following every post tagged with **DOSE**</li>
	<li><i class="fa-li fa fa-commenting"></i> Join the group chat on <a href="https://twitter.com/hashtag/DOSE"><i class="fa fa-twitter fa-lg"></i></a> and <a href="http://huati.weibo.com/k/DOSE"><i class="fa fa-weibo fa-lg"></i></a></li>
</ul>
