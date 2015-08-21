# Disease Ontology Semantic and Enrichment analysis

[![platform](http://www.bioconductor.org/shields/availability/devel/DOSE.svg)](http://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#archives)
[![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/DOSE.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/DOSE/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/DOSE.svg)](http://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#since)
[![post](http://www.bioconductor.org/shields/posts/DOSE.svg)](https://support.bioconductor.org/t/DOSE/)
[![commit](http://www.bioconductor.org/shields/commits/bioc/DOSE.svg)](http://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#svn_source)
[![download](http://www.bioconductor.org/shields/downloads/DOSE.svg)](http://bioconductor.org/packages/stats/bioc/DOSE.html)


 This package implements five methods proposed by Resnik, Schlicker, Jiang, Lin and Wang respectively for measuring semantic similarities among DO terms and gene products. Enrichment analyses including hypergeometric model and gene set enrichment analysis are also implemented for discovering disease associations of high-throughput biological data. 

## Authors ##

Guangchuang YU, School of Public Health, The University of Hong Kong [http://ygc.name](http://ygc.name)

## Citation ##

Please cite the following article when using `DOSE`:

```
Yu G, Wang L, Yan G and He Q.
DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis.
Bioinformatics, 2015, 31(4):608-609.
```

URL: [http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btu684](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btu684)

## License ##

All source code is copyright, under the Artistic-2.0 License.
For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Installation ##

To install:
 * the latest released version:
   `biocLite("DOSE")`
 * the latest development version:
   `install_github("GuangchuangYu/DOSE")`

## Documentation ##


+ [Why clusterProfiler fails](http://ygc.name/2014/08/07/why-clusterprofiler-fails/)
+ [NCG enrichment implemented in DOSE <- bioinfoblog.it](http://bioinfoblog.it/2015/04/ncg-enrichment-implemented-in-dose/)
+ [DisGeNET enrichment analysis via clusterProfiler](http://ygc.name/2015/05/11/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/)
+ [Enrichment Map](http://ygc.name/2014/08/03/enrichment-map/)
+ [dotplot for enrichment result](http://ygc.name/2015/06/23/dotplot-for-enrichment-result/)
+ [functional enrichment for GTEx paper](http://ygc.name/2015/08/13/functional-enrichment-for-gtex-paper/)
+ [functional enrichment analysis with NGS data](http://ygc.name/2015/08/21/functional-enrichment-analysis-with-ngs-data/)


Find out more at [http://www.bioconductor.org/packages/release/bioc/html/DOSE.html](http://www.bioconductor.org/packages/release/bioc/html/DOSE.html) and check out the [vignettes](http://www.bioconductor.org/packages/release/bioc/vignettes/DOSE/inst/doc/DOSE.pdf).

To view the vignette of `DOSE` installed in your system, start R and enter:
```r
vignette("DOSE", package="DOSE")
```

## Bugs/Feature requests ##

 - If you have any, [let me know](https://github.com/GuangchuangYu/DOSE/issues). Thx!

