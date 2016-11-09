DOSE: Disease Ontology Semantic and Enrichment analysis
=======================================================

[![releaseVersion](https://img.shields.io/badge/release%20version-3.0.1-green.svg?style=flat)](https://bioconductor.org/packages/DOSE) [![develVersion](https://img.shields.io/badge/devel%20version-3.1.1-green.svg?style=flat)](https://github.com/GuangchuangYu/DOSE) [![Bioc](http://www.bioconductor.org/shields/years-in-bioc/DOSE.svg)](https://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#since) [![total](https://img.shields.io/badge/downloads-49256/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/DOSE) [![month](https://img.shields.io/badge/downloads-2523/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/DOSE)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![codecov](https://codecov.io/gh/GuangchuangYu/DOSE/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/DOSE/) [![Last-changedate](https://img.shields.io/badge/last%20change-2016--11--09-green.svg)](https://github.com/GuangchuangYu/DOSE/commits/master) [![commit](http://www.bioconductor.org/shields/commits/bioc/DOSE.svg)](https://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#svn_source) [![GitHub forks](https://img.shields.io/github/forks/GuangchuangYu/DOSE.svg)](https://github.com/GuangchuangYu/DOSE/network) [![GitHub stars](https://img.shields.io/github/stars/GuangchuangYu/DOSE.svg)](https://github.com/GuangchuangYu/DOSE/stargazers)

[![platform](http://www.bioconductor.org/shields/availability/devel/DOSE.svg)](https://www.bioconductor.org/packages/devel/bioc/html/DOSE.html#archives) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/DOSE.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/DOSE/) [![Linux/Mac Travis Build Status](https://img.shields.io/travis/GuangchuangYu/DOSE/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/DOSE) [![AppVeyor Build Status](https://img.shields.io/appveyor/ci/Guangchuangyu/DOSE/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/DOSE) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-green.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-dose/README.html)

This package implements five methods proposed by *Resnik*, *Schlicker*, *Jiang*, *Lin* and *Wang* respectively for measuring semantic similarities among DO terms and gene products. Enrichment analyses including hypergeometric model and gene set enrichment analysis are also implemented for discovering disease associations of high-throughput biological data.

[![Twitter](https://img.shields.io/twitter/url/https/github.com/GuangchuangYu/DOSE.svg?style=social)](https://twitter.com/intent/tweet?hashtags=DOSE&url=http://bioinformatics.oxfordjournals.org/content/31/4/608)

------------------------------------------------------------------------

Please cite the following article when using `DOSE`:

***G Yu***, LG Wang, GR Yan, QY He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. ***Bioinformatics*** 2015, 31(4):608-609.

[![doi](https://img.shields.io/badge/doi-10.1093/bioinformatics/btu684-green.svg?style=flat)](http://dx.doi.org/10.1093/bioinformatics/btu684) [![citation](https://img.shields.io/badge/cited%20by-21-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=16627502277303919270) [![Altmetric](https://img.shields.io/badge/Altmetric-36-green.svg?style=flat)](https://www.altmetric.com/details/2788597)

------------------------------------------------------------------------

For details, please visit our project website, <https://guangchuangyu.github.io/DOSE>.

-   [Documentation](https://guangchuangyu.github.io/DOSE/documentation/)
-   [Featured Articles](https://guangchuangyu.github.io/DOSE/featuredArticles/)
-   [Feedback](https://guangchuangyu.github.io/DOSE/#feedback)

### Citation

[![citation](https://img.shields.io/badge/cited%20by-21-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=16627502277303919270)

         +-+---------+---------+--------+---------+--------+---+
      12 +                                                 *   +
         |                                                     |
    11.5 +                                                     +
         |                                                     |
      11 +                                                     +
         |                                                     |
    10.5 +                                                     +
      10 +                                                     +
         |                                                     |
     9.5 +                                                     +
         |                                                     |
       9 + *                                                   +
         +-+---------+---------+--------+---------+--------+---+
         2015     2015.2    2015.4   2015.6    2015.8    2016   

### Download stats

[![download](http://www.bioconductor.org/shields/downloads/DOSE.svg)](https://bioconductor.org/packages/stats/bioc/DOSE/) [![total](https://img.shields.io/badge/downloads-49256/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/DOSE) [![month](https://img.shields.io/badge/downloads-2523/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/DOSE)

    ggtree v1.7.3  For help: https://guangchuangyu.github.io/ggtree

    If you use ggtree in published research, please cite:
    Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution 2016, doi:10.1111/2041-210X.12628

         +-----------+---------------+---------------+---------------+--------------+---------------+--+
    2500 +                                                                                        *    +
         |                                                                                             |
         |                                                                                 *  **       |
         |                                                                                  *          |
    2000 +                                                                                *     *      +
         |                                                                                             |
         |                                                                                             |
         |                                                                                             |
         |                                                                                             |
    1500 +                                                                            * *              +
         |                                                                    * *  *                   |
         |                                                                        * *  *               |
         |                                                                     *                       |
    1000 +                                                                ***                          +
         |                                                       *     **                              |
         |                                                **   *  **  *                                |
         |                                  *     *      *   **      *                                 |
     500 +                               **  * * * * **                                                +
         |               ** *     * *         *        *                                               |
         |                   ** **   ** *                                                              |
         |          ** **                                                                              |
       0 +   * *** *                                                                                   +
         +-----------+---------------+---------------+---------------+--------------+---------------+--+
                   2012            2013            2014            2015           2016            2017
