# DOSE 3.15.3

+ allow passing additional parameters to fgsea (2020-10-09, Fri)
  - <https://github.com/YuLab-SMU/DOSE/pull/40>
+ add `termsim` and `method` slots to `compareClusterResult`, `enrichRestul` and `gseaResult`
  - <https://github.com/YuLab-SMU/DOSE/pull/39>

# DOSE 3.15.2

+ update [NCG](http://ncg.kcl.ac.uk/download.php#) and [DGN](https://www.disgenet.org/downloads) data (2020-10-09, Thu)

# DOSE 3.14.0

+ Bioconductor 3.11 release

# DOSE 3.13.2

+ fixed issue caused by R v4.0.0 (2020-03-12, Thu)
  - length > 1 in coercion to logical
  - <https://github.com/YuLab-SMU/DOSE/pull/32>

# DOSE 3.13.1

+ remove `S4Vectors` dependencies (2019-12-19, Thu)
+ extend `setReadable` to support `compareClusterResult` (2019-12-02, Mon)
+ add `gene2Symbol`, `keytype` and `readable` slots for `compareClusterResult`
+ move `compareClusterResult` class definition from `clusterProfiler` (2019-11-01, Fri)

# DOSE 3.12.0

+ Bioconductor 3.10 release

# DOSE 3.11.2

+ ignore `universe` and print a message if users passing accidentally passing wrong input (2019-10-24, Thu)
  - <https://github.com/YuLab-SMU/clusterProfiler/issues/217>
+ gene with minimal ES value (NES < 0) will be reported in `core_enrichment` (2019-07-31, Wed)

# DOSE 3.11.1

+ `build_Anno` now compatible with `tibble` (2019-05-28, Tue)

# DOSE 3.10.0

+ Bioconductor 3.9 release

# DOSE 3.9.4

+ export `parse_ratio` (2019-03-29, Tue)

# DOSE 3.9.4

+ bug fixed of `get_enriched` (2019-01-14, Mon)
  - <https://github.com/GuangchuangYu/clusterProfiler/issues/177>

# DOSE 3.9.2

+ mv enrichment vignettes to [clusterProfiler-book](https://yulab-smu.github.io/clusterProfiler-book) (2019-01-10, Thu)

# DOSE 3.9.1

+ `asis` parameter in `[.enrichResult` and `[.gseaResult` (2018-12-24, Mon)
  - <https://github.com/GuangchuangYu/enrichplot/issues/17>

# DOSE 3.8

+ Bioconductor 3.8 release

# DOSE 3.7.1

+ S3 accessor methods only return enriched terms. (2018-06-20, Wed)
