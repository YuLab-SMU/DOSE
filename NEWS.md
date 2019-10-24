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
