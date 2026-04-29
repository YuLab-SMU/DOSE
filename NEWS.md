# DOSE 4.6.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# DOSE 4.5.1

+ mv `get_organism` to 'GOSemSim' (2025-12-07, Sun)
+ update enrichment functions to use 'enrichit' (2025-12-05, Fri)
+ mv `enricher_internal` and `GSEA_internal` to the 'enrichit' package (2025-12-05, Fri)
  - renamed to `ora_gson` and `gsea_gson` with updates
  - mv helper functions to 'enrichit'

# DOSE 4.4.0

+ Bioconductor RELEASE_3_22 (2025-11-01, Sat)

# DOSE 4.2.0

+ Bioconductor RELEASE_3_21 (2025-04-17, Thu)

# DOSE 4.0.0

+ Bioconductor RELEASE_3_20 (2024-10-30, Wed)

# DOSE 3.99.1

+ return NULL in GSEA if not genes can be mapped (2024-08-26, Mon, ReactomePA#43)

# DOSE 3.31.4

+ remove `mpoSim()` and `hopSim()` (2024-08-21, Wed)

# DOSE 3.31.3

+ remove `gseMPO()` and `gseHPO()`, instead using e.g., `gseDO(ont="HPO")` (2024-08-15, Thu)
+ remove `enrichMPO()` and `enrichHPO()`, instead using e.g., `enrichDO(ont="HPO")`
+ unify API and removing the usages of MPO.db and HPO.db 
+ update DO-gene mapping data (2024-08-13, Tue)
+ remove `HDO.db` and use new API in GOSemSim (v>=2.31.1) (2024-08-13, Tue)
+ use `yulab.utils::yulab_msg()` for startup message (2024-07-26, Fri)

# DOSE 3.31.2

+ add `FoldEnrichment`, `RichFactor` and `zScore` in ORA result (2024-06-13, Thu)

# DOSE 3.31.1

+ fixed bug in `options(enrichment_force_universe=TRUE)` (2024-05-16, Thu)
  - <https://github.com/YuLab-SMU/clusterProfiler/issues/283#issuecomment-2077869230>

# DOSE 3.30.0

+ Bioconductor RELEASE_3_19 (2024-05-15, Wed)

# DOSE 3.29.2

+ fix bugs in `get_ont_info()` and `get_dose_data()`: wrong object name and wrong data type (2023-11-30, Thu)

# DOSE 3.29.1

+ mv 'MPO.db' and 'HPO.db' from 'Imports' to 'Suggests' and fixed bugs (2023-11-18, Sat)

# DOSE 3.28.0

+ Bioconductor RELEASE_3_18 (2023-10-25, Wed)

# DOSE 3.27.3

+ update `TERM2NAME()` to return term if corresponding name not found. (2023-10-09, Mon)

# DOSE 3.27.2

+ use 'MPO.db' and 'HPO.db' to support phenotype ontology for mouse and human (2023-06-30, Fri)

# DOSE 3.27.1

+ `options(enrichment_force_universe = TRUE)` will force enrichment analysis to intersect the `universe` with gene sets (2023-05-03, Wed) 
+ use `inherits` to judge the class of objects (2022-11-20, Sun)
+ test whether slot in `GSON` object is NULL (e.g., `GSON@keytype`) when assigning it to enrichment result (2022-11-07, Mon)

# DOSE 3.26.0

+ Bioconductor RELEASE_3_17  (2023-05-03, Wed)

# DOSE 3.24.0

+ Bioconductor RELEASE_3_16 (2022-11-02, Wed)

# DOSE 3.23.3

+ replace `DO.db` to `HDO.db` (2022-10-7, Fri)  
+ add values of  `organism`, `keytype` and `setType` for `GSEA_internal()` (2022-09-21, Wed)
+ add values of  `organism`, `keytype` and `ontology` for `enricher_internal()` (2022-09-21, Wed)
+ move `inst/extdata/parse-obo.R` to `HDO.db` package (2022-08-29, Mon) 
+ rename `qvalues` to `qvalue` in `gseaResult` object (2022-08-29, Mon)

# DOSE 3.23.2

+ Support `GSON` object in `GSEA_internal()` (2022-06-08, Wed)

# DOSE 3.23.1

+ Support `GSON` object in `enricher_internal()` (2022-06-06, Mon)

# DOSE 3.22.0

+ Bioconductor 3.15 release

# DOSE 3.21.2

+ enable `setReadable` for compareCluster(GSEA algorithm) result(2021-12-13, Mon)
+ update the default order of GSEA result (2021-12-09, Thu)
  - if p.adjust is identical, sorted by `abs(NES)`


# DOSE 3.21.1

+ upate DisGeNET and NCG data (2021-11-14, Sun)
  - DisGeNET v7: 21671 genes, 30170 diseases and 1134942 gene-disease associations
    - 194515 variants, 14155 diseases and 369554 variant-disease associations
  - NCG v7: 3177 cancer genes, 130 diseases and 6095 gene-disease associations

# DOSE 3.20.0

+ Bioconductor 3.14 release

# DOSE 3.19.4

+ update `clusterProfiler` citation (2021-09-30, Thu)
+ upate error message of `enricher_internal` (2021-9-3, Fri)

# DOSE 3.19.3

+ upate DisGeNET and NCG data (2021-8-16, Mon)

# DOSE 3.19.2

+ bug fixed, change 'is.na(path2name)' to 'all(is.na(path2name))' (2021-06-21, Mon)

# DOSE 3.19.1

+ add `dr` slot to `compareClusterResult`, `enrichRestul` and `gseaResult`(2021-5-21, Fri)

# DOSE 3.18

+ Bioconductor 3.13 release

# DOSE 3.17

+ support setting seed for fgsea method if e.g. `gseGO(seed = TRUE)` (2020-10-28, Wed)
  - <https://github.com/YuLab-SMU/DOSE/issues/45>
  
# DOSE 3.16.0

+ Bioconductor 3.12 release (2020-10-28, Wed)

# DOSE 3.15.4

+ update `setReadable` and `geneInCategory` methods for `compareClusterResult` object (2020-10-12, Mon)

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
