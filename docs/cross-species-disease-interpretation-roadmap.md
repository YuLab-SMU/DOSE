# Cross-species disease interpretation roadmap

Related execution plan: see `docs/deeprare-inspired-dose-development-plan.md` for the evidence-first, traceable reasoning, and optional LLM-assisted interpretation layer inspired by DeepRare.

## Positioning

DOSE should evolve from a disease ontology enrichment and semantic similarity
package into a cross-species disease interpretation framework. The central
scientific axis is:

```text
human disease <-> human disease genes <-> human phenotypes
        <-> human-mouse orthologs
mouse disease models <-> mouse genes <-> mouse phenotypes
```

The goal is not to become a single-cell, spatial transcriptomics, or mouse
database browser package. Those data types should be supported only through a
thin input adapter layer: users bring gene vectors, ranked gene lists, or named
gene-set lists derived from their own workflows, and DOSE maps those results
onto human disease, mouse model, phenotype, and orthology evidence.

The publishable story should therefore be:

> DOSE provides phenotype-aware cross-species disease interpretation that helps
> researchers connect human disease genetics with mouse models, explain disease
> similarity, and generate testable model or mechanism hypotheses from real
> omics-derived gene sets.

## Researcher needs

### Need 1: select and interpret mouse models for a human disease

Researchers often know a human disease of interest but need to choose, compare,
or justify mouse models. Existing resources can list models, but they usually do
not provide an integrated, explainable ranking that combines disease genes,
human-mouse orthology, human phenotypes, mouse phenotypes, and disease ontology
context.

DOSE should answer:

- Which mouse models are available for this human disease?
- Which models best capture the human disease phenotype profile?
- Which human disease genes and mouse orthologs support each model?
- Which disease phenotypes are covered or missing in each model?
- Which models are known annotations, and which are plausible candidate models?

### Need 2: infer human disease relevance from mouse data

Mouse researchers often have mouse genes, mouse phenotypes, or model-derived
gene signatures and want to know which human diseases they may inform.

DOSE should answer:

- Which human diseases are implicated by a mouse gene set or phenotype profile?
- Which human orthologs drive the signal?
- Are the implicated diseases supported by existing mouse disease model
  annotations?
- Are there model-supported candidate human disease genes that lack strong human
  evidence?

### Need 3: explain disease similarity across genes, phenotypes, and models

Existing semantic similarity methods can score disease or phenotype terms, but a
single score does not explain why two diseases are related.

DOSE should answer:

- Are two diseases similar because of ontology structure, shared genes, shared
  HPO phenotypes, shared mouse models, or shared MP phenotypes?
- Which specific genes and phenotypes explain the similarity?
- Does a shared mouse model provide a mechanistic bridge between the diseases?

### Need 4: connect user-derived gene sets to cross-species disease context

Users increasingly start from bulk RNA-seq, single-cell RNA-seq, spatial
transcriptomics, proteomics, or screening results. DOSE should not reimplement
those analysis pipelines, but it should accept their common outputs.

DOSE should accept:

- a gene vector
- a named ranked gene vector
- a named list of gene sets, for example cell types, clusters, contrasts, or
  spatial regions

DOSE should return:

- enriched or prioritized human diseases
- implicated mouse models
- driving genes and orthologs
- HPO and MP phenotype evidence
- evidence source summaries

This keeps real data integration practical without changing the package into an
omics workflow framework.

## Current DOSE capabilities

The current package already has several important foundations.

### Ontology semantic similarity

Current APIs:

- `doseSim()` / `doSim()` for disease ontology term similarity
- `geneSim()` for gene similarity based on ontology annotations
- `clusterSim()` and `mclusterSim()` for gene cluster similarity
- `computeIC()` and internal `semdata()` support through `GOSemSim`

Current ontologies:

- HDO / DO
- HPO
- MPO

### Enrichment and GSEA

Current APIs:

- `enrichDO()`
- `gseDO()`
- `enrichDisease()`
- `gseDisease()`
- `enrichNCG()`
- `gseNCG()`

Current data construction scripts:

- `data-raw/create-HDO.R`
- `data-raw/create-HPO.R`
- `data-raw/create-MPO.R`
- `data-raw/create-NCG.R`

HDO currently uses Alliance Genome Resource human disease annotations:

```text
DISEASE-ALLIANCE_HUMAN.tsv.gz -> DBObjectSymbol + DOID -> Entrez gene -> HDO.sqlite
```

HPO uses `phenotype_to_genes.txt`. MPO uses MGI phenotype annotations. NCG is a
cancer-specific gene set source.

### Useful implementation patterns

Existing code already supports:

- ontology-backed SQLite construction through `create_sqlite()`
- ontology-to-gene and gene-to-ontology retrieval through internal helpers
- `gson` construction for enrichment through `get_dose_data()`
- lazy cached annotation data in package environment

These patterns should be reused for cross-species data rather than inventing a
separate storage system.

## Core gaps

### Data gaps

DOSE does not yet have a unified data layer for:

- human-mouse orthologs
- mouse genotype or allele disease model annotations
- mouse model to MP phenotype profiles
- human disease to HPO phenotype profiles in a model-ranking-ready format
- cross-species disease-model evidence tables
- model-level entities distinct from gene-level or term-level entities

### Method gaps

DOSE currently supports ontology similarity and gene-set enrichment, but not:

- mouse model prioritization for a human disease
- human disease inference from mouse genes or MP terms
- explanatory decomposition of disease similarity
- evidence-aware disease interpretation combining multiple sources
- candidate disease model or candidate disease gene discovery

### Interface gaps

DOSE accepts gene vectors and ranked vectors for enrichment/GSEA, but not yet:

- a general interpretation API that returns both enrichment and cross-species
  evidence
- named gene-set lists from real workflows
- model-centric result classes or tidy result tables
- helper functions that expose ortholog and model evidence directly

## Data sources

### Required first wave

#### MGI

MGI should be the primary source for mouse disease model interpretation.

Target reports:

- `MGI_DO.rpt`: mouse gene to DO disease associations
- `MGI_Geno_DiseaseDO.rpt`: mouse genotypes with MP phenotype and DO disease
  annotations
- `MGI_DiseaseGeneModel.rpt`: human gene-centered mouse disease model report
- `MGI_DiseaseMouseModel.rpt`: human disease and mouse models by genotype
- `HOM_MouseHumanSequence.rpt`: mouse-human homology classes
- `HOM_ProteinCoding.rpt`: conservative protein-coding orthologs
- `MGI_GenePheno.rpt` or equivalent MP gene phenotype report

Primary use:

- disease -> mouse model lookup
- mouse gene -> human ortholog mapping
- model -> MP phenotype profile construction
- disease -> model -> gene/ortholog evidence

#### HPO

Required for human disease phenotype profiles.

Target files:

- HPO ontology OBO, already used by DOSE
- phenotype-to-gene, already used by DOSE
- disease-to-phenotype annotations, for example `phenotype.hpoa`

Primary use:

- disease -> HPO profile
- human disease phenotype coverage
- disease profile similarity against mouse model phenotypes

#### MPO / MP

Already used by DOSE, but should be extended for model-level phenotype profiles.

Primary use:

- mouse gene/model -> MP profile
- MP semantic similarity
- MP phenotype explanation for mouse models

#### HDO / DO

Already used by DOSE.

Primary use:

- disease ontology similarity
- existing human disease enrichment
- disease ID layer for MGI disease annotations

### Strong second wave

#### MONDO

MONDO should become the preferred disease normalization layer if cross-resource
integration grows beyond DO IDs.

Primary use:

- map DO, OMIM, Orphanet, EFO, and other disease IDs
- support cross-database disease unification
- improve disease-level interoperability

#### GenCC

Useful for clinical validity evidence on human gene-disease relationships.

Primary use:

- add confidence labels to human disease genes
- distinguish established disease genes from candidate genes

#### Open Targets

Useful for broader human target-disease evidence.

Primary use:

- score human disease-gene evidence
- add genetics, literature, known-drug, and other evidence channels
- support evidence-aware disease interpretation

#### IMPC

Useful for systematic knockout mouse phenotype data.

Primary use:

- expand mouse gene -> MP phenotype coverage
- discover candidate mouse models when curated disease model annotations are
  absent

### Lower priority for this story

These can be added later but should not drive the first cross-species paper:

- GWAS Catalog: valuable, but gene assignment from variants is not direct
- JensenLab DISEASES: useful enrichment source, less central to mouse model
  interpretation
- CTD: useful but automation and licensing details need care
- DisGeNET, OMIM, COSMIC: scientifically useful but limited by license or
  redistribution constraints

## Data products

Processed data should be generated by `data-raw` scripts and GitHub Actions,
then hosted on `gh-pages` or another release channel. Heavy processing
dependencies should not be package runtime dependencies.

### Storage convention

Use two canonical runtime formats.

#### Ontology-backed data: SQLite

If a data product contains an ontology or is used for semantic similarity, it
should be stored in the same SQLite style as current HDO/HPO/MPO data.

This applies to:

- HDO / DO
- HPO
- MPO / MP
- future MONDO
- future EFO, if supported
- any ontology that needs term metadata, parent-child relationships, ancestors,
  information content, or semantic similarity

Expected content:

```text
term table
relationship / graph table
ontology metadata
optional term-to-gene table
optional term-to-all-gene table
```

Design rule:

If a resource has a DAG and will be passed to semantic similarity functions, it
should be built with the ontology SQLite pipeline rather than stored as a plain
gene-set file.

#### Enrichment-only data: gson

If a data product is only used as a gene-set collection for ORA or GSEA, it
should be stored as `gson`.

This applies to:

- NCG
- JensenLab DISEASES, if added
- GenCC disease-gene sets, if used only for enrichment
- Open Targets disease-gene sets, if used only for enrichment
- GWAS Catalog trait-gene sets, if added
- custom mouse model gene sets when no ontology semantics are required

Expected content:

```text
gsid2gene
gsid2name
species
keytype
gsname
version
accessed_date
source metadata
```

Design rule:

If a resource does not provide ontology structure, or if DOSE will not use that
structure for semantic similarity, it should be exposed as `gson` and routed
through the existing `enrichDisease()` / `gseDisease()` style machinery.

#### Auxiliary evidence tables: compressed TSV or RDS

Some cross-species evidence does not fit either ontology SQLite or enrichment
gene-set formats. These tables should remain separate, versioned auxiliary data
products.

Examples:

- ortholog mappings
- mouse genotype/model metadata
- model-to-phenotype evidence
- disease-to-model evidence
- reference and evidence source tables
- score component tables

Default format:

```text
*.tsv.gz
```

Use RDS only if the object structure is materially richer than a table and has a
stable R-side schema.

Design rule:

Auxiliary evidence tables should not be forced into `gson`, because they are not
only term-to-gene mappings. They should instead be loaded by dedicated helper
functions and joined into interpretation results.

### Minimal first-wave data products

#### Ortholog table

File:

```text
mouse_human_ortholog.tsv.gz
```

Columns:

```text
mouse_gene_id
mouse_symbol
human_gene_id
human_symbol
orthology_type
homology_class_id
source
```

Notes:

- `orthology_type` should distinguish `one_to_one`, `one_to_many`,
  `many_to_many`, or `unknown`.
- Default analyses should use one-to-one orthologs unless the user opts into
  broad mapping.

#### Mouse disease model table

File:

```text
mouse_disease_model.tsv.gz
```

Columns:

```text
model_id
model_label
disease_id
disease_name
mouse_gene_id
mouse_symbol
human_gene_id
human_symbol
allele_id
allele_symbol
genotype_id
genotype_label
strain_background
source
reference_id
```

Notes:

- Prefer a source MGI genotype or model identifier for `model_id`.
- In the canonical processed table, `model_id` is required, non-missing, and
  immutable within a data release. If a source does not provide a suitable
  identifier, the build script must create a deterministic namespaced ID from
  the source genotype identifier and record that derivation in metadata.
- One row represents one `model_id`-disease-gene-allele association. Phenotypes
  remain in the separate model phenotype table; joins must not multiply a
  model's score by its number of genes, alleles, references, or phenotype rows.

#### Mouse model phenotype table

File:

```text
mouse_model_phenotype.tsv.gz
```

Columns:

```text
model_id
mp_id
mp_name
evidence_code
sex
zygosity
source
reference_id
```

#### Human disease phenotype table

File:

```text
human_disease_phenotype.tsv.gz
```

Columns:

```text
disease_id
disease_name
hpo_id
hpo_name
qualifier
evidence_code
onset
frequency
sex
modifier
aspect
source
reference_id
```

`qualifier` must be retained from `phenotype.hpoa`. Negated annotations such as
`NOT` must never enter a positive phenotype profile. They should be retained as
conflicting or absent evidence so that explanations can distinguish an
explicitly absent phenotype from an unannotated phenotype.

#### Human disease gene evidence table

File:

```text
human_disease_gene_evidence.tsv.gz
```

Columns:

```text
disease_id
disease_name
human_gene_id
human_symbol
source
evidence_type
evidence_score
confidence_label
```

Sources can include current HDO/Alliance annotations first, then GenCC and Open
Targets later.

#### Mouse gene phenotype table

File:

```text
mouse_gene_phenotype.tsv.gz
```

Columns:

```text
mouse_gene_id
mouse_symbol
mp_id
mp_name
source
evidence_code
```

This supports user input that is gene-based rather than model-based.

### Relational and release contract

Every processed release must include a machine-readable manifest containing:

- data release and schema versions
- source URLs, source release identifiers, access dates, and licenses
- SHA-256 checksums for source and processed files
- row counts, primary keys, foreign keys, and expected join cardinalities
- the build script commit and ontology versions

Loaders must verify checksums, use atomic cache writes, and fail with an
actionable message when a cached file does not match its manifest. Unit tests
should use small versioned fixtures; package checks must not depend on network
availability or on mutable latest-release URLs.

### Optional indexes

Precompute indexes only if runtime cost becomes high:

- disease -> model list
- model -> MP profile list
- disease -> HPO profile list
- disease -> human gene list
- human gene -> mouse ortholog list
- MP/HPO term ancestor closures

## Analysis functions

### `rankMouseModels()`

Purpose:

Prioritize mouse models for a human disease.

Prototype:

```r
rankMouseModels <- function(
    disease,
    disease_id_type = c("auto", "DOID", "MONDO", "OMIM", "Orphanet", "EFO"),
    mode = c("lookup", "discovery"),
    orthology = c("one_to_one", "all"),
    evidence = c("disease_annotation", "gene_overlap", "phenotype_similarity"),
    weights = NULL,
    min_phenotype_matches = 1,
    top = 50
)
```

Inputs:

- a disease ID or disease name
- optional ID type
- orthology strictness
- scoring components and weights
- `lookup` mode for evidence-supported model selection or `discovery` mode for
  evaluating candidate models without using the target disease-model annotation

`evidence` selects channels to expose in the returned object. `weights` applies
only to validated numeric discovery components; it cannot assign a weight to a
known disease-model annotation.

Outputs:

A tidy result with:

```text
disease_id
disease_name
model_id
model_label
score
known_disease_annotation
gene_overlap_score
phenotype_similarity_score
phenotype_coverage
mouse_gene
human_ortholog
matched_hpo_terms
matched_mp_terms
missing_hpo_terms
source
references
```

Scoring:

- known disease annotation: reported as evidence in `lookup` mode, but never
  used as a discovery score or as a feature in discovery benchmarks
- gene overlap score: overlap between disease genes and model human orthologs
- phenotype similarity score: HPO disease profile versus MP model profile
- phenotype coverage: proportion of disease HPO profile matched by model MP
  profile

First implementation can use gene overlap and known model annotations. Phenotype
similarity can be added as a second iteration if HP-MP cross-ontology mapping is
not ready.

### `inferHumanDisease()`

Purpose:

Infer human disease relevance from mouse genes, human genes, MP terms, or
derived gene sets.

Prototype:

```r
inferHumanDisease <- function(
    x,
    input = c("auto", "human_gene", "mouse_gene", "mp", "hpo"),
    orthology = c("one_to_one", "all"),
    ontology = c("HDO", "MONDO"),
    method = c("ora", "gsea", "profile_similarity"),
    ...
)
```

Supported inputs:

- human gene vector
- mouse gene vector
- named ranked vector
- MP or HPO term vector
- named list of gene sets

Outputs:

```text
query
disease_id
disease_name
score_or_pvalue
adjusted_pvalue
matched_genes
mouse_orthologs
matched_models
matched_hpo_terms
matched_mp_terms
evidence_sources
```

This function should be the main entry point for user-derived gene sets from
bulk, single-cell, spatial, proteomics, or screen data.

### `explainDiseaseModel()`

Purpose:

Explain why a mouse model is relevant to a human disease.

Prototype:

```r
explainDiseaseModel <- function(
    disease,
    model,
    include = c("genes", "orthologs", "phenotypes", "evidence", "references")
)
```

Outputs:

- matched disease genes
- mouse genes and human orthologs
- HPO disease profile
- MP model profile
- matched or semantically similar HPO-MP phenotype pairs
- missing HPO phenotypes
- supporting references and data sources

This should be optimized for interpretability rather than statistical testing.

### `explainDiseaseSimilarity()`

Purpose:

Explain similarity between two human diseases.

Prototype:

```r
explainDiseaseSimilarity <- function(
    disease1,
    disease2,
    components = c("ontology", "genes", "phenotypes", "models")
)
```

Outputs:

```text
disease1
disease2
overall_similarity
ontology_similarity
gene_jaccard
phenotype_similarity
shared_model_score
shared_genes
shared_hpo_terms
shared_mouse_models
shared_mp_terms
```

This extends the current semantic similarity story with evidence decomposition.

### `geneDiseaseProfile()`

Purpose:

Summarize cross-species evidence for a candidate gene.

Prototype:

```r
geneDiseaseProfile <- function(
    gene,
    organism = c("human", "mouse"),
    orthology = c("one_to_one", "all")
)
```

Outputs:

- human diseases associated with the gene
- clinical validity evidence, when available
- mouse orthologs
- mouse MP phenotypes
- mouse disease models involving the gene
- candidate diseases suggested by mouse phenotype similarity

## Canonical public API

### `interpretDisease()`

`interpretDisease()` is the single high-level entry point for gene-, ranked
gene-, gene-set-, and phenotype-based workflows. Other analysis functions are
advanced entry points and must return the same result class.

```r
interpretDisease <- function(
    x,
    input = c("auto", "gene", "ranked_gene", "gene_set_list", "phenotype"),
    organism = c("human", "mouse"),
    analysis = c("disease", "mouse_model", "both"),
    ontology = "auto",
    id_type = "auto",
    orthology = c("one_to_one", "all"),
    explain = c("none", "template", "llm"),
    top = 50,
    ...
)
```

User-facing behavior:

- `character` input is a gene vector unless every value has an `HP:` or `MP:`
  prefix; prefixed vectors are deterministically recognized as phenotypes.
- a named numeric vector is a ranked gene vector and a named list is a set of
  gene queries. Empty, duplicated, unrecognized, or mixed identifier inputs
  produce actionable validation messages rather than silent coercion.
- `organism = "human"` is the default. Gene symbols and Entrez IDs do not encode
  species, so `organism = "auto"` is intentionally not supported for genes.
- `ontology = "auto"` chooses the valid knowledge sources for the requested
  organism, input, and analysis. An explicit ontology is accepted only when it
  is compatible with that combination.
- `id_type = "auto"` recognizes ontology prefixes, Ensembl prefixes, and
  all-numeric Entrez IDs. Other non-empty character genes are treated as symbols
  within the declared organism; mixed identifier systems require an explicit
  supported key type or produce an actionable error.
- `explain = "none"` never loads an LLM dependency. Template and LLM
  explanations are derived from the same structured evidence rows.

Compatibility and conversion rules:

| Input | Organism | Direct annotation | Cross-species path |
|---|---|---|---|
| genes or ranked genes | human | HDO, HPO, NCG | human gene -> mouse ortholog -> MPO/model |
| genes or ranked genes | mouse | MPO | mouse gene -> human ortholog -> HDO/HPO/NCG |
| HPO terms | human | HPO disease profiles | HPO -> curated HP-MP mapping -> mouse model |
| MP terms | mouse | MPO model profiles | MP -> curated MP-HP mapping -> human disease |

An invalid combination, such as mouse genes with direct HDO enrichment or HPO
terms declared as mouse input, must fail before analysis. Ortholog conversion
must be recorded as evidence and must never be inferred merely by relabeling an
enrichment result's organism field.

### Advanced functions and accessors

`rankMouseModels()`, `inferHumanDisease()`, `explainDiseaseModel()`,
`explainDiseaseSimilarity()`, and `geneDiseaseProfile()` remain public for
focused workflows. They share validation and scoring code with
`interpretDisease()` and return `doseInterpretResult`; they do not define
parallel result classes.

Required inspection methods are:

- `as.data.frame()` for the ranked result table
- `evidence()` for filterable long-form evidence
- `summary()` and `show()` for query, source, and uncertainty summaries
- `plot()` or focused plotting helpers that consume the same object

Canonical accessor and explanation signatures:

```r
evidence <- function(
    x,
    target_id = NULL,
    type = NULL,
    source = NULL,
    direction = NULL
)

explainDisease <- function(
    x,
    top = 5,
    method = c("template", "llm"),
    provider = NULL,
    model = NULL,
    language = c("en", "zh"),
    style = c("brief", "research", "cautious")
)
```

`interpretDisease(explain = ...)` is a convenience path that calls the same
explanation implementation. `explainDisease()` is useful when explanation is
requested after inspecting or filtering an existing result; it must not rerun
or modify the analysis.

## Canonical result contract

Use one Bioconductor-style S4 class, `doseInterpretResult`, from the first public
release. It contains `result`, `evidence`, `query`, `sources`, `parameters`, and
`explanation` slots. All constructors validate the following invariants:

- each result row is uniquely identified by `query_id`, `target_type`, and
  `target_id`
- `score` is accompanied by `score_type`, `score_direction`, and visible score
  components; p-values are `NA` for non-inferential rankings
- evidence rows have stable `evidence_id` values and link to a valid result row
- source release, schema version, and data checksum are retained in `sources`
- unsupported, missing, conflicting, and ambiguous evidence remain distinct

Required result fields:

```text
query_id, target_id, target_name, target_type, rank,
score, score_type, score_direction, pvalue, p.adjust,
evidence_count, top_evidence_type, source_count
```

Required evidence fields:

```text
evidence_id, query_id, target_id, target_type, evidence_type,
direction, feature_id, feature_name, component_score,
source, source_record_id, reference_id, note
```

The class contract is normative. The DeepRare-inspired execution plan and all
vignettes must reference this section rather than restating a second signature
or result schema.

## Visualization

Visualization should support interpretation, not become the main feature.

First wave:

- dotplot for diseases or models
- barplot for score components
- upset-style or cnet-style gene/model overlap using existing ecosystem tools
- heatmap of disease-model scores

Second wave:

- evidence network:

```text
disease -> HPO terms
disease -> human genes -> mouse orthologs -> mouse model -> MP terms
```

Avoid hard dependencies on Seurat, Giotto, or SpatialExperiment. Vignettes can
show how users pass marker or DEG lists from those tools.

## Implementation plan

### Phase 0: design and data audit

Deliverables:

- confirm exact MGI report columns and stable identifiers
- confirm HPO disease annotation file and disease ID mapping
- audit the existing `data-raw/mh-mapping.R` prototype and select a versioned,
  licensed MP-HP SSSOM mapping source
- define supported ID systems for phase 1
- document data licenses and redistribution constraints
- decide whether phase 1 uses only DO IDs or includes MONDO normalization
- freeze the canonical API, compatibility matrix, result schema, primary keys,
  and score semantics defined in this roadmap

Acceptance criteria:

- all first-wave data sources can be downloaded by script
- all processed data products have stable schemas
- fixtures demonstrate that negated HPO annotations are not treated as support
- manifest validation detects checksum, schema, and foreign-key failures
- invalid organism, ontology, and identifier combinations fail before analysis
- package runtime dependencies remain close to current dependencies

### Phase 1: core data layer

Deliverables:

- `data-raw/create-MGI-ortholog.R`
- `data-raw/create-MGI-mouse-model.R`
- `data-raw/create-HPO-disease-phenotype.R`
- processed TSV or RDS files hosted externally
- loader functions with local caching

Internal helpers:

```r
get_ortholog_data()
get_mouse_model_data()
get_mouse_model_phenotype()
get_human_disease_phenotype()
get_human_disease_gene_evidence()
```

Acceptance criteria:

- disease -> models lookup works
- model -> MP terms lookup works
- human gene -> mouse ortholog lookup works
- mouse gene -> human ortholog lookup works
- data version and source metadata are recorded
- loaders verify the release manifest and cache files atomically
- joins preserve one score contribution per intended evidence entity

### Phase 2: first user-facing cross-species functions

Deliverables:

- `rankMouseModels()`
- `inferHumanDisease()`
- `explainDiseaseModel()`
- `geneDiseaseProfile()`

Initial scoring can be gene/evidence based:

```text
score = weighted_sum(
    disease_gene_overlap,
    ortholog_support
)
```

Known model annotations are evidence labels, not numeric discovery features.
Raw phenotype count must not be used as a score because it rewards annotation
density. Phase 3 may add normalized phenotype similarity and information-content
weighted coverage after the cross-ontology mapping is validated. Every component
must define its range, missing-value behavior, normalization, and weight before
scores from different queries are compared.

Acceptance criteria:

- examples run without large runtime dependencies
- results include evidence columns, not just scores
- functions work with both human and mouse gene IDs
- all public analysis functions return `doseInterpretResult`

### Phase 3: phenotype-aware scoring

Deliverables:

- disease HPO profile support
- model MP profile support
- phenotype profile similarity score
- phenotype coverage report

Open technical question:

HPO and MP are separate ontologies. Cross-ontology phenotype matching requires
one of:

- curated HP-MP mappings
- Monarch/uPheno integration
- lexical and ontology-mediated approximation as a fallback

Preferred approach:

- use established cross-species phenotype mappings when available
- keep fallback matching explicitly marked as approximate

Acceptance criteria:

- `rankMouseModels()` can rank by phenotype-aware score
- `explainDiseaseModel()` reports matched and missing phenotype evidence
- benchmark reports the effect of phenotype-aware ranking relative to gene-only
  ranking with uncertainty, including a null or negative result

### Phase 4: thin real-data adapter

Deliverables:

- `interpretDisease()`
- support for named gene-set lists
- vignette examples for:
  - bulk DEG list
  - single-cell cell-type markers
  - spatial region marker genes

Design rule:

No direct dependency on single-cell or spatial packages in core package code.
Examples can show optional extraction code, but the API should accept plain R
vectors and lists.

Acceptance criteria:

- a user can pass `list(CellTypeA = genes, CellTypeB = genes)`
- output includes per-query disease and mouse model interpretation
- results share one schema across input sources; numeric scores are compared
  only when `score_type`, component definitions, and normalization match

### Phase 5: evidence expansion

Deliverables:

- GenCC support for clinical validity evidence
- Open Targets support for human target-disease evidence
- MONDO support for cross-resource disease ID normalization
- optional IMPC support for systematic knockout phenotype expansion

Acceptance criteria:

- evidence sources are clearly tracked
- default results remain interpretable
- users can filter by evidence type or confidence

## Benchmark and publication plan

### Benchmark 1: known model recovery

Task:

Given a human disease, rank candidate mouse models and test whether held-out MGI
disease models are recovered near the top. The annotation being predicted must
not be available to candidate generation, scoring, feature construction, or
parameter tuning.

Evaluation design:

- prefer a temporal split: build features from release N and evaluate model
  annotations first appearing in release N+1
- if temporal snapshots are unavailable, use disease-stratified annotation
  holdout and remove the held-out disease-model edges before deriving features
- define the candidate universe before evaluation and include all eligible
  models, not only annotated positives or models returned by direct lookup
- tune weights on separate diseases or nested folds; report results by disease
  and annotation density, with bootstrap confidence intervals
- deduplicate at `disease_id`-`model_id` before computing metrics

Baselines:

- gene overlap only
- ontology similarity only
- phenotype count only
- direct MGI lookup, reported separately as database coverage rather than as a
  discovery ranking baseline

Metrics:

- AUROC
- AUPRC
- top-k recall
- mean reciprocal rank

### Benchmark 2: disease inference from mouse gene sets

Task:

Given genes from held-out mouse disease models, recover the annotated human
disease. The target model-disease edge and any gene-disease evidence derived
from that edge must be removed before constructing the query and disease
features. Prefer the same temporal snapshots as Benchmark 1; otherwise use
model- and disease-grouped folds so related rows cannot cross the split.

Baselines:

- human ortholog disease enrichment only
- MP enrichment only
- random ortholog-matched gene sets

Metrics:

- top-k disease recovery
- rank of true disease
- calibration by gene-set size

### Benchmark 3: phenotype-aware improvement

Task:

Compare model ranking with and without HPO-MP phenotype profile similarity.
Both variants must use the same held-out labels, candidate universe, gene
features, and tuning split. HP-MP mappings and phenotype annotations derived
from the evaluated disease-model edge must be excluded or documented as a
potential circularity.

Claim rule:

Claim improvement only if the held-out comparison and confidence intervals
support it. Otherwise report where phenotype evidence helps interpretation,
where it has no measurable ranking benefit, and how annotation coverage limits
the result.

### Case studies

Pick diseases where human disease genes, HPO profiles, and mouse models are
well-annotated. Good candidates include:

- cystic fibrosis
- Alzheimer disease
- Parkinson disease
- diabetes mellitus
- muscular dystrophy
- retinal degeneration
- inflammatory bowel disease

Each case study should show:

- ranked mouse models
- score component decomposition
- matched human genes and mouse orthologs
- matched and missing phenotype terms
- connection to a user-derived gene set when available

## Success criteria

The upgrade should be considered successful when DOSE can answer these questions
with reproducible data and interpretable output:

- Given a human disease, which mouse models are most relevant and why?
- Given mouse genes or phenotypes, which human diseases are implicated?
- Given two diseases, what genes, phenotypes, and models explain their
  similarity?
- Given a user-derived gene set from any omics workflow, what cross-species
  disease-model hypotheses does it support?

The package should still feel like DOSE: ontology-aware, enrichment-capable, and
bioconductor-friendly. The new contribution is the cross-species interpretation
layer that connects disease ontology analysis to mouse model biology.
