# DeepRare-inspired DOSE development plan

## Why this plan exists

DeepRare is useful for DOSE less because it uses an LLM, and more because it treats structured medical knowledge as the main diagnostic substrate. The useful idea for DOSE is therefore not to build an AI doctor. The useful idea is to make disease, phenotype, gene, ortholog, and model evidence traceable, computable, and explainable.

The current cross-species roadmap defines the scientific direction: phenotype-aware disease interpretation across human disease, human genes, HPO, mouse orthologs, mouse models, and MPO. This plan defines the execution layer inspired by DeepRare:

- evidence-first result objects
- traceable reasoning tables
- high-level interpretation functions
- optional LLM summarisation through `aisdk`
- strict separation between computed evidence and generated prose

The guiding rule is simple: DOSE computes the evidence; an LLM may explain it, but must not invent it.

## Product position

DOSE should evolve from enrichment and semantic similarity into a knowledge-grounded disease interpretation package.

The package should help users answer:

- Which diseases are supported by my genes, ranked genes, phenotypes, or gene-set lists?
- Which genes and phenotypes drive each disease ranking?
- Which evidence source supports each result?
- Which mouse models or orthologs connect the result to cross-species evidence?
- What is missing, ambiguous, or only weakly supported?

This is deliberately not clinical diagnosis. The user-facing language should stay in the research interpretation space: disease interpretation, model prioritisation, evidence summary, and hypothesis generation.

## Design principles

### 1. Evidence before explanation

Every user-facing explanation should come from a structured evidence table. If a claim cannot be traced to a row in that table, it should not appear in the explanation.

### 2. Small interface, deep implementation

Users should not need to manually compose enrichment, semantic similarity, ortholog mapping, phenotype profiles, and model evidence. A small number of high-level functions should hide that complexity.

### 3. LLM as an adapter, not the engine

`aisdk` support should be optional. It should receive only structured DOSE output and return prose. It should not rank diseases, fetch untracked external facts, or change evidence scores.

### 4. Keep core package checks stable

The core package should work without network access, API keys, or LLM packages.
LLM features should live behind optional dependencies. Only live provider smoke
tests are skipped when unavailable; deterministic adapter contract tests always
run.

### 5. Preserve current DOSE strengths

Existing ontology-backed SQLite data, `gson` enrichment data, semantic similarity methods, and lazy cached annotation patterns should be reused rather than replaced.

## Proposed interfaces

### `interpretDisease()`

High-level entry point for real-world disease interpretation. Its normative
signature, input inference rules, organism-ontology compatibility matrix, and
return contract are defined in the **Canonical public API** section of
`docs/cross-species-disease-interpretation-roadmap.md`. This execution plan must
not maintain a second signature.

```r
interpretDisease(genes, organism = "human")
interpretDisease(ranked_genes, organism = "mouse")
interpretDisease(hpo_terms, input = "phenotype", organism = "human") # reserved until phenotype gates pass
```

Initial scope:

- character vector of genes
- named numeric ranked gene vector
- named list of gene sets
- HPO or MPO term vector only after the main roadmap's disease-mapping,
  phenotype-mapping, coverage, and provenance gates pass; before that,
  phenotype input is not a public `interpretDisease()` path and must return an
  actionable "not yet enabled" error if explicitly requested against the
  reserved canonical signature

All paths return the canonical `doseInterpretResult` with:

- ranked disease or model result table
- long evidence table
- query metadata
- source metadata
- parameters
- optional explanation text

`interpretDisease()` should call lower-level functions. It should not duplicate enrichment or similarity logic.

Until the cross-species target contract is implemented in Phase 4,
`target = "both"` is reserved but unavailable. Early releases must fail with an
actionable error instead of silently degrading to disease-only ranking.

### `evidence()`

Accessor for traceable evidence. Its signature is defined in the main roadmap.

```r
evidence(result, target_id = "DOID:example", type = "gene")
```

Expected filters:

- disease ID
- mouse model ID
- evidence type
- evidence source

This should be the main inspection interface for users who want to verify the reasoning.

### `explainDisease()`

Human-readable explanation layer. Its signature is defined in the main roadmap;
this plan specifies only the adapter behavior.

```r
explainDisease(result, method = "template", style = "research")
```

Rules:

- `method = "template"` must work without `aisdk`.
- `method = "llm"` should require optional `aisdk` configuration.
- LLM prompts must include only structured evidence generated by DOSE.
- The prompt should explicitly prohibit adding unsupported facts.
- The returned object should record provider, model, evidence hash, prompt hash,
  response hash, generation parameters, source versions, and generation time.

### `as.data.frame()`, `summary()`, `show()`

These methods should make the result useful without LLM support.

`summary()` should report:

- query type
- top diseases or models
- strongest evidence types
- missing evidence, if available
- data source versions

## Result object

The implementation uses the single S4 `doseInterpretResult` contract defined in
the main roadmap. A temporary list class is not an accepted public alternative,
because migrating it later would create two method and serialization contracts.
Constructors may use internal lightweight tables while building a result, but
all exported analysis functions return the validated S4 object.

The canonical result fields, evidence fields, `evidence_type` vocabulary,
`direction` vocabulary, and provenance-chain fields live only in the main
roadmap. This execution plan must consume those controlled vocabularies rather
than maintain a second list.

## Development phases

### Phase 0: schema and object design

Deliverables:

- implement the canonical S4 `doseInterpretResult` constructor and validity
  checks from the main roadmap
- implement, rather than redefine, the canonical `result` and `evidence` schemas
- implement `as.data.frame()`, `summary()`, `show()`, and `evidence()`
- write unit tests for object construction and evidence filtering

Acceptance criteria:

- a toy result can be created without external data
- evidence rows can be filtered by target, type, and source
- summary output is useful without LLM support
- no new required runtime dependency is introduced
- invalid organism, ontology, input, and ID combinations return actionable
  errors before enrichment is called
- `input = "phenotype"` and `target = "both"` return actionable "not yet
  enabled" errors until their scientific and data contracts are active

### Phase 1: gene-based disease interpretation

Deliverables:

- implement the first `interpretDisease()` path for human gene vectors
- reuse existing `enrichDisease()` and `get_dose_data()` computation only after
  validating species and ontology compatibility; do not infer species from a
  relabeled result slot
- convert enrichment hits into `result` rows
- convert driving genes into `evidence` rows
- support HDO disease targets only in this phase; existing HPO and NCG
  enrichment remain lower-level analyses and must not be relabeled as disease
  targets

Acceptance criteria:

- `interpretDisease(c("1234", "5678"), ontology = "HDO")` returns a structured result
- `as.data.frame()` returns the ranked disease table
- `evidence()` returns driving gene rows
- results are consistent with direct `enrichDisease()` output
- every returned `target_type` is a human disease or mouse model; ontology terms
  remain evidence features unless linked through a validated target profile
- Phase 1 does not expose phenotype input as a public working path and rejects
  `target = "both"` with an actionable error

### Phase 2: ranked genes and gene-set lists

Deliverables:

- support named numeric vectors through existing GSEA machinery
- support named lists of gene sets, such as clusters or cell types
- add `query_id` to distinguish multiple input sets
- return comparable result and evidence tables across query sets

Acceptance criteria:

- a named ranked vector uses `gseDisease()` internally
- a named list returns one combined result table with query labels
- no single-cell or spatial package is imported
- examples can use plain R vectors and lists

### Phase 3: template explanation

Deliverables:

- implement `explainDisease(method = "template")`
- produce short evidence-grounded explanations
- include strongest support, missing evidence, and ambiguity when present
- support English and Chinese templates if practical

Acceptance criteria:

- explanation contains only facts present in `result`, `evidence`, or `sources`
- output works offline
- tests verify that unsupported evidence is not mentioned

### Phase 4: cross-species evidence integration

Deliverables:

- connect this plan to the validated cross-species roadmap data layer
- add ortholog evidence rows
- add mouse model evidence rows
- add model phenotype evidence rows when available
- let `rankMouseModels()` and `inferHumanDisease()` return `doseInterpretResult`;
  explanation helpers consume these objects and do not rerun analysis
- depend on the main roadmap's frozen disease ID mapping, model entity,
  MP-HPO mapping, release manifest, and provenance contracts
- treat the canonical mouse-model target as a genotype- or allele-defined
  disease-model entity; phenotype profiles remain evidence linked to that target

Acceptance criteria:

- `rankMouseModels()` can expose evidence through `evidence()`
- `inferHumanDisease()` can expose ortholog and model support
- `interpretDisease(..., target = "both")` can combine disease and model evidence
- score components remain visible rather than collapsed into a black-box rank
- negated phenotypes cannot appear as positive evidence
- benchmark fixtures show that held-out labels and all evidence derived from the
  same provenance chain are absent from discovery features
- low phenotype-mapping coverage produces an unavailable component, not zero

### Phase 5: LLM explanation through `aisdk`

Deliverables:

- implement an optional provider adapter behind a small internal interface
- define a strict prompt template
- serialize validated structured evidence as JSON or a compact markdown table
- record provider, model, evidence hash, prompt hash, response hash, generation
  parameters, response metadata, and generation time
- provide a deterministic fake provider for offline contract tests

Prompt constraints:

```text
You may only use the evidence supplied below.
Do not add external facts.
Do not upgrade weak evidence to strong evidence.
State uncertainty when evidence is incomplete or ambiguous.
```

Acceptance criteria:

- package checks pass without `aisdk`
- offline tests validate request serialization, response parsing, metadata,
  provider errors, and unsupported-fact rejection without credentials
- optional live smoke tests are skipped when `aisdk` or credentials are absent;
  core adapter tests are never skipped
- generated claims can be traced to evidence IDs, and an unverifiable claim
  fails validation or is explicitly marked as generated and unsupported

## Relationship to the existing roadmap

The existing cross-species roadmap should remain the main scientific roadmap. It already defines:

- cross-species positioning
- MGI, HPO, MPO, HDO, MONDO, GenCC, Open Targets, and IMPC data sources
- canonical data products
- `rankMouseModels()`
- `inferHumanDisease()`
- `explainDiseaseModel()`
- `explainDiseaseSimilarity()`
- `geneDiseaseProfile()`
- benchmarks and publication plan

This document should be treated as the execution plan for the evidence-first and LLM-assisted interpretation layer.

Its phase numbers describe software work, not permission to bypass scientific
gates. Phase 0-3 can wrap existing HDO analysis while the main roadmap's Phase
0-1 data work proceeds. This plan's Phase 4 cannot begin until those data
contracts pass. Phase 5 depends on a stable evidence contract and is not part of
the first paper's core scientific claim.

Recommended linkage:

- keep `docs/cross-species-disease-interpretation-roadmap.md` as the main roadmap
- add a short reference from that roadmap to this plan
- avoid duplicating the cross-species data source details here

## Suggested implementation order

1. Add result and evidence schemas.
2. Add `evidence()` accessor.
3. Wrap existing `enrichDisease()` output into `interpretDisease()` for gene vectors.
4. Add named gene-set list support.
5. Add offline template explanation.
6. Integrate MGI and cross-species evidence only after the main roadmap's
   scientific gates pass.
7. Add the optional `aisdk` explanation adapter and offline provider tests.

This order keeps the first pull requests small. It also ensures the package gains value before any LLM dependency is introduced.

## Non-goals

- Do not build a clinical rare disease diagnosis tool.
- Do not let an LLM rank diseases or modify scores.
- Do not require API keys for core package functionality.
- Do not import single-cell, spatial, or workflow-specific packages into core DOSE.
- Do not force auxiliary evidence tables into `gson` when they are not gene sets.

## Deferred scope decisions

- Phase 1 exposes HDO disease targets. HPO enters `interpretDisease()` only
  through validated disease profiles after the phenotype gates pass. NCG stays
  in `enrichNCG()`/`gseNCG()` until a canonical disease mapping exists.
- When the adapter is implemented, `aisdk` is declared in `Suggests` and loaded
  with optional runtime detection; it is never a required dependency.
- English template explanations are required first. Chinese templates may be
  added after the evidence vocabulary and snapshot tests are stable; LLM output
  is not the only permitted route to future Chinese support.
