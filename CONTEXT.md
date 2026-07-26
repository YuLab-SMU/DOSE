# DOSE Domain Context

## Glossary

### Disease interpretation

Evidence-grounded ranking and explanation of human diseases or mouse models from
gene, ranked-gene, gene-set, or phenotype inputs. This is a research
interpretation workflow, not a clinical diagnosis workflow.

### Target

A ranked entity returned by an interpretation workflow. In the current domain,
a target is either a human disease or a genotype-/allele-defined mouse disease
model. Phenotype terms are evidence features, not targets.

### Result row

One ranked target for one query. A result row is the summary view consumed by
`as.data.frame()`, `summary()`, and plotting helpers.

### Evidence row

The atomic trace record that explains why a target appears in the result. An
evidence row links one query, one target, one evidence type, and one source
record. Evidence rows can represent support, missing, conflict, or ambiguity.

### Explanation

A human-readable summary generated from structured result and evidence rows.
Explanation is derivative output. It does not define scores, create evidence, or
change rankings.

### Direct annotation

An evidence path that uses annotations native to the declared input organism and
target space, without ortholog or phenotype mapping.

### Cross-species path

An evidence path that crosses between human and mouse through explicit ortholog
or phenotype mapping. The conversion itself is part of the evidence and must
remain visible to the user.

### Mouse model

A genotype- or allele-defined mouse disease-model entity used as a target in
model prioritization. It is not synonymous with a gene, phenotype term, disease
term, or phenotype profile.

### Ontology compatibility

The rule that only organism-compatible ontologies may be used as direct
annotation sources. Invalid organism and ontology combinations fail before
analysis. Cross-species use requires an explicit recorded conversion path.
