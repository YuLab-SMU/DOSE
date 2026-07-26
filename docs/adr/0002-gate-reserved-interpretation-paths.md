# ADR 0002: Gate reserved interpretation paths and model targets

## Status

Accepted

## Context

The canonical `interpretDisease()` signature already reserves phenotype input
and `target = "both"`, but the data and scientific contracts behind those paths
arrive later than the first public implementation. The planning documents also
need one stable answer for what a mouse-model target actually is.

## Decision

- `input = "phenotype"` stays in the canonical signature, but it is not a
  working public path until the phenotype mapping, coverage, and provenance
  gates pass. Early calls fail with an actionable "not yet enabled" error.
- `target = "both"` stays in the canonical signature, but it is not available
  before the Phase 4 cross-species target contract is implemented. Early calls
  fail with an actionable "not yet enabled" error instead of silently degrading
  to disease-only output.
- A canonical mouse-model target is a genotype- or allele-defined disease-model
  entity. Phenotype profiles are evidence attached to that target, not
  alternative targets.

## Consequences

This keeps one stable public contract without pretending unfinished paths are
already valid. It also prevents later drift between genotype-level model
ranking and profile-level evidence summaries.
