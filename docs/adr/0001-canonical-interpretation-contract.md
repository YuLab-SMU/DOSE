# ADR 0001: Canonical interpretation contract

## Status

Accepted

## Context

DOSE is moving from enrichment and semantic similarity toward disease
interpretation and cross-species model prioritization. Two planning documents now
exist:

- `docs/cross-species-disease-interpretation-roadmap.md`
- `docs/deeprare-inspired-dose-development-plan.md`

This creates an immediate design risk: the execution plan can accidentally
introduce a second public signature, a second result schema, or a second
explanation contract. That would make implementation, testing, serialization,
and downstream methods drift apart before the first public release.

The package also already enforces organism and ontology compatibility before
analysis. Any new interpretation interface must preserve that discipline rather
than hiding cross-species conversion behind relabeled result metadata.

## Decision

DOSE uses one public interpretation contract.

- `interpretDisease()` is the single high-level public entry point.
- The normative public signature, compatibility matrix, and result schema live in
  `docs/cross-species-disease-interpretation-roadmap.md`.
- `docs/deeprare-inspired-dose-development-plan.md` is an execution plan and must
  not define a second signature or a second public result schema.
- All public analysis functions return the same S4 class,
  `doseInterpretResult`, from the first public release.
- `evidence()`, `summary()`, `show()`, and explanation helpers operate on that
  same result object rather than on parallel list or temporary public classes.
- Explanation is always downstream of computed evidence. Template and LLM
  explanations may summarize evidence, but they do not create facts, alter
  scores, or rerun analysis.
- Cross-species conversion remains explicit evidence. It is never represented as
  a silent organism relabeling step.

## Consequences

### Positive

- Public API design stays stable across disease, phenotype, and mouse-model
  workflows.
- Tests, methods, and serialization target one object contract instead of two.
- Documentation can distinguish "scientific roadmap" from "execution plan"
  without semantic drift.
- LLM support remains optional and subordinate to structured evidence.

### Negative

- The execution plan must avoid convenient shorthand that looks like a new
  signature.
- Early implementation work must commit to S4 sooner than a temporary list-based
  prototype would.
- New result fields or accessor arguments now require contract-level review
  rather than ad hoc doc updates.

## Follow-up

- Keep `CONTEXT.md` aligned with interpretation vocabulary as terms are refined.
- When the execution plan gives examples, treat them as usage examples only and
  keep them compatible with the canonical signature.
- If explanation traceability needs a claim-level schema later, record that as a
  new ADR instead of quietly extending this one.
