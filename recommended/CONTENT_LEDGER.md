# QR-UOV content ledger

This ledger controls the aggressive consolidation of `QRUOV_main.tex`.
The original source and the mechanical split remain available for comparison.

## Central thesis

QR-UOV admits two complementary structural reductions. The direct route gets
stronger as the vinegar dimension grows, while the directional route gets
stronger as it shrinks. One attack therefore always reduces to a square
nonlinear system in at most `m - 4` variables. Under the stated symbolic
solver model, target separation places the Level I cost below its category
target.

## Preservation rules

- The abstract was protected until explicit author approval. Slice 9 revises
  it after the attack spine and claim boundaries stabilized.
- Preserve every distinct idea and claim from the introduction. Its structure
  may change after the body has stabilized.
- Keep every exact algebraic reduction needed for the central thesis in the
  main paper.
- Keep the assumptions beside every conditional cost claim.
- Move unique proof material to the supplement instead of deleting it.
- Keep enough implementation detail to show the public input boundary and the
  final verifier check.
- Remove repeated summaries, revision history, audit narration, and repeated
  statements of the same limitation.
- Record any deleted result below with a reason and a pointer to the stronger
  result that replaces it.

## Main body decisions

| Current section | Main paper | Supplement | Reason |
| --- | --- | --- | --- |
| Introduction | Rewritten in Slice 4 and tightened in Slice 9 | Original remains in the canonical source | The introduction now builds the scheme and attack context, asks the parameter question, and states the fixed-output theorem in descriptive language before introducing technical details. |
| Related work | Keep in shorter form | Full current text remains recoverable from the original | The main paper needs positioning, not a catalog. |
| From a recovered graph to signatures on any message | Keep the model, graph identity, public map, certificate chain, and bounds used later | Keep the full pull-back and random-matrix calculations in the new certificate appendix | The attacks need the exact interface and bounds, but not every derivation. |
| Direct forgery in `m - 4` variables | Keep target normalization, the constructive isotropic-subspace reduction, the usable-root bound, the verified attack, and the fixed-output theorem | Keep the local distribution, dimension-count, Jacobian, second-moment, preprocessing, and toy-execution details in `Direct-forgery proofs and preprocessing` | This is one half of the central thesis. The main paper now follows the attack in execution order. |
| Equivalent-key recovery from a planted root | Keep the public directional system, planted-root identity, exact graph completion, rank-density statement, and derivative-based oil-space recovery | Keep the multi-direction fallback, Wronskian proof, deterministic lists, exact list counts, and optimality argument in `Deterministic fallbacks for equivalent-key recovery` | This is the other half of the central thesis. The main paper now exposes the exact certificate chain from one root to a verified signing key. |
| Why one root with invertible Jacobian is enough | Keep the localization polynomial, its planted-root value, and the resulting solver hypotheses | Keep the commutative-algebra proof in `Why the local Jacobian condition suffices` | The solver needs the exact local implication, not its proof in the main flow. |
| Finding roots over the QR-UOV base field | Keep the solver target, the two characteristic-127 replacements, cost, success probability, and verified recovery theorem | Keep the product start, separator, lifting identities, and proof crosswalk in the finite-characteristic appendix | The main paper uses the solver as an attack component. |
| Dense quadratic batching and flattened lifting | Keep the exact cost consequence and corrected terminal precision | Keep the construction and proof in `Concrete solver derivations and secondary ledgers` | This supports the ledger but is not part of the structural attack. |
| Sparse working field models | Keep the selected field degrees and one cost envelope | Keep irreducibility certificates and multiplication proofs in `Concrete solver derivations and secondary ledgers` | These details belong with the full ledger. |
| Target-separation headline | Keep | Keep the full derivation in `Concrete solver derivations and secondary ledgers` | This gives the candidate Level I result. |
| Collisions among rejected roots do not raise the degree | Keep the key lemma and its consequence | Move the full parametric degree proof | This is the main solver refinement. |
| Local separation | Keep the local failure bound and simple factor result | Move the valuation proof and counting details | The local condition explains the cost improvement. |
| Base-field power test | Keep the exact field-membership lemma | Move only supporting detail if space requires it | This is short and removes a probabilistic error term. |
| Sharpened root bound | Keep the final bound in the security table | Move the derivation | The derivation is a ledger detail. |
| Finite ledger | Keep the final Level-I total, principal components, and a short sensitivity statement | Keep component calculations and secondary tables in `Concrete solver derivations and secondary ledgers` | Readers need the result and its robustness, not every arithmetic line. |
| Reduced parameter implementation | Merge its claim boundary into `Evidence and limitations` | Keep the full original section, including the trace table and regression details | The implementation shows composition, not full parameter cost. |
| Validation and remaining risk | Merge into one limitations subsection | Keep the detailed audit record | The current section repeats claims made elsewhere. |
| Assessment | Remove as a standalone section | Preserve in the original source | It repeats the introduction, security section, and conclusion. |
| Combined security and parameter implications | Keep the common cost table, memory boundary, and fixed-output design consequence | Move repeated attack descriptions | This is where the paper should state its design consequence. |
| Executable audits | Replace with one artifact paragraph | Move the inventory to the supplement or repository README | The current inventory reads as a release report. |
| Conclusion | Rewrite only after the new body is stable | No | The conclusion should state the exact structural result and the conditional cost result once. |

## Flow architecture after Slice 3

The main paper now has 11 top-level sections. Sections 1 through 7 establish
the model, the two structural routes, localization, and the generic root
solver. Section 8 turns those results into the concrete candidate estimate in
this order:

1. Candidate target-separation bound.
2. Dense-quadratic batching and flattened lifting.
3. Sparse working-field models.
4. Characteristic polynomials without global separation.
5. Local separation at the target point.
6. Exact one-sketch Frobenius descent.
7. Root success probability.
8. Finite Level-I ledger.

Sections 9 through 11 then separate evidence and limitations, parameter
implications, and the conclusion. Slice 3 changed hierarchy and order only.
It did not delete or move technical content to the supplement.

## Existing appendices

| Appendix | Action |
| --- | --- |
| Independent projective smooth core route | Keep in the supplement. |
| Finite field solver interfaces and verified restart | Keep in the supplement. |
| Expanded key certificate calculations | Keep in the supplement. |
| Characteristic-127 solver details | Keep the finite-field start, random separator, Newton--Hensel identities, renormalization, and proof crosswalk. |
| Concrete solver derivations and secondary ledgers | Keep the complete batching, sparse-field, characteristic-polynomial, target-separation, Frobenius, root-probability, and ledger derivations. |
| Deterministic fallbacks for equivalent-key recovery | Keep in the supplement. |
| Direct-forgery proofs and preprocessing | Keep in the supplement. |
| Official seeded expansion | Keep in the supplement. |
| Optional identification with the planted oil space | Keep in the supplement. |
| Alternative elimination and low space directions | Keep in the supplement. |
| Reproducibility and executable attack description | Keep in the supplement, then merge with moved audit material. |

## Deletion log

No mathematical result has been deleted.

- Slice 1 replaced the standalone main-paper sections `Validation and
  remaining risk`, `Assessment`, and `Executable audits` with `Evidence and
  limitations`.
- The exact algebraic claims, toy-scale evidence, and four conditions on the
  full-parameter estimate remain in the main paper.
- The complete removed text, including every audit name and numerical result,
  is preserved in Supplement G, `Archived validation, assessment, and audit
  inventory`.
- The self-assessment and manuscript-history sentences no longer carry the
  main argument. The sensitivity tables and `Evidence and limitations` carry
  their substantive claims.
- Slice 2 moved the complete `Reduced-parameter end-to-end instantiation`
  section to Supplement F. The main paper now keeps only the public input and
  verifier boundary, the experiment boundary, the successful run count, and
  the warning that the prototype does not implement the full-scale schedule.
- The moved prose and table data match the canonical source. Only the section
  label and the table placement option changed. No experiment detail or table
  entry was deleted.
- Slice 2 corrected the short run count. Seed `20260731` and regression seeds
  `20260000` through `20260004` give six distinct successful $q=13$ runs, not
  five. The separate $q=29$ run also succeeded.
- Slice 3 grouped the eight concrete solver and cost sections under one
  top-level section and moved the existing candidate theorem to its first
  subsection. No mathematical result or supporting detail was deleted.
- Slice 4 replaced the introduction's proof tour with a cited context, question,
  thesis, and four-part contribution overview. The detailed claim-status table
  moved to `Evidence and limitations`; every table row remains in the main
  paper. The original introduction remains in `QRUOV_strict_final.tex`.
- Slice 5 replaced the 1,977-word expanded-key section with a 1,054-word
  certificate chain. The main paper still defines both expansion models,
  states every public check, and retains all signing-search bounds used later.
  The complete pull-back calculation, rank distribution, pair-line count, and
  determinant argument now appear in `Expanded-key certificate calculations`
  in the supplement. No mathematical result or numerical bound was deleted.
- Slice 6 replaced the 3,202-word direct-forgery section with a shorter attack
  narrative. The main paper still contains target normalization, the exact
  fixed-key seed bound, the constructive Chevalley--Warning theorem, the
  `m - 4` projective reduction, the explicit root probability, the verified
  attack, and the fixed-output theorem. The complete local-randomness proof,
  extension count, Jacobian argument, second-moment calculation, preprocessing
  ledger, and degree-three toy execution now appear in `Direct-forgery proofs
  and preprocessing`. No mathematical result or numerical bound was deleted.
- Slice 7 replaced the 1,616-word directional-recovery section with a 724-word
  certificate chain. The main paper still defines the public directional
  system, proves that `Gc` is its planted root, gives the exact linear systems
  that recover every column of `G`, states the rank probability used by the
  verified theorem, and retains the derivative-kernel route. The complete
  multi-direction fallback, Wronskian construction, deterministic list bounds,
  numerical counts, and optimality argument now appear in `Deterministic
  fallbacks for equivalent-key recovery`. No mathematical result or numerical
  bound was deleted.
- Slice 8 replaced the main-paper proof of the localization lemma with a short
  proof map. The main paper still defines the localization polynomial, states
  the exact regular-sequence and radical-prefix conclusion, evaluates the
  determinant at the planted root, and identifies the two hypotheses of the
  locally closed solver. The complete Jacobian-criterion and
  Cohen--Macaulay proof now begins the finite-field solver appendix. No
  mathematical statement was deleted.
- Slice 9 revised the abstract and introduction after explicit author
  approval. The abstract still explains UOV, QR-UOV's quotient-field
  compression, the fixed-output question, both reductions, the three concrete
  estimates, the reduced-parameter experiment, and every major limitation.
  The introduction retains the NIST context, UOV history, structured-key
  lineage, attack landscape, question, theorem, two contributions, concrete
  rows, evidence boundary, and organization. It now states the theorem as a
  bound on variables in a square nonlinear system instead of leading with the
  manuscript's internal term `solver core`. No citation, claim, numerical
  result, or limitation was deleted.
- Slice 10 reduced the finite-characteristic solver to its attack-facing
  contract and moved 1,062 words of construction and proof detail to the
  supplement. It also reduced the concrete-cost section to the refinements
  that determine the final ledger and moved 2,909 words of derivations and
  secondary tables to `Concrete solver derivations and secondary ledgers`.
  Every theorem statement, exact probability, selected field, numerical cost,
  sensitivity boundary, and unresolved assumption used by the main claim
  remains in the main paper. No mathematical result or numerical row was
  deleted.
- Slice 13 audited the complete main paper in reading order. It added plain
  definitions where the compressed draft assumed solver-specific vocabulary,
  repaired transitions between the scheme, attacks, solver, evidence, and
  security synthesis, and removed reader-facing shorthand such as `solver
  core`, `directional recovery`, `global separation`, and `Frobenius descent`
  where a descriptive phrase suffices. It also corrected two symbols used
  before definition, one reused symbol, one duplicated phrase, and an
  overstatement about which attack the implementation executes. No formula,
  theorem, numerical result, citation, or claim limitation was removed.
