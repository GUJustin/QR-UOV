# QR-UOV content ledger

This ledger controls the aggressive consolidation of `QRUOV_main.tex`.
The original source and the mechanical split remain available for comparison.

## Central thesis

QR-UOV admits two complementary structural reductions. The direct route gets
stronger as the vinegar dimension grows, while the directional route gets
stronger as it shrinks. Their nonlinear core is therefore always bounded by
dimension `m - 4`. Under the stated symbolic solver model, the local separator
refinement places the Level I cost below its category target.

## Preservation rules

- Preserve the abstract unless a human author chooses to edit it.
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
| Introduction | Rewritten in Slice 4 | Original remains in the canonical source | The old section began with attack mechanics. The replacement builds the cryptographic context and prior attack landscape before stating the fixed-output thesis. |
| Related work | Keep in shorter form | Full current text remains recoverable from the original | The main paper needs positioning, not a catalog. |
| From the expanded key to universal forgery | Keep the model, graph identity, public map, certificate chain, and bounds used later | Keep the full pull-back and random-matrix calculations in the new certificate appendix | The attacks need the exact interface and bounds, but not every derivation. |
| Direct forgery in `m - 4` variables | Keep target normalization, the constructive isotropic-subspace reduction, the usable-root bound, the verified attack, and the fixed-output theorem | Keep the local distribution, dimension-count, Jacobian, second-moment, preprocessing, and toy-execution details in `Direct-forgery proofs and preprocessing` | This is one half of the central thesis. The main paper now follows the attack in execution order. |
| Graph direction and derivative kernel recovery | Keep the planted root, graph completion, and success statements | Move detailed density and condenser proofs | This is the other half of the central thesis. |
| Jacobian localization | Keep the local condition and resulting solver hypothesis | Move the commutative algebra proof | The solver needs the result, not the full proof in the main paper. |
| Finite characteristic root recovery | Keep the solver interface, cost, success probability, and verified recovery theorem | Move the product start, separator, lifting identities, proof crosswalk, and proofs | The main paper uses the solver as an attack component. |
| Dense quadratic batching and flattened lifting | Keep one cost consequence | Move the construction and proof | This supports the ledger but is not part of the structural attack. |
| Sparse working field models | Keep the chosen field degrees in the cost table | Move the irreducibility and multiplication proofs | These details belong with the full ledger. |
| Local separator headline | Keep | Move the full derivation as needed | This gives the candidate Level I result. |
| Characteristic polynomials without global separation | Keep the key lemma and its consequence | Move the full parametric degree proof | This is the main solver refinement. |
| Local separation | Keep the local failure bound and simple factor result | Move the valuation proof and counting details | The local condition explains the cost improvement. |
| Frobenius descent | Keep the exact descent lemma | Move only supporting detail if space requires it | This is short and removes a probabilistic error term. |
| Sharpened root bound | Keep the final bound in the security table | Move the derivation | The derivation is a ledger detail. |
| Finite ledger | Keep one combined table and a short sensitivity statement | Move component calculations and secondary tables | Readers need the result and its robustness, not every arithmetic line. |
| Reduced parameter implementation | Merge its claim boundary into `Evidence and limitations` | Keep the full original section, including the trace table and regression details | The implementation shows composition, not full parameter cost. |
| Validation and remaining risk | Merge into one limitations subsection | Keep the detailed audit record | The current section repeats claims made elsewhere. |
| Assessment | Remove as a standalone section | Preserve in the original source | It repeats the introduction, security section, and conclusion. |
| Combined security and parameter implications | Keep and merge with the scissors theorem and ledger | Move repeated attack descriptions | This is where the paper should state its design consequence. |
| Executable audits | Replace with one artifact paragraph | Move the inventory to the supplement or repository README | The current inventory reads as a release report. |
| Conclusion | Rewrite only after the new body is stable | No | The conclusion should state the exact structural result and the conditional cost result once. |

## Flow architecture after Slice 3

The main paper now has 11 top-level sections. Sections 1 through 7 establish
the model, the two structural routes, localization, and the generic root
solver. Section 8 turns those results into the concrete candidate estimate in
this order:

1. Candidate local-separator bound.
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
