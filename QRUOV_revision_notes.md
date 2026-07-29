# QR-UOV paper: revision and audit notes

## Bottom line

The paper's central structural attack survived the audit and is now substantially stronger and cleaner. The revised paper proves that one regular public direction determines the entire hidden oil space, gives compact deterministic direction and signing lists, and wraps the polynomial-system solvers in exact verification and finite-invocation restart logic.

This pass also found several substantive problems in the earlier draft. The most important were a mismatch between the fast solver's exponent parameter and its required field size, an unjustified internal-restart step, an invalid conditional-uniformity statement in the seeded/distributional analysis, and the transfer of explicit failure constants to an older theorem that did not state them. These have all been removed or repaired.

The final paper is 48 pages. It compiles without unresolved references, missing citations, duplicate labels, overfull boxes, or underfull boxes. The companion numerical audit uses only the Python standard library and passes frozen assertions for the headline values.

## Stronger results added

### 1. One regular direction determines the whole hidden graph

After normalizing the hidden oil space to

```text
O = {(-Gc,c) : c in K^O},
```

the public direction system `P_i(-x,c)=0` has planted root `x=Gc`. Translation by that root gives the exact identity

```text
P_i(-(Gc+z),c) = z^T A_i z - 2(B_i c)^T z.
```

If the planted derivative matrix has full column rank, one recovered root `Gc` determines every other column of `G` through public linear systems. The revised theorem is key by key: it needs only local nonsingularity at the planted root, not uniqueness of the root, a global smoothness assumption, or the projective discriminant analysis.

### 2. Optimal Wronskian condenser and much shorter deterministic lists

The previous square reduction contemplated enumerating row subsets. The revision introduces an explicit one-parameter Wronskian condenser. For every full-rank `m x V` matrix, a public list of exactly

```text
V(m-V)+1
```

condenser values contains an invertible `V x V` core. A matching lower bound shows that this list length is optimal for the stated one-parameter family.

At Levels I, III, and V, the condenser lists have lengths 105, 153, and 307. In the public-random expanded-key model, testing only the first 3, 4, and 4 coordinate directions therefore yields compact deterministic lists of

```text
315, 612, 1228
```

square systems. Their rank-exception logarithms are at most

```text
-188.694, -251.593, -335.457.
```

Scanning all coordinate directions gives the stronger exact exceptional-set logarithms

```text
-1132.167, -1635.352, -2935.248.
```

### 3. Basis-first signing certificate

The signing search is now both simpler and stronger. Fixed basis vinegar assignments produce independent uniform square oil matrices over the base field. As a result:

- the expected basis scan uses at most `1.0079995` tests;
- testing only 22, 31, and 40 fixed assignments gives failure logarithms below `-153.502`, `-216.298`, and `-279.095`;
- testing all basis assignments gives `-362.823`, `-530.280`, and `-711.692`;
- the designated pair-line fallback gives `-544.821`, `-796.277`, and `-1068.687`.

Combining the compact direction list with the short basis signing list gives a deterministic public certificate list whose total failure logarithms remain below `-153.502`, `-216.298`, and `-279.095`, without assuming independence between the recovery and signing failures.

### 4. Rigorous finite-invocation fast solver theorem

The fast projective route now links the solver exponent `epsilon` to the field-size condition required by the cited theorem. For the paper's primary choice `epsilon=0.3`, the minimum extension degrees are

```text
s = 8, 12, 15,
```

and one finite invocation succeeds with probability at least

```text
0.992337625602,
0.999997944684,
0.927725586161.
```

Using one additional extension degree, `s=9,13,16`, raises these lower bounds to

```text
0.999999996259,
0.999999999998,
0.999999964716.
```

Every solver call has a finite operation schedule. Exceptional divisions, degree mismatches, incomplete representations, and failed consistency checks return `fail`; they do not trigger an unbounded hidden loop. Exact public verification prevents false acceptance.

### 5. Cleaner independent cost cross-checks

The paper now keeps three solver interfaces distinct:

1. the modern nearly quadratic finite-field theorem;
2. the older straight-line-program arithmetic bound, used only as an epsilon-free growth diagnostic;
3. the locally closed-set theorem, specialized independently to both the projective and graph-local systems.

The displayed rows are explicitly described as normalized soft-O indicators, not measured gates or implementation claims. The conclusions are strongest at Levels III and V. Level I remains structurally important but more sensitive to hidden constants and representation choices.

## Definite errors and proof-boundary problems corrected

1. **Reversed affine-chart statement.** The earlier draft said to choose an affine chart avoiding all projective roots. The correct construction chooses the hyperplane at infinity to avoid the roots, so that the complementary affine chart contains them.

2. **Exponent/field mismatch in the fast solver.** The required extension field depends on the chosen exponent parameter. Treating the field list as independent of that parameter made the stated theorem interface incorrect. The field condition and the cost exponent are now linked throughout.

3. **Unjustified internal restart.** A Las Vegas guarantee under regularity hypotheses cannot be invoked as an unrestricted inner loop on an irregular outer draw. The revised wrapper makes one checked finite invocation and restarts only at the outer level.

4. **Unsupported failure constants on the classical SLP line.** Exact finite-invocation constants from the modern theorem were being attached to an older arithmetic theorem that did not state them. The SLP calculation is now clearly labeled an arithmetic cross-check and carries no borrowed probability claim.

5. **False conditioning on the secret graph.** One distributional proposition conditioned on `G` while the model explicitly allowed correlation between `G` and the symmetric blocks. That conditional uniformity need not hold. The revised proposition uses the correct marginal law on a fixed transverse slice; the attack and density proof do not need conditional uniformity given `G`.

6. **Norm-selector probability overstated.** A selector failure probability was stated as exactly `1/|L|`; the correct general statement is at most `1/|L|`.

7. **Lower-bound rounding errors.** Two displayed success lower bounds were rounded upward in their final digit. They are now conservatively rounded downward.

8. **Stale memory accounting.** The memory section used field degrees and numerical values from an older solver calculation. It now uses the corrected minimum fields. Even one dense degree-`2^V` vector requires about `0.0819 EiB`, `2.06 x 10^6 EiB`, and `1.73 x 10^14 EiB`; storing roughly `V` coordinate arrays is still larger.

9. **Nonexistent reproducibility claims.** The earlier text named a large repository of checks and a toy experiment that were not present in the uploaded material. Those claims were removed. The paper now describes one actual companion arithmetic audit and gives executable pseudocode for both attack routes.

10. **Signing-list and certificate ambiguity.** The earlier wording blurred a chosen pair, an invertible pair, and the deterministic fallback. The revised signing theorem separates the basis-first event, the short fixed prefix, the all-basis event, and the pair-line fallback.

11. **Overstated independence between solver routes.** The solver results are separate theorem interfaces and specializations in the Kronecker lineage, not wholly unrelated algorithmic families. The related-work and cost sections now say this directly.

## Readability and organization changes

- Rewrote the title, abstract, introduction, related work, cost interpretation, and conclusion.
- Added UOV and equivalent-key background before the algebraic machinery.
- Put the core graph-direction identity near the start and gave the attack as four concrete steps.
- Added a claim-status table distinguishing exact algebra, key-density statements, asymptotic solver theorems, normalized exponents, and memory costs.
- Added a terminology paragraph for regular roots, projective slices, geometric/Kronecker resolutions, key-by-key statements, and soft-O notation.
- Moved lower-priority implicit extraction and underslicing investigations to appendices.
- Replaced repeated defensive qualifications with one focused section explaining what the cost tables establish and what they do not.
- Added a compact executable attack description that makes all restart and verification boundaries explicit.
- Compressed and balanced the bibliography so the paper ends cleanly on page 48 rather than spilling one reference onto an almost empty final page.

## Verification performed

- Clean `latexmk` build of the final 48-page PDF.
- 169 labels, all unique; no unresolved cross-references.
- 29 bibliography entries, all cited; no missing citations.
- No LaTeX, package, overfull-box, or underfull-box warnings in the final log.
- PDF preflight: openable, unencrypted, not scanned, 48 consistent A4 pages, embedded/subset fonts, 45 outline entries.
- Rendered every page and visually inspected contact sheets plus the title page, introduction table, main recovery theorems, cost tables, memory section, reproducibility appendix, and bibliography.
- Independent high-precision recomputation of the direction-rank, signing, field-size, finite-invocation, cost-indicator, bad-key, seeded-transfer, and memory values.
- Frozen assertions in `QRUOV_numerical_audit.py` pass.
- Primary-source audit of the regularity/radicality assumptions, finite-field solver interface, locally closed-set interface, and complete-intersection discriminant reference.

## Remaining limitations

1. **No production solver implementation or benchmark.** The solver costs are theorem-backed asymptotic bit-complexity indicators, not measured gates or wall-clock estimates.
2. **The direct representations require astronomical memory.** The paper identifies smaller possible outputs, but no proved low-space construction of the required quotient algebra is known here.
3. **Fixed SHAKE/AES density remains a separate issue.** The key-by-key theorem is unconditional for every key satisfying its algebraic conditions. High-probability seeded-key statements are proved in the stated ideal-XOF and ideal-cipher models, not for the fixed public primitives by generic statistical replacement.
4. **No replacement parameter set is proposed.** The revised paper establishes a strong reason for a concrete parameter review, especially at Levels III and V; choosing replacements requires implementation-level solver accounting and comparison with all other attack families.
5. **This is an intensive mathematical audit, not machine-checked formal verification.** The exact algebra and numerical ledger were checked independently, but the full paper has not been formalized in a proof assistant.
