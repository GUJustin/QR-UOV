# Revision and correctness notes

## Strict full-break status

No end-to-end attack in this bundle completes on an official QR-UOV parameter set. The paper therefore makes no practical full-break claim and contains no claim about a looming or expected breakthrough. Only proved reductions, reproducible symbolic cost diagnostics, and executable toy-model checks are reported.

This revision strengthens the paper in two exact ways:

1. The target-orthogonal direct-forgery residual is solved projectively, removing one homogeneous radial variable.
2. Every auxiliary-field multiplication is charged through an executable recursive schoolbook/Karatsuba circuit rather than a schoolbook convolution envelope.

The resulting corrected Boolean-gate diagnostics are:

| Level | Projective direct forgery | Directional equivalent-key recovery | Best row | Target | Best margin below target |
|---|---:|---:|---:|---:|---:|
| I | 151.639 | 154.826 | 151.639 | 143 | -8.639 |
| III | 203.507 | 206.177 | 203.507 | 207 | 3.493 |
| V | 260.389 | 260.857 | 260.389 | 272 | 11.611 |

A positive final-column value means that the modeled attack row lies below the category target. The level-specific conclusion is:

- Level I is not conservatively broken by these rows.
- Level III is below target through two independent reductions: the projective direct row by 3.493 bits and the equivalent-key row by 0.823 bits.
- Level V is below target through both reductions by more than 11 bits.

These values are reproducible symbolic-algorithm diagnostics, not measured implementation costs. The direct symbolic representations remain physically unrealistic in memory.

## Main strengthening 1: exact projective direct forgery

### Affine residual before this revision

After target normalization and target-orthogonal slicing with three forms, the earlier direct attack constructed an `r = m - 3` dimensional subspace `U` and solved

```text
H_1(z) = ... = H_{r-1}(z) = 0,
H_0(z) = 1
```

in `r` affine variables. This gave residual dimensions 51, 75, and 102. Those affine dimensions were already present in Hashimoto-based analysis.

### Projective quotient

The first `r - 1` equations are homogeneous. Over an odd field, affine roots modulo `z ~ -z` are in bijection with projective points `[z]` satisfying

```text
H_1(z) = ... = H_{r-1}(z) = 0
```

for which `H_0(z)` is a nonzero square. The last equation is therefore used only after root extraction: square testing and one field inversion recover the affine scale.

The affine Jacobian is invertible exactly when the corresponding projective zero is nonsingular. Euler's identity supplies the equivalence: the gradients of the zero-target forms annihilate the radial vector, while the gradient of `H_0` evaluates to 2 on it.

A uniformly random nonzero chart form `L` contains any fixed projective root with exact probability

```text
1 - (q^(r-1) - 1)/(q^r - 1) > 1 - 1/q.
```

Thus the square solver dimensions become:

```text
50, 74, 101.
```

The existing QR-UOV-specific second-moment theorem already proves the existence of an eligible nonsingular affine root with probability at least

```text
0.326336998767509...
```

and hence proves the same projective existence event. The extra random-chart normalization costs less than 0.011405 bits at `q = 127`.

### Consequences for the working fields

The projective direct solver uses geometric degrees `2^50`, `2^74`, and `2^101`. Its audited working-field degrees are now:

```text
15, 22, 30.
```

The Level-I field degree drops from 16 to 15. The new degree-15 sparse modulus is

```text
X^15 + X^12 - X + 1,
```

and the companion audit certifies its irreducibility over `F_127` with Rabin's criterion.

## Main strengthening 2: executable Karatsuba field circuits

For a length-`d` coefficient convolution, the audit compares schoolbook multiplication with a recursively realized Karatsuba plan. With

```text
a = ceil(d/2), b = floor(d/2),
```

the recursive candidate has counts

```text
M_d = 2 M_a + M_b,
A_d = 2 A_a + A_b + 4d - 4.
```

The `4d - 4` additions cover forming the two block sums, subtracting the low and high products from the middle product, and combining shifted outputs. The implementation executes the same uneven-split recursion and checks its products against schoolbook convolution.

Sparse reduction uses only additions and subtractions because all nonleading modulus coefficients are `+1` or `-1`. The selected official-field plans are:

| Degree | Base-field multiplications | Base-field additions | Nonleading modulus terms |
|---:|---:|---:|---:|
| 15 | 81 | 332 | 3 |
| 18 | 153 | 404 | 2 |
| 22 | 225 | 480 | 2 |
| 24 | 243 | 512 | 3 |
| 30 | 243 | 1112 | 2 |

Relative to the previous sparse schoolbook envelope, this saves about 0.604 to 1.037 bits on the directional rows and 0.674 to 1.037 bits on the direct rows.

## Sharper fixed-output security scissors

Let the quotient degree be `ell`, the output count be `m`, and the vinegar count be `V`. The two attack families now give:

- directional equivalent-key core `V` when `V <= m - ell - 1`;
- projective target-orthogonal core `m - ell - 1` when `V >= m - ell`.

Because these cover every integer `V`, for every nonzero target over an odd field,

```text
r_attack(V) <= m - ell - 1.
```

For degree-three QR-UOV, the fixed-output ceiling is therefore

```text
m - 4,
```

not `m - 3`. Increasing only the vinegar count cannot move the combined algebraic bottleneck past this value. A repair must increase `m` as well or alter the construction so that one of the reductions no longer applies.

This is an algebraic-core theorem. It does not claim that runtime is determined only by the core dimension or that every QR-UOV-like parameterization is impossible.

## Revised parameter implications

The companion audit searches quotient-compatible `(V,m)` pairs using the same padded ledger for both attack families.

| Level | Current `(V,m)` | First pair reaching target | First pair reaching target plus 16 bits |
|---|---:|---:|---:|
| I | (52,54) | (52,54) | (57,60) |
| III | (76,78) | (78,81) | (87,90) |
| V | (102,105) | (108,111) | (117,120) |

The table is an attack-specific lower-bound diagnostic, not a complete replacement proposal. Other attacks, key sizes, signatures, failure probabilities, and implementation costs must be recomputed.

## Correctness checks and corrections

This pass also made the following precision and presentation changes:

1. The projective theorem and the scissors theorem now state the required odd-characteristic hypothesis explicitly.
2. The paper no longer describes the projective step with an unqualified novelty slogan; it distinguishes the prior affine residual dimensions from the additional projective quotient analyzed here.
3. PDF metadata now contains the full paper title.
4. The claim-status table now describes the exact eligible-projective-root event rather than using the older affine-root label.
5. The directional audit output is labeled according to its actual purpose.

The earlier major corrections remain in force: reconstruction precision is `2*C0`, not `C0`; sparse auxiliary fields are certified; every returned candidate is checked against the unchanged public map; and density claims for seeded SHAKE/AES expansion are confined to the stated ideal-primitive models with explicit transfer terms.

## Executable validation

### Combined audit

`QRUOV_scissors_audit.py` recomputes:

- Chevalley-Warning dimension slacks;
- projective cores 50, 74, and 101;
- fixed-key seed-oil intersection bounds;
- the exact nonsingular-root lower bound;
- chart and separator probabilities;
- recursive Karatsuba plans;
- direct and directional cost rows;
- memory envelopes;
- retuning thresholds; and
- Rabin certificates for every sparse modulus used in the official and threshold rows.

### Directional audit

`QRUOV_directional_audit.py` verifies:

- 28 Vandermonde start roots;
- 120 characteristic-127 lifted paths;
- 360 flattened bivariate products;
- 48 recursive Karatsuba convolutions against schoolbook multiplication; and
- 200 complete toy equivalent-key instances.

Of the 200 toy instances, 199 had a full-rank planted direction and all 199 recovered the exact graph. The remaining direction was correctly rejected for retry. No nonplanted root produced a verified key.

### Exact degree-three direct-forgery toy

`QRUOV_direct_qr_toy.py` executes the actual projective reduction over a small degree-three QR-UOV model. In the frozen 200-key run:

- 199 keys yielded a verified forgery within twelve attempts;
- every returned candidate passed independent base-field and extension-field evaluation;
- all selected seed points were off oil; and
- 185 of the 200 final seed spaces were disjoint from the hidden oil space.

The toy experiment validates attack logic, not full-parameter performance.

## Build and visual QA

- The revised paper has 65 pages.
- The source changes 648 lines and deletes 242 lines relative to `QRUOV_optimal.tex`.
- There are 239 unique labels and no missing or duplicate references.
- All 34 bibliography entries are cited, with no missing citations.
- LaTeX reports no package, citation, reference, overfull-box, or underfull-box warnings.
- All fonts are embedded and subset.
- All 65 pages were rendered at 160 dpi and visually inspected.
- The previous and revised PDFs were also rendered and compared page by page.

## Remaining limitation

The projective quotient and Karatsuba circuit materially improve the modeled work factors, but they do not solve the memory problem. One terminal flattened coefficient array still requires approximately

```text
2^49.41, 2^98.51, 2^153.40 EiB
```

for the three direct rows. The corresponding directional arrays are even larger. The paper therefore does not claim a practical full-parameter implementation, and it does not speculate about unimplemented attack variants.
