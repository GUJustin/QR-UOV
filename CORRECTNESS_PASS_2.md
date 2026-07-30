# Second adversarial correctness pass

## Verdict

I did not find a fatal contradiction in the local-separator mechanism. The central determinant/polar identity, local `C-1` bad-condition count, and simple-factor Frobenius descent remain intact. I did find two genuine defects in the presentation/accounting, both repaired in this bundle:

- the specialization proof claimed too much about escaping-branch numerator terms; and
- the ledger omitted one small leaf-formation cost during the final coordinate rerun.

The corrected Level-I primary total is `140.383995063`, rather than `140.383963623`. The change is `0.00003144` bits. The attack remains `2.616004937` bits below the inherited `143`-bit target.

## Correctness repairs

### 1. Escaping numerator terms do not necessarily vanish

For an escaping branch, common valuation normalization can leave a nonzero numerator contribution. The correct reason it is harmless is algebraic divisibility: that contribution contains every finite target factor, so it is a multiple of the finite specialized eliminant. It therefore vanishes modulo that eliminant and at each finite simple root.

Both the local upgrade and the preceding split-before-lift note now state and prove the corrected result. `toy_infinite_branch_specialization.py` gives an exact example where the escaping term survives but is divisible by the target eliminant.

### 2. One-time coordinate separator leaves

The successful coordinate rerun reconstructs all base-field branch series and evaluates a coordinate product tree at the accepted separator value. It must also reform `u_s=lambda.x_s`. The original local ledger charged the branch rerun and all long products, but not these pointwise extension-field operations.

The new explicit term is:

- Level I: `2^124.898162...` gates;
- Level III: `2^175.158104...` gates;
- Level V: `2^230.587052...` gates.

It is far below the dominant terms but is now included rather than hidden in a miscellaneous cushion.

### 3. Finite-field primitive-element circularity

The proof no longer says merely “extend to an infinite field.” It observes that, over the algebraic closure of the generic fiber, every pairwise difference `Lambda.(x_s-x_r)` is a nonzero polynomial in independent `Lambda` variables. Thus the generic linear element is primitive over `K'(t)(Lambda)` without requiring any `K'`-rational global separator. The determinant coefficients already lie over the original field, so the parameter-degree bound descends.

### 4. Sharpened root probability

The `0.3690869822...` probability is now derived explicitly. If `M` counts all affine normalized roots and `N` the nonsingular ones, the parent manuscript gives

`E[M]=mu`, `E[M^2]=2mu+mu^2-(q-1)mu q^(-r)`, and `E[N]=mu*pi`.

Pairing `z` with `-z` and applying the second Bonferroni inequality gives

`p_root >= mu*pi/2 - mu^2/8 + (q-1)mu q^(-r)/8`.

The expression is increasing on the relevant interval, so `mu=1-1/q` is a valid substitution. As a conservative cross-check, replacing this sharpening by the parent's established `0.326336998767509...` lower bound gives Level-I total `140.531460701`, still `2.468539299` bits below target.

### 5. Forbidden offset sign

For `psi=H0/(c+lambda.x)`, the forbidden values of `c` are `-lambda(P)`. Earlier prose called them “the lambda values.” The sign does not change any cardinality or probability calculation, but the statements are now exact.

## Updated totals

| Level | Primary total | Margin | Parent-root fallback | Fallback margin |
|---|---:|---:|---:|---:|
| I | 140.383995063 | 2.616004937 | 140.531460701 | 2.468539299 |
| III | 191.027133323 | 15.972866677 | 191.169964592 | 15.830035408 |
| V | 247.125438947 | 24.874561053 | 247.265866395 | 24.734133605 |

The independent no-import verifier reproduces every primary component to below `5e-10` bits.

## Sensitivity that matters most

The dominant unresolved finite-constant issue is online multiplication. Reoptimizing the extension degree at Level I gives:

| Multiplier on the whole online envelope | Best degree | Total | Margin |
|---:|---:|---:|---:|
| 1 | 8 | 140.383995 | 2.616005 |
| 2 | 8 | 141.117483 | 1.882517 |
| 4 | 9 | 141.951549 | 1.048451 |
| 8 | 9 | 142.856481 | 0.143519 |
| 16 | 9 | 143.806490 | -0.806490 |

This is reassuring but not unlimited robustness. A fully explicit online schedule is still required before the Level-I claim is publication-grade.

Rational reconstruction is less fragile: the primary Level-I total remains `141.8146` when its multiplier is increased from `8` to `128`, and remains below target through multiplier `256`.

## Remaining correctness blockers, ranked

1. **Exact theorem matching for the saturated dominant algebra.** Check the parent's ideal, saturation, generic-fiber degree, and normalization against Schost's characteristic-polynomial degree theorem line by line.
2. **Rational-filter graph.** Reprove that graph closure and saturation by `c+lambda.X` have degree at most `2C0` and introduce no additional parameter pole.
3. **Explicit online schedule.** Convert the dyadic envelope into a finite algorithm whose product count is literally the one charged, or implement a reduced-parameter prototype and certify the count.
4. **Rational reconstruction and factorization constants.** The sensitivity is strong, but the exact circuit-level lower-order terms remain modeled.
5. **Memory and data movement.** The arithmetic count assumes extreme streaming and says nothing practical about an astronomical computation.
6. **Seeded transfer.** The local-separator step does not repair the parent's ideal-XOF/ideal-cipher versus fixed-primitive caveat.

## Recommended claim language

The defensible current statement is:

> Under the inherited public-random model and Boolean-gate accounting, the local-separator split-before-lift algorithm has a reproducible Level-I diagnostic of `2^140.384` gates. The result survives the parent's weaker root-probability bound and an eightfold increase in the claimed online-product envelope, but still requires specialist verification of the parametric degree theorem, graph specialization, and finite online schedule.

Calling this a measured, practical, or unconditional fixed-seed break would still be incorrect.
