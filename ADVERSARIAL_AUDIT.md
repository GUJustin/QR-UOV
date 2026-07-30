# Adversarial audit - local-separator split-before-lift candidate

## Current status

No fatal algebraic defect is known after two adversarial passes. The v2 primary diagnostics are:

- Level I: `140.383995063` gates, margin `2.616004937` bits;
- Level III: `191.027133323`, margin `15.972866677` bits;
- Level V: `247.125438947`, margin `24.874561053` bits.

These are symbolic Boolean-gate diagnostics under inherited assumptions, not implementation measurements.

## Arguments that survived review

### Generic characteristic polynomial at a colliding separator

The eliminant is the characteristic polynomial

`det(U I - sum Lambda_j M_{x_j})`.

Specializing `Lambda=lambda` is valid whether or not `lambda` separates all branches. Over the rational-function field in independent `Lambda` variables, the generic linear element separates the distinct points of the finite etale generic fiber. Schost's parameter-degree bound can therefore be applied after transcendental scalar extension without choosing a large-field separator. Differentiation in `Lambda_j` gives coordinate and scalar interpolation polars without increasing the parameter degree.

### Local bad-condition count

If `I` branches escape and `F` remain finite, counted with multiplicity, then `I+F=C`. Relative to a fixed eligible target point, there are at most `I` leading-valuation cancellation conditions and `F-1` target-collision conditions. Each has at most a `1/|K'|` fraction of `K'`-coefficient vectors, even when its coefficient vector lives over an algebraic closure. The total failure bound is `(C-1)/|K'|`.

### Exact simple-factor descent

For `b1(P)=lambda^(q^{-1}).P`, one has `b1(P)^q=lambda.P^q`. If a point represented by a simple target factor passes `b1(P)^q=lambda.P` but is not base-field rational, then it and its Frobenius conjugate are distinct represented target points with the same separator value. The factor would be repeated, a contradiction.

### One coordinate rerun

Every scalar-passing simple factor is already a genuine base-field target point, and the rational oil filter is exact at that point. Thus failed trials perform no coordinate rerun and the successful execution performs exactly one. The v2 ledger now separately charges both the branch rerun and its separator-leaf formation.

## Defects found and fixed

1. **Infinite-branch proof overstatement.** Escaping numerator terms can survive common normalization. They are harmless because they are divisible by the finite eliminant, not because each vanishes. Both notes and a new exact toy now reflect this.
2. **Omitted coordinate leaf cost.** Added explicit one-time formation of `lambda.x_s` during coordinate recovery. It changes Level I by about `3.14e-5` bits.
3. **Root-bound derivation omitted.** The second-Bonferroni calculation is now written out and independently checked.
4. **Filter-offset sign.** Corrected the forbidden set to `{-lambda(P)}`.
5. **Finite-field proof wording.** Replaced an informal base-extension statement by a direct generic-`Lambda` primitive-element argument.

## Independent checks

- `verify_local_separator_ledger.py` imports neither audit module and reproduces all primary totals to below `5e-10` bits.
- Rabin irreducibility certificates pass for all used sparse extension moduli.
- `toy_local_separator_specialization.py` checks specialization of the generic characteristic polynomial at a colliding linear form.
- `toy_exact_simple_frobenius.py` exhaustively checks the simple-factor Frobenius collision criterion over `F_9/F_3`.
- `toy_infinite_branch_specialization.py` checks the corrected nonvanishing-but-divisible infinite-branch case.
- `verify_root_bonferroni.py` checks the sharpened root bound and monotonicity numerically.

## Robustness

- Parent root lower bound fallback at Level I: `140.531460701`.
- Rational-reconstruction multiplier `128`: approximately `141.815` under the sharpened root bound.
- Whole online-envelope multiplier `8`: `142.856480953` after reoptimization.
- Whole online-envelope multiplier `16`: `143.806489657`, above target.

The online schedule is therefore the most important finite-constant obligation.

## Items still requiring independent verification before any public break claim

1. Match the exact saturated dominant homotopy ideal and generic fiber to Schost's theorem.
2. Audit the rational-filter graph closure, saturation, degree `2C0`, and parameter denominators.
3. Supply a literal finite online-multiplication schedule or a count-certified prototype.
4. Audit coefficientwise rational reconstruction, factorization, squarefree extraction, and public verification costs.
5. Preserve the parent's seeded-instance and fixed-primitive caveats exactly.
6. Treat memory and data movement separately from Boolean arithmetic.

## Bottom line

The second pass made the claim more defensible and slightly more expensive. It did not overturn the candidate Level-I crossing. The current result deserves specialist computational-algebra review, but it is still premature to call it an unconditional or practical QR-UOV break.
