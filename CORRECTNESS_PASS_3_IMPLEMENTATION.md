# Correctness pass 3: reduced end-to-end implementation

## Question addressed

A reviewer objected that the artifact checked isolated lemmas and ledgers without ever instantiating the claimed attack. This pass asks the narrow, executable question:

> Can the exact local-separator pipeline be composed into a finite program that forges a signature on a fresh reduced QR-UOV instance, without using the secret key in the attack?

The answer is yes.

## What was implemented

The new code implements a reduced QR-UOV-style scheme over `F_(q^3)` with compact secret transform and trace pullback, honest signing, and unchanged public verification. The attack routine receives only the public matrices, message, and salt and runs:

1. target normalization;
2. anchored Chevalley-Warning seed and isotropic-subspace construction;
3. projective chart and residual affine core;
4. product-start homotopy and coefficientwise lifting of every branch;
5. local-separator recombination for the characteristic polynomial, Frobenius sketch, and rational oil filter;
6. rational reconstruction and valuation-safe specialization;
7. simple-root and exact Frobenius descent;
8. one coordinate rerun;
9. affine scaling; and
10. verification under the original public map.

The post-acceptance exhaustive enumeration is a cross-check, not a solver.

## Canonical result

For `q=13, ell=3, V=3, O=2`, seed `20260731`, the implementation:

- generated an honest signature accepted by the verifier;
- produced a forgery in the first attack record and first salt;
- passed the unchanged public verification map;
- reran publicly from the serialized transcript without the secret key;
- lifted all four start branches to precision 33 and checked all four homotopy residuals;
- performed one coordinate rerun and 31 rational reconstructions; and
- matched exhaustive base- and extension-field root enumeration after acceptance.

The successful transcript is `implementation/qruov_local_separator_forgery_run.json`.

## Regression

Five consecutive deterministic seeds `20260000` through `20260004` all produced publicly valid forgeries and all five public-only replays succeeded. A separate `q=29` run also forged and replayed successfully.

These runs show algorithmic composition and guard against a one-transcript accident. They are not used to infer a full-parameter success probability.

## Correctness checks embedded in the implementation

- Extension-field public/secret consistency on random points during key generation.
- Exact lifted-branch residual check modulo `t^P` for every initial branch.
- Direct verification of every rational reconstruction against all known Taylor coefficients.
- Equality of the recomputed coordinate-rerun eliminant with the first specialized eliminant.
- Exact residual-system and scale-form verification of reconstructed coordinates.
- Verification under both the transformed target-normalized public map and the original public map.
- Exhaustive enumeration of the tiny affine core over `F_q` and `F_(q^2)` after acceptance.
- Public-only transcript replay.

## What the implementation proves

It refutes the narrow objection that no finite instance of the claimed attack had been constructed. The branch lifting, recombination, specialization, Frobenius descent, coordinate recovery, and public verification now compose in executable code.

## What it does not prove

- It does not establish the Level-I Boolean-gate estimate.
- It does not implement the asymptotically fast full-scale multiplication schedule.
- It does not test official full parameters.
- It does not establish seeded-distribution transfer.
- It does not resolve memory traffic.
- It does not replace specialist review of the parametric characteristic-polynomial degree theorem on the exact saturated dominant algebra.

## Assessment

The reviewer’s original criticism was fair for the previous artifact but is no longer accurate for the updated bundle. A more precise remaining criticism would be:

> The artifact now instantiates and publicly verifies the complete attack at reduced parameters, but it does not yet instantiate the fast full-parameter schedule underlying the quoted `2^140.384` gate projection.

That criticism remains valid and is stated explicitly in the paper.
