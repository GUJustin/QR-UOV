# Reduced-parameter end-to-end local-separator forgery

## Outcome

The bundle now contains an executable reduced-parameter instantiation of the attack claimed in `QRUOV_local_separator_upgrade.pdf`. It generates a fresh QR-UOV-style key pair, verifies an honest signature, constructs a forgery without passing the secret key to the attack routine, and verifies the forgery using the unchanged public verification map.

The implementation is intentionally small enough to run and exhaustively cross-check. It is evidence that the individual algebraic components compose into one finite attack. It is **not** evidence that the Level-I gate estimate is accurate, practical, or memory-feasible.

## Reduced QR-UOV instance

The default instance uses

- `q = 13`, extension degree `ell = 3`;
- compact dimensions `V = 3`, `O = 2`;
- expanded dimensions `v = 9`, `m = o = 6`, `n = 15`;
- residual affine core dimension `a = m - 4 = 2`;
- `C = 2^a = 4`, `C0 = 2^(a-1)(a+2) = 8`, and lift precision `P = 4C0+1 = 33`.

The key generator follows the QR-UOV algebraic construction:

1. work over `K = F_(q^3) = F_q[z]/(z^3-z-1)`;
2. sample `m` symmetric `N x N` central matrices over `K` with zero oil-oil block;
3. use the compact extension-linear secret transform `S = [[I,G],[0,I]]`;
4. form `S^T F_i S` over `K`; and
5. pull each form to `F_q` through the trace pairing.

This is an algebraically faithful reduced QR-UOV instance, but it does not reproduce the official byte serialization, seed expansion, or exact hash-domain conventions. The simplified SHAKE256 hash-to-field interface preserves the public relation that the attack must invert.

## Instantiated attack path

`implementation/qruov_local_separator_forgery.py` performs the following path end to end:

1. normalize the chosen hash target to `(0,...,0,1)` by an invertible output transformation;
2. construct an anchored common totally isotropic subspace for three transformed forms using the Chevalley-Warning procedure;
3. choose an affine projective chart and form the `a=m-4` residual quadrics and scale form `H0`;
4. build the Vandermonde product start system with all `C=2^a` known branches;
5. lift every branch coefficientwise along `(1-t)g+t f` to precision `P`;
6. choose a local separator in `F_(q^2)`, form the filter series, and recombine branches by a terminal product tree into `q`, `v_b1`, and `v_psi`;
7. rationally reconstruct the coefficient functions and specialize them valuation-safely at `t=1`;
8. factor the tiny target eliminant, require a simple root, and apply the exact one-sketch Frobenius test and nonzero-square filter;
9. rerun the lifted branches once to reconstruct the accepted coordinates;
10. recover the affine scale and verify the candidate against both the transformed and original public maps.

The exhaustive enumeration of the tiny residual core is performed **after** the attack has accepted a candidate. It is used only as an independent regression check that the target eliminant contains every extension-field root and that the accepted simple separator value has a unique fiber. It is not used to find the forgery.

## Canonical successful run

Run:

```bash
cd implementation
python3 qruov_local_separator_forgery.py
python3 verify_qruov_local_separator_run.py qruov_local_separator_forgery_run.json
```

The fixed seed `20260731` produced:

| Quantity | Result |
|---|---:|
| Honest signature accepted | yes |
| Forgery produced | yes |
| Forgery accepted by public verifier | yes |
| Salt attempts | 1 |
| Attack records | 1 |
| Target eliminant degree | 4 |
| Branch lifts | 8 (4 initial + 4 coordinate rerun) |
| Lifted branch residual checks | 4 |
| Rational reconstructions | 31 |
| Coordinate reruns | 1 |
| Simple target roots accepted by Frobenius | 1 |
| Exact candidates verified | 1 |
| Exhaustive `F_q` residual roots | 2 |
| Exhaustive `F_(q^2)` residual roots | 2 |
| Python wall time | 0.415 seconds |

The public transcript contains the public matrices, message, salt, forged signature, isotropic subspace, chart, affine residual system, separator, specialized eliminants, accepted root, recovered coordinates, and operation trace.

`verify_qruov_local_separator_run.py` reads only this public transcript, reconstructs the residual system, reruns the complete split-before-lift/local-separator core with the recorded separator, and verifies the final signature. The replay succeeds.

## Fixed regression suite

`run_qruov_local_separator_regression.py` ran five consecutive deterministic seeds `20260000` through `20260004` under bounded retry budgets. Results:

- five of five runs produced publicly valid forgeries;
- five of five public-only certificate replays succeeded;
- the largest salt count was two;
- the slowest individual Python run took 10.252 seconds;
- total suite time, including the public replays, was 30.958 seconds.

This is a regression suite, not a statistical estimate of the full attack success probability.

A supplementary run over `q=29` with the same dimensions also succeeded and publicly replayed. Its target eliminant had degree three and the independent exhaustive check found one base-field root and three `F_(q^2)` roots.

## Files

- `implementation/qruov_local_separator_forgery.py`: reduced scheme and complete attack.
- `implementation/QRUOV_cw_helper.py`: exact base-field linear algebra and Chevalley-Warning helpers.
- `implementation/verify_qruov_local_separator_run.py`: public-only replay checker.
- `implementation/run_qruov_local_separator_regression.py`: fixed multi-instance regression.
- `implementation/qruov_local_separator_forgery_run.json`: canonical public transcript.
- `implementation/qruov_local_separator_regression.json`: regression summary.
- `implementation/regression_transcripts/`: five complete public transcripts.
- `implementation/qruov_local_separator_forgery_q29.json`: supplementary field-size run.

## What this resolves

The reviewer criticism that the artifact checked only isolated micro-lemmas while never instantiating the attack is no longer accurate at reduced parameters. The code supplies one composed finite algorithm and demonstrates a valid forgery from public data.

## What it does not resolve

The following remain open:

1. The implementation uses naive truncated-series and product-tree arithmetic at tiny `a=2`; it does not instantiate the asymptotically fast dyadic multiplication schedule used in the Level-I ledger.
2. It does not test the full NIST parameters, which are astronomically beyond execution.
3. It does not validate the public-random-model transfer to the official seeded key distribution.
4. It does not measure memory traffic or establish the inherited Boolean-gate model as the right practical model.
5. It does not replace the need for specialist review of the characteristic-polynomial degree bound on the exact saturated dominant homotopy algebra.

Accordingly, the implementation upgrades the work from a collection of component checks to a genuine reduced-parameter attack prototype, but the full-parameter result remains a model-based complexity claim.
