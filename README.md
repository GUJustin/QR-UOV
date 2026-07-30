# QR-UOV local-separator candidate bundle - correctness pass v2

This bundle contains the split-before-lift candidate and the stronger local-separator upgrade, together with a second adversarial correctness pass.

## Current headline

Under the inherited Boolean-gate stress model:

| Level | Strict parent | Global split | Local split v2 | Target | Margin |
|---|---:|---:|---:|---:|---:|
| I | 151.638746 | 141.125734 | **140.383995** | 143 | **2.616005 bits** |
| III | 203.506934 | 191.630088 | **191.027133** | 207 | **15.972867 bits** |
| V | 260.389091 | 247.637444 | **247.125439** | 272 | **24.874561 bits** |

The Level-I optimum uses `F_(127^8)` with the Rabin-certified sparse modulus `X^8 + X^4 - 1`.

This remains a finite, model-based candidate attack. It is not a measured implementation and should not yet be presented as a public break.

## What the second pass changed

1. **Corrected target specialization.** An escaping-branch numerator term need not vanish after common valuation normalization. The correct statement is that any surviving such term is divisible by the entire finite target eliminant, so it vanishes at each finite simple root. The new toy `toy_infinite_branch_specialization.py` checks this exact phenomenon.
2. **Added an omitted charge.** The successful coordinate rerun must reform `u_s=lambda.x_s`. The new one-time leaf-formation term has log-cost `124.898` at Level I and changes the total by only `0.000031` bits.
3. **Made the finite-field characteristic-polynomial argument cleaner.** Generic independent coefficients separate the geometric points over the rational-function field, so Schost's degree argument can be applied after transcendental scalar extension without assuming a large finite ground field.
4. **Proved the sharpened root bound in the note.** The `0.3690869822...` value follows directly from the parent's first and second moments by second Bonferroni. `verify_root_bonferroni.py` checks the numerical substitution and monotonicity.
5. **Corrected offset wording.** The forbidden filter offsets are `-lambda(P)`, not `lambda(P)`; the counting bound is unchanged.
6. **Added conservative sensitivities.** Using the parent's already-proved root lower bound gives Level-I total `140.531461`. Multiplying the complete online-product envelope by eight gives `142.856481`; multiplying it by sixteen gives `143.806490`.

## Main mechanism

The global split note paid for one linear form to separate all `C=2^a` homotopy branches. The upgrade instead:

- keeps the full characteristic polynomial even when the chosen linear form collides;
- uses characteristic-polynomial polars for the scalar interpolation numerators;
- requires only valuation genericity at infinity and isolation of one canonical eligible target point, at most `C-1` bad linear conditions;
- samples the filter offset outside the known start forbidden values;
- accepts only simple target factors; and
- uses the exact test `(lambda^(q^{-1}).P)^q = lambda.P`, for which a non-base point would force a Frobenius collision and a repeated factor.

This removes the global `C^2/|K'|` separator event, the second Frobenius numerator, the probabilistic false-positive term, the abort cap, and repeated coordinate recovery.

## Principal files

- `QRUOV_local_separator_upgrade.pdf/.tex` - corrected technical note.
- `QRUOV_local_separator_audit.py/.json` - primary finite ledger with sensitivity data.
- `verify_local_separator_ledger.py` - independent no-import reconstruction of every primary component.
- `CORRECTNESS_PASS_2.md` - findings and remaining blockers from this pass.
- `ADVERSARIAL_AUDIT.md` - cumulative adversarial audit.
- `toy_infinite_branch_specialization.py` - exact check of the corrected specialization argument.
- `toy_local_separator_specialization.py` - colliding characteristic-polynomial specialization toy.
- `toy_exact_simple_frobenius.py` - exhaustive `F_9/F_3` exact-descent toy.
- `verify_root_bonferroni.py` - root-bound numerical check.
- `QRUOV_split_before_lift_note.pdf` - preceding branchwise solver note, retained as supporting background; its TeX source is intentionally omitted because it is not part of the paper itself.

## Reproduce

```bash
python3 QRUOV_local_separator_audit.py
python3 verify_local_separator_ledger.py
python3 verify_root_bonferroni.py
python3 toy_infinite_branch_specialization.py
python3 toy_exact_simple_frobenius.py
python3 toy_local_separator_specialization.py
```

The other symbolic toys require SymPy. Rebuild the paper note with TeX Live:

```bash
pdflatex QRUOV_local_separator_upgrade.tex
pdflatex QRUOV_local_separator_upgrade.tex
```

Only the TeX source for the paper itself is included in this bundle.

## Remaining blockers

The result still needs specialist verification that the characteristic-polynomial parameter-degree theorem applies exactly to the parent's saturated dominant homotopy algebra and to the rational-filter graph. The dyadic online-product envelope also needs a line-by-line finite schedule or a prototype; the new sensitivity table shows that an eightfold error is absorbed but a sixteenfold error is not. Seeded-instance transfer, memory/data movement, and the inherited Boolean-gate model remain separate caveats.
