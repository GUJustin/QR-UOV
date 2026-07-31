# QR-UOV integrated e-print bundle

This is the full integrated manuscript. The overwrite-ready paper source is:

- `QRUOV_strict_final.tex`

It preserves the original full paper and integrates the local-separator attack, corrected Level-I/III/V ledgers, reduced-parameter end-to-end forgery implementation, public replay experiment, limitations, and updated conclusions.

Build:

```bash
pdflatex QRUOV_strict_final.tex
pdflatex QRUOV_strict_final.tex
```

The `implementation/` directory contains the reduced-parameter attack and public replay scripts. The full-parameter cost remains a model-based extrapolation as stated in the manuscript.
