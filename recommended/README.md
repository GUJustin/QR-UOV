# Recommended QR-UOV paper split

This directory develops a shorter main paper while preserving the canonical
source and PDF in the parent directory. The revision is being made one
reviewable slice at a time.

- `QRUOV_main.tex` contains the protected original abstract and the revised
  main text.
- `QRUOV_supplement.tex` contains the original appendices and technical
  material moved out of the main paper.
- `CONTENT_LEDGER.md` records the central thesis, preservation rules, and the
  disposition of every section.
- `output/pdf/` contains the built PDFs and LaTeX build files.
- `tmp/pdfs/` is reserved for rendered page images used during visual checks.

Run `make` from this directory to build both documents and resolve references
between them. The original source and PDF in the parent directory are not used
as build outputs and remain unchanged.
