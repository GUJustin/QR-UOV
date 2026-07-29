# Cross-paper framing update

The QR-UOV paper now presents the MAYO, QR-UOV, and SNOVA manuscripts as three concurrent, self-contained, technically independent papers in one structural-cryptanalysis program.

Changes:

1. Added an introduction paragraph that states the shared high-level lesson and explicitly says that none of the papers is required for the correctness of the others.
2. Added a QR-UOV-specific comparison at the start of Related Work distinguishing projective extension-field slicing and graph completion from the MAYO polar quotient/Hodge-contraction attacks and the SNOVA symmetric-square quotient.
3. Added two mutually cited placeholder bibliography records with the identical author set and order: Justin Thaler, Kai Zhang, Scott Kominers, Quang Dao.
4. Did not use formal "companion paper," serial-title, or Part I/II/III branding.
5. Isolated the two provisional titles in `\MAYOPaperTitle` and `\SNOVAPaperTitle` macros near the top of `main.tex` for easy replacement.

No theorem, numerical estimate, or attack claim was changed.
