LATEXMK = latexmk
OUTDIR = output/pdf
FLAGS = -pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir=$(OUTDIR)

.PHONY: all main supplement

all:
	mkdir -p $(OUTDIR)
	$(LATEXMK) $(FLAGS) QRUOV_main.tex
	$(LATEXMK) $(FLAGS) QRUOV_supplement.tex
	$(LATEXMK) $(FLAGS) QRUOV_main.tex
	$(LATEXMK) $(FLAGS) QRUOV_supplement.tex

main:
	mkdir -p $(OUTDIR)
	$(LATEXMK) $(FLAGS) QRUOV_main.tex

supplement:
	mkdir -p $(OUTDIR)
	$(LATEXMK) $(FLAGS) QRUOV_supplement.tex
