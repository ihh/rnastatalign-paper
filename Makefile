all: paper.pdf.open

include Makefile.rectest

LATEX := pdflatex
BIBTEX := bibtex

PRODUCTION := paper s1 s2 s3 s4

all: $(addsuffix .pdf,$(PRODUCTION)) latex.zip

all-open: $(addsuffix .pdf.open,$(PRODUCTION))

TEXFILES := paper.tex paper.bbl algorithm2e.sty
latex.zip: $(TEXFILES)
	zip $@ $(TEXFILES)

%.clean:
	rm $*.aux $*.bbl $*.blg $*.log $*.pdf

%.pdf: %.tex
	test -e $*.aux && rm $*.aux || eval
	$(LATEX) $*
	-$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*

%.open: %
	open $<

.SECONDARY:
