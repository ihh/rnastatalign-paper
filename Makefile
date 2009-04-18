all: paper.pdf.open

include Makefile.rectest

LATEX := pdflatex
BIBTEX := bibtex

%.pdf: %.tex
	test -e $*.aux && rm $*.aux || eval
	$(LATEX) $*
	-$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*

%.open: %
	open $<

.SECONDARY:
