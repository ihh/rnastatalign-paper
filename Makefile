all: paper.pdf.open

include Makefile.rectest

LATEX := pdflatex
BIBTEX := bibtex

PRODUCTION := paper s1 s2 s3 s4

all: $(addsuffix .pdf,$(PRODUCTION))

all-open: $(addsuffix .pdf.open,$(PRODUCTION))

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
