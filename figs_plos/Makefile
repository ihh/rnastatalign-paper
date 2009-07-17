LATEX := pdflatex

TEXFIGS := inside insidetransemit cyk cyktransemit cyktraceback outside outsidetransemit

all: $(addsuffix .pdf,$(TEXFIGS))

clean: $(addsuffix .clean,$(TEXFIGS))

%.clean:
	rm $*.aux $*.bbl $*.blg $*.log $*.pdf

%.pdf: %.tex algdefs.tex
	test -e $*.aux && rm $*.aux || eval
	$(LATEX) $*
	$(LATEX) $*

%.open: %
	open $<

.SECONDARY:
