LATEX := pdflatex

TEXFIGS := inside insidetransemit cyk cyktransemit cyktraceback outside outsidetransemit

all: $(addsuffix .eps,$(TEXFIGS))

clean: $(addsuffix .clean,$(TEXFIGS))

%.clean:
	rm $*.aux $*.bbl $*.blg $*.log $*.pdf

%.eps: %.pdf
	convert $< $@

%.pdf: %.tex algdefs.tex
	test -e $*.aux && rm $*.aux || eval
	$(LATEX) $*
	$(LATEX) $*

%.open: %
	open $<

.SECONDARY:
