LATEX := latex
DVIPS := dvips
CONVERT := convert

TEXFIGS := inside insidetransemit cyk cyktransemit cyktraceback outside outsidetransemit

all: $(addsuffix .eps,$(TEXFIGS))

clean: $(addsuffix .clean,$(TEXFIGS))

%.clean:
	rm $*.aux $*.bbl $*.blg $*.log $*.pdf $*.eps $*.cmyk.eps

%.cmyk.eps: %.tex algdefs.tex
	test -e $*.aux && rm $*.aux || eval
	$(LATEX) $*
	$(LATEX) $*
	$(DVIPS) -Tpdf -E -o $@ $*

%.eps: %.cmyk.eps
	$(CONVERT) -colorspace RGB $< $@ || eval

%.open: %
	open $<

.SECONDARY:
