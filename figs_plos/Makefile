PDFLATEX := pdflatex
LATEX    := latex
DVIPS    := dvips
CONVERT  := convert

TEXFIGS := inside insidetransemit cyk cyktransemit cyktraceback outside outsidetransemit

all: $(addsuffix .eps,$(TEXFIGS))

clean: $(addsuffix .clean,$(TEXFIGS))

FIGS := $(addprefix fig,1 2 3 4 5 6 7)
numbered: $(addsuffix .dvi,$(FIGS)) $(addsuffix .eps,$(FIGS)) $(addsuffix .cmyk.eps,$(FIGS)) $(addsuffix .pdf,$(FIGS)) $(addsuffix .fig-cleanup,$(FIGS))

%.clean:
	rm $*.aux $*.bbl $*.blg $*.log $*.pdf $*.eps $*.cmyk.eps

%.pdf: %.tex algdefs.tex
	test -e $*.aux && rm $*.aux || eval
	$(PDFLATEX) $*
	$(PDFLATEX) $*

%.cmyk.eps: %.tex algdefs.tex
	test -e $*.aux && rm $*.aux || eval
	$(LATEX) $*
	$(LATEX) $*
	$(DVIPS) -Tpdf -E -o $@ $*

%.eps: %.cmyk.eps
	$(CONVERT) -colorspace RGB $< $@ || eval

%.open: %
	open $<

%.fig-cleanup:
	rm $*.tex $*.aux $*.log

fig1.%: inside.%
	cp $< $@

fig2.%: insidetransemit.%
	cp $< $@

fig3.%: cyk.%
	cp $< $@

fig4.%: cyktransemit.%
	cp $< $@

fig5.%: cyktraceback.%
	cp $< $@

fig6.%: outside.%
	cp $< $@

fig7.%: outsidetransemit.%
	cp $< $@

.SECONDARY:
