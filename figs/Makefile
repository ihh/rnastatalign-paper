
%.xrate.stock: %.stock
	cat $< | calc_ss_cons.pl | perl -pe 'if(/SS_cons/){s/\./_/g}' | xrate -ar -g $(DARTDIR)/grammars/pfold-mix80.eg >$@
