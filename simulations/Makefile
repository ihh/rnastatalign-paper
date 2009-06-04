
# there has to be a better way of doing this...
#MAKEFILE := /mnt/nfs/users/rbradley/projects/rnastatalign-paper/simulations/Makefile
#MAKEFILE := /Users/rbradley/papers/rnastatalign-paper/simulations/Makefile
MAKEFILE := $(CURDIR)/Makefile

EVOLSAYER = $(DARTDIR)/perl/evolsayer.pl
EVOLSAYER_ARGS = -sub $(DARTDIR)/grammars/pfold-mix80.eg -looplen 9 -stemlen 5 -stem 0.01

QSUB := /opt/sge/bin/lx24-amd64/qsub
QSUB_ARGS := -cwd -V -l qname=himem_nodes.q -P stemloc

USE_SGE := 

# Run evolsayer to create all simulated alignments of the form .reference.stock
# as well as indiegram input files of the form .input.stock.
alignments/all:
	@echo "Making $@"
	@./run-evolsayer.pl $(@D)
	touch $@

# Run indiegram on a single simulated alignment.
# stem: '0.01/1'
inferred/%.indiegram.stock: alignments/%.input.stock
	@echo "Making $@"
	$(DARTDIR)/bin/indiegram $< -log 7 -logfile $(subst .stock,.log,$@) >$@

# Run indiegram on all simulated alignments.
inferred/all: alignments/all
	@echo "Making $@"
	@mkdir -pv $(@D)
	@for directory in `find $(<D)/ -mindepth 1 -type d ! -iname "*cvs*" -print | xargs -n 1 basename`;\
	do\
		mkdir -p $(@D)/$$directory;\
		for file in `find $(<D)/$$directory -iname "*.input.stock" -print | xargs -n 1 -I{} basename {} .input.stock`;\
		do\
			if (! test -s $(@D)/$$directory/$$file.indiegram.stock && ! test -e $(@D)/$$directory/$$file.indiegram.submission);\
			then\
				rm -fv $(@D)/$$directory/$$file.indiegram.stock $(@D)/$$directory/$$file.indiegram.submission;\
				if (test -n "$(USE_SGE)");\
				then\
					echo "farming out $(@D)/$$directory/$$file.indiegram.stock (himem queue)...";\
					$(QSUB) $(QSUB_ARGS) -N m-$$directory-$$file -b y -o $(@D)/$$directory/$$file.indiegram.sge-out -e $(@D)/$$directory/$$file.indiegram.sge-err\
					 make -f $(MAKEFILE) $(@D)/$$directory/$$file.indiegram.submission >$(@D)/$$directory/$$file.indiegram.sge-jid;\
				else\
					make -f $(MAKEFILE) $(@D)/$$directory/$$file.indiegram.submission;\
				fi;\
			else\
				echo "found $(@D)/$$directory/$$file.indiegram.stock, not re-building...";\
			fi;\
		done;\
	done
	@touch $@

# Wrapper rule.
# Touches itself to prevent re-execution in case of job failure.
# Sets a .failed flag if an alignment job doesn't run successfully
# (the || means that the second half of the command gets executed only if
# the first doesn't succeed, ie has a non-zero exit status).
# If the job ran successfully, then discard all SGE and logging information.
# Stores timing information.
# Stem is e.g. $(WONGAU)/family/myseqs/myseqs.fsa.gf1.
inferred/%.submission:
	@touch $@
	@echo "Job started running on host $$HOSTNAME at `date`." >$(subst .submission,.host,$@)
	make -f $(MAKEFILE) $(subst .submission,.stock,$@) || touch $(subst .submission,.failed,$@)
	if (test -s $(subst .submission,.stock,$@) && ! test -e $(subst .submission,.failed,$@));\
	then\
		make -f $(MAKEFILE) $(subst .submission,.timing,$@);\
	fi

# Get approximate time (in seconds) required to run job.
# Uses the difference in modification times between the 
# .sge-host file, which is created when the job starts running,
# and the .stock file, which is last modified when the job finishes.
# Stem is e.g. $(WONGAU)/family/myseqs.fsa.gf1.
inferred/%.timing: 	inferred/%.stock\
			inferred/%.host
	@echo "Making $@"
	@echo $$[ $(shell stat --format "%Y" $<) - $(shell stat --format "%Y" $(word 2,$^)) ] >$@

# Remove .submission tags for jobs with a .failed tag.
# Also cleans 'all' tags.
# The submission tags indicate that a job (may) be currently running,
# and serve to prevent re-submission of identical jobs.
# This rule must be called before failed jobs can be re-submitted.
inferred/clean-failed-submission-tags:
	@echo "Making $@"
	@find $(@D) -type f -name "all" | xargs -n 50 rm -fv
	@find $(@D) -type f -name "*.failed" | perl -ne'chomp; $$_ =~ s/\.failed//; print "$$_.failed $$_.submission ";' | xargs -n 50 rm -fv

# Remove all .submission tags.
# Also cleans 'all' tags.
inferred/clean-all-submission-tags:
	@echo "Making $@"
	@read -p "You are about to clean all submission tags.  Do you want to continue? [yes/No]: ";\
	if [ "$$REPLY" = "yes" ];\
	then\
		echo "Removing files...";\
		find $(@D) -type f -name "all" | xargs -n 50 rm -fv;\
		find $(@D) -type f -name *.submission -print0 | xargs -0 -n 50 rm -fv;\
		find $(@D) -type f -name *.failed -print0 | xargs -0 -n 50 rm -fv;\
	else\
		echo "Aborting operation";\
	fi

# Find all failed jobs.
inferred/find-failed:
	@find $(@D) -type f -name *.failed | perl -pe's/\.failed//;'

# Find all empty Stockholm alignments.
# This rule is slow!
inferred/find-empty:
	@find $(@D) -type f -name *.stock | perl -MStockholm -ne'$$file = $$_; chomp $$file; $$stk=Stockholm->from_file("$$file"); if (@{$$stk->seqname} == 0) { print "$$file\n";}'

# Remove all empty Stockholm alignments.
# Cleans corresponding .submission tags as well.
inferred/clean-empty:
	@echo "Making $@"
	@find $(@D) -type f -name *.stock | perl -MStockholm -ne'$$file = $$_; chomp $$file; $$stk=Stockholm->from_file("$$file"); if (@{$$stk->seqname} == 0) { $$file =~ s/\.stock//; print "$$file.stock $$file.submission\n"; }' | xargs -n 50 rm -fv


# Run handalign (TKF91 model) on a single simulated alignment.
# stem: '0.01/1'
inferred/%.tkf.stock: alignments/%.reference.stock
	@echo "Making $@"
	./run-tkf.pl $< >$@

# Run handalign (TKF91 model) on all simulated alignments.
inferred/tkf:
	make -f $(MAKEFILE) $(subst alignments,inferred,$(subst reference,tkf,$(wildcard alignments/*/*.reference.stock)))


# Run handalign (long indel model) on a single simulated alignment.
# stem: '0.01/1'
inferred/%.long.stock: alignments/%.reference.stock
	@echo "Making $@"
	./run-long.pl $< >$@

# Run handalign (long indel model) on all simulated alignments.
inferred/long:
	make -f $(MAKEFILE) $(subst alignments,inferred,$(subst reference,long,$(wildcard alignments/*/*.reference.stock)))



# Prepare a stemloc input alignment (indiegram format is not strict flush Stockholm, as unaligned input sequences are not separated by "# STOCKHOLM 1.0" delimiter lines)
alignments/%.stemloc.input.stock: alignments/%.input.stock
		@echo "Converting $< (indiegram Stockholm dialect) into $@ (pared-down Stockholm for stemloc)"
	cat $< | perl -ne 'if(/#=/){if(/#=GR/){print;print"\n//\n\n# STOCKHOLM 1.0\n"unless++$$seqs>=3}}else{print}' >$@

# Run stemloc on a single simulated alignment.
# stem: '0.01/1'
inferred/%.stemloc.stock: alignments/%.stemloc.input.stock
	@echo "Making $@"
	stemloc -na -1 --global --noshow-intermediates --multiple $< >$@

# Run stemloc on all simulated alignments.
inferred/stemloc:
	make -f $(MAKEFILE) $(subst alignments,inferred,$(subst reference,stemloc,$(wildcard alignments/*/*.reference.stock)))


# Run stemlocama on a single simulated alignment.
# stem: '0.01/1'
inferred/%.stemlocama.stock: alignments/%.stemloc.input.stock
	@echo "Making $@"
	stemloc -na -1 --global --noshow-intermediates --ama $< >$@

# Run stemlocama on all simulated alignments.
inferred/stemlocama:
	make -f $(MAKEFILE) $(subst alignments,inferred,$(subst reference,stemlocama,$(wildcard alignments/*/*.reference.stock)))



# Analysis.
inferred/%.indiegram.analysis:	alignments/%.reference.stock\
				inferred/%.indiegram.stock\
				inferred/%.tkf.stock\
				inferred/%.stemloc.stock\
				inferred/%.stemlocama.stock
	@echo "Making $@"
	./compare.pl $^ >$@

# Collate results.
inferred/%/indiegram.analysis: $(addsuffix .indiegram.analysis,$(addprefix inferred/%/,$(shell ./gseq 1 25)))
	@echo "Making $@"
	rm -fv $@
	@for i in $^;\
	do\
		perl -e 'print "$*\t", `cat @{[shift]}`' $$i >>$@;\
	done

# Collate results.
inferred/indiegram.analysis: inferred/all inferred/tkf inferred/long inferred/stemloc inferred/stemlocama
	@echo "Making $@"
	perl -e '@s=qw(AlignAcc AlignSn AlignPPV);print join ("\t", @s, "StructureReconstruction"); for $$a (qw(TKF Long Stemloc StemlocAma)) { print map ("\t$$_$$a", @s) } print "\n"' >$@
	@for directory in `find $(<D)/ -mindepth 1 -type d ! -iname "*cvs*" -print | xargs -n 1 basename`;\
	do\
		make -f $(MAKEFILE) $(@D)/$$directory/indiegram.analysis;\
		cat $(@D)/$$directory/indiegram.analysis >>$@;\
	done



.SECONDARY:

