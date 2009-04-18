#!/usr/bin/perl -w


# The following one-liner plots postprob against probability of correctness
# (the third column is the error in the estimate of prob-of-correctness;
#  actually it's the standard deviation of the beta distribution)

# Replace "fields 4 8" by "fields 4 10" to get the postprob estimates for the naive grammar.
# Replace "cat" by "grep pfoldF" to get stats for paired columns only.

#  cat Published.xval-results | fields 4 8 | perl -e '$mul=shift;while(<>){($c,$p)=split;$bin=int($p/$mul+.5);++$c{$bin}if$c;++$t{$bin}}for$bin(sort{$a<=>$b}keys%t){$a=$c{$bin};$b=$t{$bin}-$c{$bin};$p=$a/($a+$b);$perr=sqrt($a*$b/(($a+$b)*($a+$b)*($a+$b+1)));print$bin*$mul," $p $perr\n"}' .1



# The following one-liners break down incorrect basepair reconstructions
# The first number printed is the number of incorrect reconstructions; the second, the number of these which were reconstructed as noncanonical basepairs

# For the naive grammar:
# grep "n= 0" Published.xval-results | perl -e 'while(<>){if(/pfoldF ..\/..\/(..)/){++$c{lc$1};++$t}}print$t," ",$t-$c{ua}-$c{au}-$c{cg}-$c{gc}-$c{gu}-$c{ug},"\n"'

# For the main grammar:
# grep "g= 0" Published.xval-results | perl -e 'while(<>){if(/pfoldF ..\/(..)\/../){++$c{lc$1};++$t}}print$t," ",$t-$c{ua}-$c{au}-$c{cg}-$c{gc}-$c{gu}-$c{ug},"\n"'


# The following one-liners break down incorrect basepair reconstructions by postprob

# grep "g= 0" Published.xval-results | perl -e '$mul=shift;while(<>){if(/pfoldF ..\/(..)\/.. .*gpostprob= (\S+)/){($chars,$pp)=($1,$2);next if$chars=~/[\-\.]/;$bin=int($pp/$mul+.5);++$cb{$bin}->{lc$chars};++$t{$bin}}}for$bin(sort{$a<=>$b}keys%cb){%c=%{$cb{$bin}};$bad=$t{$bin}-$c{ua}-$c{au}-$c{cg}-$c{gc}-$c{gu}-$c{ug};print$bin*$mul," ",$bad/$t{$bin},"\n"}' .1

# grep "n= 0" Published.xval-results | perl -e '$mul=shift;while(<>){if(/pfoldF ..\/..\/(..) .*npostprob= (\S+)/){($chars,$pp)=($1,$2);next if$chars=~/[\-\.]/;$bin=int($pp/$mul+.5);++$cb{$bin}->{lc$chars};++$t{$bin}}}for$bin(sort{$a<=>$b}keys%cb){%c=%{$cb{$bin}};$bad=$t{$bin}-$c{ua}-$c{au}-$c{cg}-$c{gc}-$c{gu}-$c{ug};print$bin*$mul," ",$bad/$t{$bin},"\n"}' .1



use Stockholm::Database;
use Newick;

my $maxcols = 500;  # alignments with more columns than this are ignored
my $maxrows = 400;  # alignments with more rows than this are ignored

my $gram = $ENV{DARTDIR} . "/grammars/pfold.eg";
my $naive = $ENV{DARTDIR} . "/grammars/rev.eg";
my $trials = 1;
my @seqname;
my $usage = "Usage: $0 [-g <covariant RNA grammar>] [-n <naive point-substitution grammar>] [-t trials] [-m maxcols] [-r maxrows] [-seqname seqname]  <Stockholm database>\n";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-g") { defined($gram = shift) or die $usage }
    elsif ($arg eq "-n") { defined($naive = shift) or die $usage }
    elsif ($arg eq "-t") { defined($trials = shift) or die $usage }
    elsif ($arg eq "-m") { defined($maxcols = shift) or die $usage }
    elsif ($arg eq "-r") { defined($maxrows = shift) or die $usage }
    elsif ($arg eq "-seqname") { my $seqname; defined($seqname = shift) or die $usage; push @seqname, $seqname }
    else { push @argv, $arg }
}
die $usage unless @argv == 1;
my ($filename) = @argv;

my $stockdb = Stockholm::Database->from_file ($filename);
my %seqname2stock;
for my $stock (@$stockdb) {
    for my $seqname (@{$stock->seqname}) {
	die "Duplicate seqname $seqname" if exists $seqname2stock{$seqname};
	$seqname2stock{$seqname} = $stock;
    }
}

my %seen;
for (my $trial = 1; $trial <= $trials; ++$trial) {

    # choose an alignment & row
    my ($stock, $seqname);

    # name(s) pre-specified?
    if (@seqname) {
	$seqname = shift @seqname;
	$stock = $seqname2stock{$seqname};
	die "Sequence name '$seqname' unknown" unless defined $stock;

    } else {
	# choose a random alignment & row
	$stock = $stockdb->[int rand @$stockdb];
	$seqname = $stock->seqname->[int rand scalar @{$stock->seqname}];
    }

    warn "Trial $trial of $trials: sequence $seqname\n";

    if ($stock->columns > $maxcols) {
	warn "...too many columns; skipping\n";
	--$trial;
	next;
    }

    if (@{$stock->seqname} > $maxrows) {
	warn "...too many rows; skipping\n";
	--$trial;
	next;
    }

    if ($seen{$seqname}) {
	warn "...already seen; skipping\n";
	--$trial;
	next;
    }
    ++$seen{$seqname};

    # get true row
    my $true = $stock->seqdata->{$seqname};

    # clone alignment; replace sequence with Felsenstein wildcards
    my $copy = $stock->copy;
    $copy->seqdata->{$seqname} =~ s/[^\-\.]/\*/g;

    # for speed, remove all single-stranded columns with gaps in reference sequence
    my @dropcols;
    my $ungapped = "";
    my $ss = $stock->gc_SS_cons;
    for my $col (0..$stock->columns - 1) {
	my $ss_char = substr ($ss, $col, 1);
	my $row_char = substr ($true, $col, 1);
	if ($ss_char ne '<' && $ss_char ne '>' && ($row_char eq '.' || $row_char eq '-')) {
	    push @dropcols, $col;
	} else {
	    $ungapped .= substr ($true, $col, 1);
	}
    }
    warn "Dropping columns (@dropcols)\n";
    $copy->drop_columns (@dropcols);

    # create temporary filenames
    my $tmproot = $filename;
    $tmproot =~ s/.*?([^\/]+)$/$1/;

    my ($tmpstock, $tmprec, $tmpnaive)
	= map ("$tmproot.$_",
	       qw(xval-stk xval-rec xval-naive));

    # save alignment to temporary file
    $copy->to_file ($tmpstock);

    # reconstruct
    syswarn ("xrate -g $gram -ar -arpp $tmpstock -log 6 >$tmprec");

    syswarn ("xrate -g $naive -ar -arpp $tmpstock -log 6 >$tmpnaive");

    my $rec = Stockholm->from_file ($tmprec);
    my $naif = Stockholm->from_file ($tmpnaive);

    my (%total, %correct);

    # read postprobs from $naif reconstruction
    my %naive_postprob;
    for my $gs (@{$naif->gs_ancrec_CYK_MAP_PP->{$seqname}}) {
	if ($gs =~ /State (\S+) columns \((\d+)\).*chars \((.)\).* postprob (\S+)/) {
	    my ($state, $col, $char, $postprob) = ($1, $2, $3, $4);
	    $naive_postprob{$col} = $postprob;
	}
    }

    # read postprobs from $rec reconstruction
    my $recstr = $rec->gr_ancrec_CYK_MAP->{$seqname};
    my $naifstr = $naif->gr_ancrec_CYK_MAP->{$seqname};
    for my $gs (@{$rec->gs_ancrec_CYK_MAP_PP->{$seqname}}) {
	if ($gs =~ /State (\S+) columns \(([^\)]+)\).*chars \(([^\)]+)\).* postprob (\S+)/) {
	    my ($state, $cols, $spc_chars, $postprob) = ($1, $2, $3, $4);
	    my @cols = split /\s+/, $cols;
	    my (@real_chars, @chars, @naive_chars);
	    my $naive_postprob = 1;
	    for my $col (@cols) {
		push @real_chars, substr ($ungapped, $col-1, 1);
		push @chars, substr ($recstr, $col-1, 1);
		push @naive_chars, substr ($naifstr, $col-1, 1);
		$naive_postprob *= $naive_postprob{$col} if defined($naive_postprob{$col});
	    }

	    my ($chars, $real_chars, $naive_chars)
		= map (uc(join("",map($_ =~ /[tT]/ ? "u" : $_, @$_))),
		       \@chars,
		       \@real_chars,
		       \@naive_chars);

	    print
		"$seqname $state $real_chars/$chars/$naive_chars g= ",
		$chars eq $real_chars ? 1 : 0,
		" n= ",
		$naive_chars eq $real_chars ? 1 : 0,
		" gpostprob= ", $postprob,
		" npostprob= ", $naive_postprob,
		"\n";
	}
    }
}

sub syswarn {
    my ($cmd) = @_;
    warn "Executing $cmd\n";
    system $cmd;
}
