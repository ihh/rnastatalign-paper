#!/usr/bin/perl -w


# The following one-liner plots postprob against probability of correctness
# (the third column is the error in the estimate of prob-of-correctness;
#  actually it's the standard deviation of the beta distribution)

# Replace "fields 4 8" by "fields 4 10" to get the postprob estimates for the naive grammar.
# Replace "cat" by "grep pfoldF" to get stats for paired columns only.

#  cat Hammerhead_1.results | fields 4 8 | perl -e '$mul=shift;while(<>){($c,$p)=split;$bin=int($p/$mul+.5);++$c{$bin}if$c;++$t{$bin}}for$bin(sort{$a<=>$b}keys%t){$a=$c{$bin};$b=$t{$bin}-$c{$bin};$p=$a/($a+$b);$perr=sqrt($a*$b/(($a+$b)*($a+$b)*($a+$b+1)));print$bin*$mul," $p $perr\n"}' .1



# The following one-liners break down incorrect basepair reconstructions
# The first number printed is the number of incorrect reconstructions; the second, the number of these which were reconstructed as noncanonical basepairs

# For the naive grammar:
# grep "n= 0" Hammerhead_1.results | perl -e 'while(<>){if(/pfoldF ..\/..\/(..)/){++$c{$1};++$t}}print$t," ",$t-$c{ua}-$c{au}-$c{cg}-$c{gc}-$c{gu}-$c{ug},"\n"'

# For the main grammar:
# grep "g= 0" Hammerhead_1.results | perl -e 'while(<>){if(/pfoldF ..\/(..)\/../){++$c{$1};++$t}}print$t," ",$t-$c{ua}-$c{au}-$c{cg}-$c{gc}-$c{gu}-$c{ug},"\n"'


use Stockholm;
use Newick;

my $maxcols = 500;  # simulated alignments with more columns than this are junked

my $gram = $ENV{DARTDIR} . "/grammars/pfold.eg";
my $naive = $ENV{DARTDIR} . "/grammars/rev.eg";
my $trials = 1;
my $usage = "Usage: $0 [-g <covariant RNA grammar>] [-n <naive point-substitution grammar>] [-t trials] [-m maxcols] <alignment>\n";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-g") { defined($gram = shift) or die $usage }
    elsif ($arg eq "-n") { defined($naive = shift) or die $usage }
    elsif ($arg eq "-t") { defined($trials = shift) or die $usage }
    elsif ($arg eq "-m") { defined($maxcols = shift) or die $usage }
    else { push @argv, $arg }
}
die $usage unless @argv == 1;
my ($filename) = @argv;

my $stock = Stockholm->from_file ($filename);
my $tree = join "", @{$stock->gf_NH};

my $tmproot = $filename;
$tmproot =~ s/.*?([^\/]+)$/$1/;

my ($tmptree, $tmpsimgram, $tmpsimleaf, $tmprec, $tmpnaive)
    = map ("$tmproot.$_",
	   qw(nh simgram simleaf rec naive));

local *TREE;
open TREE, ">$tmptree";
print TREE $tree;
close TREE;

for (my $trial = 1; $trial <= $trials; ++$trial) {
    warn "Trial $trial of $trials\n";

    my $rand = int (rand (1234567890));
    syswarn ("simgram -rndseed $rand -g $gram -t $tmptree >$tmpsimgram");

    my $simleaf = Stockholm->from_file ($tmpsimgram);
    next if $simleaf->columns > $maxcols;

    my $newick = Newick->from_file ($tmptree);

    my %is_leaf = map (($newick->node_name->[$_] => 1), $newick->leaves);
    $simleaf->seqname ([grep $is_leaf{$_}, @{$simleaf->seqname}]);
    $simleaf->to_file ($tmpsimleaf);

    syswarn ("xrate -g $gram -ar -arpp $tmpsimleaf -log 6 >$tmprec");

    syswarn ("xrate -g $naive -ar -arpp $tmpsimleaf -log 6 >$tmpnaive");

    my $sim = Stockholm->from_file ($tmpsimgram);
    my $rec = Stockholm->from_file ($tmprec);
    my $naif = Stockholm->from_file ($tmpnaive);

    my (%total, %correct);
    for my $node (grep !$is_leaf{$_}, @{$sim->seqname}) {

	# read postprobs from $naif reconstruction
	my %naive_postprob;
	for my $gs (@{$naif->gs_ancrec_CYK_MAP_PP->{$node}}) {
	    if ($gs =~ /State (\S+) columns \((\d+)\).*chars \((.)\).* postprob (\S+)/) {
		my ($state, $col, $char, $postprob) = ($1, $2, $3, $4);
		$naive_postprob{$col} = $postprob;
	    }
	}

	# read postprobs from $rec reconstruction
	for my $gs (@{$rec->gs_ancrec_CYK_MAP_PP->{$node}}) {
	    if ($gs =~ /State (\S+) columns \(([^\)]+)\).*chars \(([^\)]+)\).* postprob (\S+)/) {
		my ($state, $cols, $spc_chars, $postprob) = ($1, $2, $3, $4);
		my @cols = split /\s+/, $cols;
		my @chars = split /\s+/, $spc_chars;
		my (@real_chars, @naive_chars);
		my $naive_postprob = 1;
		for my $col (@cols) {
		    push @real_chars, substr ($sim->seqdata->{$node}, $col-1, 1);
		    push @naive_chars, substr ($naif->gr_ancrec_CYK_MAP->{$node}, $col-1, 1);
		    $naive_postprob *= $naive_postprob{$col};
		}

		my ($chars, $real_chars, $naive_chars)
		    = map (join("",map($_ eq "t" ? "u" : $_, @$_)),
			   \@chars,
			   \@real_chars,
			   \@naive_chars);

		print
		    "$node $state $real_chars/$chars/$naive_chars g= ",
		    $chars eq $real_chars ? 1 : 0,
		    " n= ",
		    $naive_chars eq $real_chars ? 1 : 0,
		    " gpostprob= ", $postprob,
		    " npostprob= ", $naive_postprob,
		    "\n";
	    }
	}
    }
}

sub syswarn {
    my ($cmd) = @_;
    warn "Executing $cmd\n";
    system $cmd;
}
