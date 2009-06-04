#!/usr/bin/perl -w

# imports
use Stockholm;
use FindBin qw($Bin);

# read command-line args
my ($true_alignment) = @ARGV;

# create temp filenames
my $temp_alignment = "$true_alignment.handalign-input";

# define constants
my $tiny = .00001;  # default length for dummy branch
my $handalign = "handalign";  # path to handalign binary
my $loopmodel = "$Bin/loopmodel.eg";  # substitution model for handalign (should be equivalent to loop-submodel of pfold model used by evolsayer.pl to generate alignments)

# load alignment
my $true = Stockholm->from_file ($true_alignment);
warn "True alignment\n";
warn $true->to_string;

# read stuff from alignment
my %seq = %{$true->seqdata};
my $treestr = join "", @{$true->gf_NH};
my $args = join "", @{$true->gf_PARAM};

# extract branch lengths & sequences
my (%branchlen, %seqlen);
for my $branch (qw(X Y Z)) {
    if ($treestr =~ /$branch:([0-9\.eE\+\-]+)/) {
	$branchlen{$branch} = $1;
    } else {
	die "Couldn't find $branch in $treestr";
    }
    die "Couldn't find row $branch in alignment" unless exists $seq{$branch};
    $seq{$branch} =~ s/[\.\-_]//g;
    $seqlen{$branch} = length $seq{$branch};
}

# create rows of new alignment
my %newrow;
$newrow{'Z'} = $seq{'Z'} . ("-" x ($seqlen{'Y'} + $seqlen{'X'}));
$newrow{'Y'} = ("-" x $seqlen{'Z'}) . $seq{'Y'} . ("-" x $seqlen{'X'});
$newrow{'X'} = ("-" x ($seqlen{'Z'} + $seqlen{'Y'})) . $seq{'X'};

$newrow{'xyzroot'} = $newrow{'Z'};
$newrow{'xyzroot'} =~ s/[^\.\-_]/*/g;

$newrow{'xyroot'} = $branchlen{'X'} < $branchlen{'Y'} ? $newrow{'X'} : $newrow{'Y'};
$newrow{'xyroot'} =~ s/[^\.\-_]/*/g;

# create new alignment
my $new = Stockholm->new;

$new->seqname ([qw(xyzroot xyroot X Y Z)]);
$new->seqdata (\%newrow);

$new->add_gf ("NH", "((X:$branchlen{X},Y:$branchlen{Y})xyroot:$branchlen{Z},Z:$tiny)xyzroot;");

# extract loop & stem TKF parameters from alignment
my ($lambda_loop, $mu_loop, $lambda_stem, $mu_stem, $stem_prob);
for my $var (qw(lambda_loop mu_loop lambda_stem mu_stem stem_prob)) {
    if ($args =~ /$var=(\S+)/) { eval('$'."$var=$1") } else { die "Couldn't find $var" }
}
my $loop_len = $lambda_loop / ($mu_loop - $lambda_loop);
my $stem_len = $lambda_stem / ($mu_stem - $lambda_stem);
my $loop_seq_len = ($loop_len * (1 - $stem_prob) + 2 * $stem_len * $stem_prob) / (1 - $stem_prob * $loop_len);
my $stem_seq_len = 2 * $stem_len + $loop_seq_len;
my $gap_len = $stem_prob + $stem_seq_len * (1 - $stem_prob);

# dump alignment to file, and to stderr for debugging purposes
warn "Writing the following alignment to $temp_alignment\n";
warn $new->to_string;

local *TEMP;
open TEMP, ">$temp_alignment";
print TEMP $new->to_string;
close TEMP or die "Couldn't write $temp_alignment: $!";

# run handalign
#my $command = "$handalign -l $loop_seq_len -d $mu_loop -g $gap_len -m $loopmodel -s 0 -r -ub $temp_alignment -ha -hca '-lstdc++'";   # let's not do the hmmoc cache for now; avoid creepy side-effects
my $command = "$handalign -l $loop_seq_len -d $mu_loop -g $gap_len -m $loopmodel -s 0 -r -ub $temp_alignment -ha -hca '-lstdc++' -hc /tmp/indiegram-benchmark-cache";   # in a hurry so put cache back (6/4/09)
warn "Running $command\n";
my @output = `$command`;
#print grep (!/root/, @output);
print @output;
