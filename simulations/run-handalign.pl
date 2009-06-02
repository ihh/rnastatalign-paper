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
$newrow{'X'} = $seq{'X'} . ("-" x ($seqlen{'Y'} + $seqlen{'Z'}));
$newrow{'Y'} = ("-" x $seqlen{'X'}) . $seq{'Y'} . ("-" x $seqlen{'Z'});
$newrow{'Z'} = ("-" x ($seqlen{'X'} + $seqlen{'Y'})) . $seq{'Z'};

$newrow{'subroot'} = $newrow{'Z'};
$newrow{'subroot'} =~ s/[^\.\-_]/*/g;

$newrow{'root'} = $branchlen{'X'} < $branchlen{'Y'} ? $newrow{'X'} : $newrow{'Y'};
$newrow{'root'} =~ s/[^\.\-_]/*/g;

# create new alignment
my $new = Stockholm->new;

$new->seqname ([qw(root subroot X Y Z)]);
$new->seqdata (\%newrow);

$new->add_gf ("NH", "((X:$branchlen{X},Y:$branchlen{Y})root:$branchlen{Z},Z:$tiny)subroot;");

# extract loop TKF parameters from alignment
my ($lambda, $mu);
if ($args =~ /lambda_loop=(\S+)/) { $lambda = $1 } else { die "Couldn't find lambda" }
if ($args =~ /mu_loop=(\S+)/) { $mu = $1 } else { die "Couldn't find mu" }
my $loop_len = $lambda / ($mu - $lambda);

# dump alignment to file, and to stderr for debugging purposes
warn "Writing the following alignment to $temp_alignment\n";
warn $new->to_string;

local *TEMP;
open TEMP, ">$temp_alignment";
print TEMP $new->to_string;
close TEMP or die "Couldn't write $temp_alignment: $!";

# run handalign
my $command = "$handalign -tkf91 -l $loop_len -d $mu -m $loopmodel -s 0 -r -ub $temp_alignment -ha -hca '-lstdc++' -hc /tmp/indiegram-benchmark-cache";
warn "Running $command\n";
my @output = `$command`;
#print grep (!/root/, @output);
print @output;