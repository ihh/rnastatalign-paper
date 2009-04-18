#!/usr/bin/perl -w

use Stockholm::Database;
use Newick;

my ($dbfile) = @ARGV;
defined($dbfile) or die "Usage: $0 <Stockholm database>\n";

my $db = Stockholm::Database->from_file ($dbfile);

my $binsize = .1;
my @bin;

for my $stock (@$db) {
    my $tree = Newick->from_string (@{$stock->gf_NH});
    my $total_len = 0;
    map ($total_len += $_, grep (defined, @{$tree->branch_length}));
    my $bin = int ($total_len / $binsize + .5);

    $bin[$bin] = Stockholm::Database->new unless defined $bin[$bin];
    push @{$bin[$bin]}, $stock;

    warn "total_len=$total_len bin=$bin entries=",@{$bin[$bin]}+0,"\n";
}

my $suffix = $dbfile;
$suffix =~ s/^.*?([^\/]+)$/$1/;

for my $bin (0..$#bin) {
    if (defined $bin[$bin]) {
	my $outfile = "bin$bin.$suffix";
	$bin[$bin]->to_file ($outfile);
    }
}
