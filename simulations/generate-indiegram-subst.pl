#!/usr/bin/perl -w

use strict;
use PhyloGram;

my $g = PhyloGram->from_file ("/Users/rbradley/projects/dart/grammars/pfold.eg");

my $chain = $g->find_chain ("LNUC");

print $chain->terminal->to_string,"\n";
my $mutateHash = $chain->mutate_hash;
my @states = map (join ("",@{$_->state->value}), $chain->find_all("initial"));


# get initial probability distribution
my %init;
foreach my $i (@states) {
  my $p = $chain->initial($i);
  if (defined $p) {
    $init{$i} = $p->prob->value;
  } else {
    warn "Setting initial prob of $i = 0.\n";
    $init{$i} = 0;
  }
}

# get rate matrix
my %rates;
foreach my $i (@states) {
  my $total = 0;
  foreach my $j (@states) {
    my $m = $chain->mutate($i,$j,$mutateHash);
    if (defined $m) {
      $rates{$i}->{$j} = $m->rate->value;
      $total += $rates{$i}->{$j};
    } else {
      warn "Setting $i -> $j = 0.\n";
      $rates{$i}->{$j} = 0;
    }
  }
  $rates{$i}->{$i} = -$total;
}

my @sortedstates = qw(aa ca ga ua ac cc gc uc ag cg gg ug au cu gu uu);

print "bkstem << \"1 16\\n\"\n";
print "<< \"a c g u A C G U b d h v B D H V\\n\"\n";
# initial dist
print "<< \"";
foreach my $i (@sortedstates) {
  print $init{$i}, " ";
}
print "\\n\"\n";
# rate matrix
foreach my $i (@sortedstates) {
  print "<< \"";
  foreach my $j (@sortedstates) {
    print $rates{$i}->{$j}, " ";
  }
  print "\\n\"\n";
}
print "\"0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n0\\n\"\n";
