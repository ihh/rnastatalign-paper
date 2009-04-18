#!/usr/bin/perl -w

use strict;
use Stockholm;

# params: chosen to agree with those hard-coded into indiegram/tkfst.cc
my $looplen = 5;
my $stemlen = 7/3.;
my $pstem = 0.1;

# minimum number of stems in ancestral structure
my $ancstems_min = 2;
my $looplength_max = 10;
my $looplength_min = 3;
my $stemlength_min = 0;
my $seqlength_min = 30;
my $seqlength_max = 70;

# programs
my $evolsayer = "$ENV{'DARTDIR'}/perl/evolsayer.pl";
my $evolsayer_args = "-sub $ENV{'DARTDIR'}/grammars/pfold.eg -sroot -looplen $looplen -stemlen $stemlen -pstem $pstem";

# tree
my ($x, $y) = (1.0, 1.0);
my ($z_min, $z_max) = (0, 2.5);
my $num_bins = 25;
my $z_increment = ($z_max - $z_min) / $num_bins;

my $num_per_bin = 25;



# output directory
my $outdir = shift;
die "Output directory not defined.\n" unless defined $outdir;

# create simulations
my $maxtries = 50000;
my $bin = 1;
for (my $z = $z_min; $z < $z_max; $z += $z_increment) {

  # create directory to hold files for this outgroup branch length
  my $cmd = "mkdir -p $outdir/$z";
  if (! -d "$outdir/$z/") {
    system ($cmd) == 0 or die "Couldn't create directory '$outdir/$z'.\n";
  }

  my $tries = 0;
  my $cnt = 1;
  while ($cnt <= $num_per_bin) {

    # create simulated alignment
    if (! -r "$outdir/$z/$cnt.reference.stock") {

      if ($tries++ > $maxtries) {
	die "Exceeded maximum number of tries.\n";
      }
      
      my $tree = get_tree ($x, $y, $z);
      open TREE, ">/tmp/evolsayer_tree.newick";
      print TREE "$tree\n";
      close TREE;
      
      # trick from http://perldoc.perl.org/functions/srand.html to get a constantly-changing seed
      my $rndseed = time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`;
      my $cmd = "cd $ENV{'DARTDIR'}/perl && $evolsayer -seed $rndseed $evolsayer_args /tmp/evolsayer_tree.newick >/tmp/evolsayer_out.stock 2>/tmp/evolsayer_out.err";
      system ($cmd) == 0 or die "Couldn't run evolsayer.pl as: $cmd.\n";
      
      # get simulated alignment
      my $stock = Stockholm->from_file ("/tmp/evolsayer_out.stock");
      
      # require that there be at least $ancstems_min ancestral stems
      my $ancstems = $stock->get_gf ("ANCSTEMS");
      if ($ancstems =~ /([0-9]+)\s*\(nonempty\)/) { $ancstems = $1; }
      else { die; }
      if ($ancstems < $ancstems_min) { next; }
      
      # requirements on the lengths of loop sequences
      my $looplengths = $stock->get_gf ("ANCLOOPLEN");
      my @looplengths = split / /, $looplengths;
      if (!@looplengths) { next; }
      if (max (@looplengths) > $looplength_max) { next; }
      if (min (@looplengths) < $looplength_min) { next; }
      
      # requirements on the lengths of stem sequences
      my $stemlengths = $stock->get_gf ("ANCSTEMLEN");
      my @stemlengths = split / /, $stemlengths;
      if (!@stemlengths) { next; }
      if (min (@stemlengths) < $stemlength_min) { next; }
      
      # make sure that the longest sequence in the simulated alignment falls in [$seqlength_min, $seqlength_max]
      my $seqlen = 0;
      map { $seqlen = max ($seqlen, ungapped_length ($_)) } values %{$stock->seqdata};
      if ($seqlen < $seqlength_min || $seqlen > $seqlength_max) { next; }
      
      # if the simulated alignment has passed all of these tests, then store it
      $cmd = "$ENV{'DARTDIR'}/bin/gc2gr-ss /tmp/evolsayer_out.stock >$outdir/$z/$cnt.reference.stock";
      system ($cmd) == 0 or die "Couldn't write file '$outdir/$z/$cnt.reference.stock'.";
      warn "Created '$outdir/$z/$cnt.reference.stock' ($cnt / $num_per_bin) for bin z = $z ($bin / $num_bins).\n";
            
    } else {
      warn "Found '$outdir/$z/$cnt.reference.stock.\n";
    }

    # now create input file for indiegram
    if (! -r "$outdir/$z/$cnt.input.stock") {
      my $stock = Stockholm->from_file ("$outdir/$z/$cnt.reference.stock");
      @{$stock->seqname} = ("X", "Y", "Z");
      foreach my $seqname (@{$stock->seqname}) {
	my $newgr = "";
	my $ss = $stock->gr->{"SS"}->{$seqname};
	if (!defined $ss) { die "Couldn't extract SS line.\n"; }
	for (my $col = 0; $col < $stock->columns; ++$col) {
	  if (!$stock->is_gap ($seqname, $col)) {
	    $newgr .= substr ($ss, $col, 1);
	  }
	}
	$newgr =~ s/-/./g;
	$stock->gr->{"SS"}->{$seqname} = $newgr;
	$stock->seqdata->{$seqname} =~ s/-//g;
      }
      my $sscons = $stock->gc->{"SS_cons"};
      delete $stock->gc->{"SS_cons"};
      $stock->to_file ("$outdir/$z/$cnt.input.stock", 100);
      warn "Created '$outdir/$z/$cnt.reference.stock' ($cnt / $num_per_bin) for bin z = $z ($bin / $num_bins).\n";
      warn "       $sscons\n"; # show SS_cons line

    } else {
      warn "Found '$outdir/$z/$cnt.input.stock.\n";
      ## increment count of alignments for this bin
      ++$cnt;
    }
 
  }
  
  ++$bin;

}




sub num_files {
  my $prefix = shift;

  my @list = glob ("$prefix*.stock");
  if (!@list) {
    return 0;
  }

  return scalar @list;

};


sub get_tree {
  my ($x, $y, $z) = @_;

  my $tree = "(X:$x,Y:$y,Z:$z)root;";

  return $tree;

}

sub max {
  my (@F) = @_;
  my $f = $F[0];
  map { if ($f < $_) { $f = $_; } } @F;
  return $f;
}

sub min {
  my (@F) = @_;
  my $f = $F[0];
  map { if ($f > $_) { $f = $_; } } @F;
  return $f;
}

sub ungapped_length {
  my ($data) = @_;

  my $numgaps = $data =~ tr/-._//;
  return (length ($data) - $numgaps);

}
