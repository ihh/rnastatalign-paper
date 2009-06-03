#!/usr/bin/perl -w

use strict;
use Stockholm;

my $tmpdir = "/tmp";

my $ref = Stockholm->from_file (shift);
my $cmp = Stockholm->from_file (shift);
my $tkf = Stockholm->from_file (shift);
my $stemloc = Stockholm->from_file (shift);
my $stemlocama = Stockholm->from_file (shift);

# indiegram adds cruft to the sequence names and prints non-Stockholm stuff that gets read back incorrectly by Stockholm.pm
# handalign and indiegram both add "root" and "subroot" lines
# so we need to tidy up the sequence name fields...
@{$ref->seqname} = @{$cmp->seqname} = @{$tkf->seqname} = ("X", "Y", "Z");

$ref->gc->{"ancestral_SS"} = $ref->gr->{"SS"}->{"root"};
delete $ref->gr->{"SS"}->{"root"};

while (my ($name, $data) = each %{$cmp->seqdata}) {
  if ($name =~ /([XYZ])\/.+/) {
    delete $cmp->seqdata->{$name};
    $cmp->seqdata->{$1} = $data;
    $cmp->gr->{"SS"}->{$1} = $cmp->gr->{"SS"}->{$name};
    delete $cmp->gr->{"SS"}->{$name};
  }
}

for my $file (qw(ref cmp tkf stemloc stemlocama)) {
    eval('$'.$file)->to_file ("$tmpdir/$file.stock");
}

my ($acc, $sn, $ppv, $tcs) = cmpalign ("$tmpdir/ref.stock", "$tmpdir/cmp.stock");
my @acc_sn_ppv_tcs = map ([cmpalign ("$tmpdir/ref.stock", "$tmpdir/$_.stock")], qw(tkf stemloc stemlocama));

my $refanc = $ref->gc->{"ancestral_SS"};
my $refseq = $ref->seqdata->{"root"};
for (my $i = 0; $i < length $refseq; ++$i) {
  if (substr ($refseq, $i, 1) ne '-') {
    if (substr ($refanc, $i, 1) eq '-') {
      substr ($refanc, $i, 1) = '.';
    }
    if (substr ($refanc, $i, 1) eq '<') {
      substr ($refseq, $i, 1) = "C";
    } elsif (substr ($refanc, $i, 1) eq '>') {
      substr ($refseq, $i, 1) = "G";
    } else  {
      substr ($refseq, $i, 1) = "A";
    }
  }
}
$refanc =~ s/-//g;
$refseq =~ s/-//g;

my $cmpanc = $cmp->gc->{"ancestral_SS"};
$cmpanc =~ s/-//g;
my $cmpseq = 'N' x length $cmpanc;
for (my $i = 0; $i < length $cmpanc; ++$i) {
  if (substr ($cmpanc, $i, 1) eq '<') {
    substr ($cmpseq, $i, 1) = "C";
  } elsif (substr ($cmpanc, $i, 1) eq '>') {
    substr ($cmpseq, $i, 1) = "G";
  } else  {
    substr ($cmpseq, $i, 1) = "A";
  }
}
$cmpseq =~ s/-//g;

my $refancstock = Stockholm->new;
my $cmpancstock = Stockholm->new;
@{$refancstock->seqname} = ("ref");
@{$cmpancstock->seqname} = ("cmp");
$refancstock->seqdata->{"ref"} = $refseq;
$cmpancstock->seqdata->{"cmp"} = $cmpseq;
$refancstock->gr->{"SS"}->{"ref"} = $refanc;
$cmpancstock->gr->{"SS"}->{"cmp"} = $cmpanc;

$refancstock->to_file ("$tmpdir/refanc.stock");
$cmpancstock->to_file ("$tmpdir/cmpanc.stock");

my $stemloc_bin = "$ENV{'DARTDIR'}/bin/stemloc";
my $cmd = "$stemloc_bin $tmpdir/refanc.stock $tmpdir/cmpanc.stock >$tmpdir/stemloc.out";
system ($cmd) == 0 or die "Couldn't run '$cmd'.\n";

my $stemlocout = Stockholm->from_file ("$tmpdir/stemloc.out");
my $ancoverlap = 0;
# catch case of inferred ancestral structure having no base-pairs
if ($stemlocout->gc->{"SS_cons"} =~ /[<>()]/) {
  $ancoverlap = $stemlocout->gc->{"SS_cons"} =~ tr/><//;
  $ancoverlap /= $refanc =~ tr/><//;
}

warn "## ", join ("\t", ("Acc", "Sn", "PPV", "ancestral_bp_overlap", map (("Acc_$_", "Sn_$_", "PPV_$_"), qw(tkf stemloc stemlocama)))), "\n";
print join ("\t", ($acc, $sn, $ppv, $ancoverlap, map (@$_[0..2], @acc_sn_ppv_tcs))), "\n";

sub cmpalign {
    my ($ref, $cmp) = @_;
    my $cmd = "cmpalign -s $ref $cmp";
    warn "Running $cmd\n";
    my $line = `$cmd`;
    warn "Output: $line" if $line =~ /\S/;
    ($acc, $sn, $ppv, $tcs) = split /\s+/, $line;
    return ($acc, $sn, $ppv, $tcs);
}
