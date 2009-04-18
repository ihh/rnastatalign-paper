#!/usr/bin/perl -w

use strict;
use Stockholm;

my $ref = Stockholm->from_file (shift);
my $cmp = Stockholm->from_file (shift);

@{$ref->seqname} = ("X", "Y", "Z");
@{$cmp->seqname} = ("X", "Y", "Z");

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

$ref->to_file ("/tmp/ref.stock");
$cmp->to_file ("/tmp/cmp.stock");

my $cmd = "cmpalign.pl /tmp/ref.stock /tmp/cmp.stock >/tmp/cmpalign.out";
system ($cmd) == 0 or die "Couldn't run command: $cmd\n";

open CMPALIGN, "</tmp/cmpalign.out" or die "Couldn't open '/tmp/cmpalign.out'\n";
my ($acc, $sn, $ppv);
while (<CMPALIGN>) {
  if (/^Acc\s+(\S+)$/) { $acc = $1; }
  if (/^Sn\s+(\S+)$/) { $sn = $1; }
  if (/^PPV\s+(\S+)$/) { $ppv = $1; }  
}
close CMPALIGN;

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

$refancstock->to_file ("/tmp/refanc.stock");
$cmpancstock->to_file ("/tmp/cmpanc.stock");

my $stemloc = "$ENV{'DARTDIR'}/bin/stemloc";
$cmd = "$stemloc /tmp/refanc.stock /tmp/cmpanc.stock >/tmp/stemloc.out";
system ($cmd) == 0 or die "Couldn't run '$cmd'.\n";

my $stemlocout = Stockholm->from_file ("/tmp/stemloc.out");
my $ancoverlap = 0;
# catch case of inferred ancestral structure having no base-pairs
if ($stemlocout->gc->{"SS_cons"} =~ /[<>()]/) {
  $ancoverlap = $stemlocout->gc->{"SS_cons"} =~ tr/><//;
  $ancoverlap /= $refanc =~ tr/><//;
}

# print "## ", join ("\t", ("Acc", "Sn", "PPV", "ancestral_bp_overlap")), "\n";
print join ("\t", ($acc, $sn, $ppv, $ancoverlap)), "\n";
