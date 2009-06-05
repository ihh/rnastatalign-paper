#!/usr/bin/perl
$tkfst = $ENV{'DARTDIR'} . "/src/indiegram/tkfst.cc";
$parse = "alignments/0.1/13.reference.parse";
open IN, "<$parse" or die;
open OUT, ">$parse.inits";
print OUT "// $parse.inits: generated by $0 from $parse and $tkfst\n\n";
while(<IN>){next unless /^\s*[^;\s]/;chomp;if(/\(?(\S+).*\->\s+(\S+)\s+(\S+)/){mysys("egrep \"init_bifurc.*$1.*$2.*$3\" $tkfst");$x=undef}else{s/^\s*\(?([^\(\)\s]+).*/$1/; if ( defined$x ) {($a,$b)=/end/?($_,$x):($x,$_);mysys("egrep '^[^M]*transition_scores\[^M]*$a\[^M]*$b' $tkfst")}$x=$_ eq"end"?undef:$_}}
sub mysys { my$c=shift; # print OUT "// $c\n",`$c`;
print OUT `$c`}
close IN;
close OUT or die;
