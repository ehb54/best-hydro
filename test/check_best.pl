#!/usr/bin/perl
# check_best.pl file.be expected.txt
#   Compares values in a BEST .be file with expected values.  Each line of
#   expected.txt is "value relative-tolerance label"; the label is the text in
#   the .be file that the value follows ("label = value" or "label value").
#   Exits nonzero if a value is missing or out of tolerance.
use strict;
use warnings;

my ( $be, $exp ) = @ARGV;
die "usage: $0 file.be expected.txt\n" unless defined $exp;
my $text = do { local $/; open my $h, '<', $be or die "$be: $!\n"; <$h> };
open my $h, '<', $exp or die "$exp: $!\n";
my $bad = 0;
while ( <$h> ) {
   next if /^\s*(#|$)/;
   my ( $want, $tol, $label ) = /^\s*(\S+)\s+(\S+)\s+(.*?)\s*$/ or die "$exp: bad line: $_";
   my ( $got ) = $text =~ /\Q$label\E\s*=?\s*([-+]?[0-9.]+(?:[EeDd][-+]?[0-9]+)?)/;
   if ( !defined $got ) {
      print "  $label: not found in $be\n";
      $bad++;
      next;
   }
   ( my $g = $got ) =~ tr/Dd/EE/;
   my $rel = abs( $g - $want ) / abs( $want );
   if ( $rel > $tol ) {
      printf "  %s: got %s, expected %s (relative difference %.2g > %g)\n", $label, $got, $want, $rel, $tol;
      $bad++;
   }
}
exit( $bad ? 1 : 0 );
