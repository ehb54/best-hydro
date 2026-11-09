#!/usr/bin/perl
# meshcmp.pl a b [tolerance]
#   Compares two BEST input meshes (as written by rcoal): the same number of
#   triangles, the same exact area and volume on the first line, and vertex
#   coordinates that agree within the tolerance (default 1e-4, in the
#   coordinates' units).  Exits nonzero if they differ.
use strict;
use warnings;

my ( $fa, $fb, $tol ) = @ARGV;
die "usage: $0 a b [tolerance]\n" unless defined $fb;
$tol //= 1e-4;
open my $a, '<', $fa or die "$fa: $!\n";
open my $b, '<', $fb or die "$fb: $!\n";
my @ha = split ' ', scalar <$a>;
my @hb = split ' ', scalar <$b>;
if ( $ha[ 0 ] != $hb[ 0 ] ) {
   print "  $fa has $ha[0] triangles, $fb has $hb[0]\n";
   exit 1;
}
for my $i ( 1, 2 ) {
   next unless defined $ha[ $i ] && defined $hb[ $i ];
   if ( abs( $ha[ $i ] - $hb[ $i ] ) > 1e-6 * abs( $ha[ $i ] ) ) {
      print "  first lines differ: @ha | @hb\n";
      exit 1;
   }
}
my ( $line, $max ) = ( 1, 0 );
while ( defined( my $la = <$a> ) ) {
   my $lb = <$b>;
   $line++;
   if ( !defined $lb ) {
      print "  $fb ends at line $line\n";
      exit 1;
   }
   my @x = split ' ', $la;
   my @y = split ' ', $lb;
   for my $i ( 0 .. $#x ) {
      my $d = abs( $x[ $i ] - $y[ $i ] );
      $max = $d if $d > $max;
   }
}
if ( defined <$b> ) {
   print "  $fb has more lines than $fa\n";
   exit 1;
}
if ( $max > $tol ) {
   printf "  coordinates differ by up to %.3g (tolerance %g)\n", $max, $tol;
   exit 1;
}
exit 0;
