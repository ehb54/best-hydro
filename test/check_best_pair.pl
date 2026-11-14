#!/usr/bin/perl
# check_best_pair.pl a.be b.be [tolerance]
#   Checks that two BEST runs on the same mesh give the same D_t, D_r, CHI and
#   [eta] (whichever a.be has): each value in b.be must be within the relative
#   tolerance (default 2e-5) of a.be's, or within one unit in the last digit
#   printed.  For the H-matrix solver against the dense one.  Exits nonzero if
#   a value is missing or differs more.
use strict;
use warnings;

my ( $fa, $fb, $tol ) = @ARGV;
die "usage: $0 a.be b.be [tolerance]\n" unless defined $fb;
$tol //= 2e-5;
my @labels = ( '1/3 Trace DTT(EXACT)', '1/3 Trace DRR(EXACT)', 'CHI', 'ETA (cm^3/g)' );
my $num = qr/[-+]?([0-9]*)\.([0-9]+)(?:[EeDd]([-+]?[0-9]+))?/;

sub slurp { my $f = shift; local $/; open my $h, '<', $f or die "$f: $!\n"; return <$h> }
my ( $ta, $tb ) = ( slurp( $fa ), slurp( $fb ) );
my ( $bad, $seen ) = ( 0, 0 );
for my $label ( @labels ) {
   my ( $va ) = $ta =~ /^\s*\Q$label\E\s*=?\s*($num)/m or next;
   $seen++;
   my ( $vb ) = $tb =~ /^\s*\Q$label\E\s*=?\s*($num)/m;
   if ( !defined $vb ) {
      print "  $label: not found in $fb\n";
      $bad++;
      next;
   }
   # One unit in the last digit of a, as printed.
   my ( undef, $frac, $exp ) = $va =~ /^$num$/;
   my $ulp = 10**( ( $exp // 0 ) - length( $frac ) );
   ( my $x = $va ) =~ tr/Dd/EE/;
   ( my $y = $vb ) =~ tr/Dd/EE/;
   my $d = abs( $x - $y );
   if ( $d > $tol * abs( $x ) && $d > 1.01 * $ulp ) {
      printf "  %s: %s against %s (relative difference %.2g > %g)\n", $label, $vb, $va, $d / abs( $x ), $tol;
      $bad++;
   }
}
if ( !$seen ) {
   print "  none of the values found in $fa\n";
   $bad++;
}
exit( $bad ? 1 : 0 );
