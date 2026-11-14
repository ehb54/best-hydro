#!/usr/bin/perl
# check_best_frame.pl a.be a.log b.be
#   Checks results that must not depend on BEST's coordinate frame.
#   a: a mesh.  b: the same mesh rotated by test/permute_mesh.pl, (x,y,z) -> (y,z,x).
#   1. Krt at the center of resistance (a.log) is symmetric, to 1e-3 of its
#      largest element.  It is symmetric there by definition.
#   2. b's center vectors are a's, rotated, to 1e-3 A.  b's 1/3 traces and
#      Dtt eigenvalues are a's, to 1e-4 (one unit in the fifth printed digit).
#      Drr eigenvalues only to 1e-3: the BE Krr is not exactly symmetric, and
#      JACOBI reads only the upper triangle, so they move by about 2e-4 when
#      the mesh is rotated.
#   Exits nonzero if a check fails.
use strict;
use warnings;

my ( $abe, $alog, $bbe ) = @ARGV;
die "usage: $0 a.be a.log b.be\n" unless defined $bbe;

sub lines { open my $h, '<', $_[0] or die "$_[0]: $!\n"; my @l = <$h>; @l }
sub after {    # the numbers on the n lines after the line matching re
   my ( $l, $re, $n, $f ) = @_;
   for my $i ( 0 .. $#$l ) {
      next unless $l->[$i] =~ $re;
      return [ map { [ split ' ', $l->[ $i + $_ ] ] } 1 .. $n ];
   }
   die "$f: no line matching $re\n";
}
sub value {
   my ( $l, $label, $f ) = @_;
   for ( @$l ) { return $1 if /\Q$label\E\s*=?\s*([-+]?[0-9.]+(?:[EeDd][-+]?[0-9]+)?)/ }
   die "$f: no $label\n";
}
sub parse {
   my $f = shift;
   my @l = lines( $f );
   return {
      urc => after( \@l, qr/Vector between Center of Resistance\/Centroid/, 1, $f )->[0],
      udr => after( \@l, qr/Vector between Centers of Diffusion\/Resistance/, 1, $f )->[0],
      dtt => [ sort { $a <=> $b } @{ after( \@l, qr/EIGENVALUES OF Dtt TENSOR \(1\/A\)/, 1, $f )->[0] } ],
      drr => [ sort { $a <=> $b } @{ after( \@l, qr/EIGENVALUES OF Drr TENSOR \(1\/A\^3\)/, 1, $f )->[0] } ],
      dt  => value( \@l, '1/3 Trace DTT(EXACT)', $f ),
      dr  => value( \@l, '1/3 Trace DRR(EXACT)', $f ),
   };
}

my $bad = 0;
my @l   = lines( $alog );
my $k   = after( \@l, qr/Symmetric Krt Tensor/, 3, $alog );
my ( $max, $asym ) = ( 0, 0 );
for my $i ( 0 .. 2 ) {
   for my $j ( 0 .. 2 ) {
      $max  = abs $k->[$i][$j] if abs $k->[$i][$j] > $max;
      $asym = abs( $k->[$i][$j] - $k->[$j][$i] ) if abs( $k->[$i][$j] - $k->[$j][$i] ) > $asym;
   }
}
if ( $asym > 1e-3 * $max ) {
   printf "  Krt at the center of resistance is not symmetric: max |Krt - Krt^T| %.3g, max |Krt| %.3g\n", $asym, $max;
   $bad++;
}

my ( $a, $b ) = ( parse( $abe ), parse( $bbe ) );
for my $v ( [ urc => 'center of resistance - centroid' ], [ udr => 'center of diffusion - center of resistance' ] ) {
   my ( $x, $y ) = ( $a->{ $v->[0] }, $b->{ $v->[0] } );
   my @rot = @$x[ 1, 2, 0 ];
   for ( 0 .. 2 ) {
      next if abs( $y->[$_] - $rot[$_] ) <= 1e-3;
      printf "  %s: rotated mesh gives (%s), expected (%s)\n", $v->[1], join( ' ', @$y ), join( ' ', @rot );
      $bad++;
      last;
   }
}
for my $q ( [ dt => '1/3 Trace DTT', 1e-4 ], [ dr => '1/3 Trace DRR', 1e-4 ], [ dtt => 'Dtt eigenvalues', 1e-4 ],
   [ drr => 'Drr eigenvalues', 1e-3 ] )
{
   my @x = ref $a->{ $q->[0] } ? @{ $a->{ $q->[0] } } : $a->{ $q->[0] };
   my @y = ref $b->{ $q->[0] } ? @{ $b->{ $q->[0] } } : $b->{ $q->[0] };
   for ( 0 .. $#x ) {
      next if abs( $y[$_] - $x[$_] ) <= $q->[2] * abs( $x[$_] );
      printf "  %s: rotated mesh gives %s, expected %s\n", $q->[1], join( ' ', @y ), join( ' ', @x );
      $bad++;
      last;
   }
}
exit( $bad ? 1 : 0 );
