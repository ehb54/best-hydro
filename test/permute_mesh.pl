#!/usr/bin/perl
# permute_mesh.pl in out
#   Rotates a BEST triangulation by the cyclic axis permutation (x,y,z) -> (y,z,x):
#   new x = old y, new y = old z, new z = old x.  The number strings are moved,
#   not reformatted, so the rotation is exact.
use strict;
use warnings;

my ( $in, $out ) = @ARGV;
die "usage: $0 in out\n" unless defined $out;
open my $h, '<', $in  or die "$in: $!\n";
open my $o, '>', $out or die "$out: $!\n";
print $o scalar <$h>;
while ( <$h> ) {
   my @t = split ' ';
   die "$in line $.: expected 9 numbers, got " . scalar( @t ) . "\n" unless @t == 9;
   print $o join( ' ', '', map { @t[ 3 * $_ + 1, 3 * $_ + 2, 3 * $_ ] } 0 .. 2 ), "\n";
}
