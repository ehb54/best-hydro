#!/usr/bin/perl
# meshcheck.pl mesh...
#   Checks that each BEST input mesh (as written by rcoal) is a closed,
#   manifold surface: every edge is shared by exactly two triangles, and no
#   triangle appears twice.  Vertices are matched by their coordinates.
#   Exits nonzero if any mesh fails.
use strict;
use warnings;

die "usage: $0 mesh...\n" unless @ARGV;
my $failed = 0;
for my $f ( @ARGV ) {
   open my $h, '<', $f or die "$f: $!\n";
   my $header = <$h>;
   my ( %vid, %edge, %tri );
   my ( $nv, $nt, $dup ) = ( 0, 0, 0 );
   while ( <$h> ) {
      my @c = split;
      next unless @c == 9;
      $nt++;
      my @t;
      for my $k ( 0 .. 2 ) {
         my $key = join ' ', map { sprintf '%.6f', $_ } @c[ 3 * $k .. 3 * $k + 2 ];
         $vid{ $key } //= $nv++;
         push @t, $vid{ $key };
      }
      $dup++ if $tri{ join ' ', sort { $a <=> $b } @t }++;
      for my $e ( [ 0, 1 ], [ 1, 2 ], [ 2, 0 ] ) {
         my ( $p, $q ) = sort { $a <=> $b } @t[ @$e ];
         $edge{ "$p $q" }++;
      }
   }
   my $open   = grep { $_ == 1 } values %edge;
   my $pinch  = grep { $_ > 2 } values %edge;
   if ( $open || $pinch || $dup ) {
      print "  $f: $nt triangles, $open open edges, $pinch edges shared by more than two triangles, $dup duplicate triangles\n";
      $failed++;
   }
}
exit( $failed ? 1 : 0 );
