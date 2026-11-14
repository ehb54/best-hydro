// BEST H-matrix solver: the self blocks (tgetgr.f) and off-diagonal blocks (triangleint.f) of GR,
// ported with the Fortran's order of operations.  Build with -ffp-contract=off.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "kernel.h"

#include <cmath>

namespace bh {

namespace {

const double kPi = 3.1415926535897932385;

// 21-point Gauss abscissas and weights on [0, 1] (triangleint.f).
const double kY[ 21 ] = { 3.1239146898052499e-3, 1.6386580716846853e-2, 3.9950332924799586e-2,
                          7.3318317708341358e-2, 0.1157800182621610,    0.1664305979012938,
                          0.2241905820563901,    0.2878289398962806,    0.3559893415987995,
                          0.4272190729195525,    0.5,                   0.5727809270804476,
                          0.6440106584012005,    0.7121710601037194,    0.7758094179436099,
                          0.8335694020987062,    0.8842199817378390,    0.9266816822916586,
                          0.9600496670752004,    0.9836134192831531,    0.9968760853101948 };
const double kW[ 21 ] = { 1.601722825777433e-2, 3.695378977085249e-2, 5.713442542685720e-2,
                          7.610011362837930e-2, 9.344442345603390e-2, 0.1087972991671484,
                          0.1218314160537285,   0.1322689386333375,   0.1398873947910731,
                          0.1445244039899701,   0.1460811336496904,   0.1445244039899701,
                          0.1398873947910731,   0.1322689386333375,   0.1218314160537285,
                          0.1087972991671484,   9.344442345603390e-2, 7.610011362837930e-2,
                          5.713442542685720e-2, 3.695378977085249e-2, 1.601722825777433e-2 };

double si( double x ) { return 2.0 * ( std::log( ( 1.0 + std::sqrt( 1.0 - x * x ) ) / x ) - std::sqrt( 1.0 - x * x ) ); }
double csi( double x ) { return 3.0 - 1.0 / x - 2.0 * x; }
double dot( const double* a, const double* b ) { return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ]; }

}  // namespace

Mesh mesh_from_best( int n, const double* v1, const double* v2, const double* v3, const double* area, double afac ) {
   Mesh m;
   m.n = n;
   m.afac = afac;
   m.v1.resize( 3 * size_t( n ) );
   m.v2.resize( 3 * size_t( n ) );
   m.v3.resize( 3 * size_t( n ) );
   m.area.assign( area, area + n );
   for ( int i = 0; i < n; ++i ) {
      for ( int c = 0; c < 3; ++c ) {
         m.v1[ 3 * size_t( i ) + c ] = v1[ size_t( c ) * n + i ];
         m.v2[ 3 * size_t( i ) + c ] = v2[ size_t( c ) * n + i ];
         m.v3[ 3 * size_t( i ) + c ] = v3[ size_t( c ) * n + i ];
      }
   }
   return m;
}

Kernel::Kernel( const Mesh& mesh, double cg ) : n_( mesh.n ), cg_( cg ), mesh_( mesh ) {
   x_.resize( 3 * size_t( n_ ) );
   self_.resize( 9 * size_t( n_ ) );
   const double afac = mesh.afac;
   const double pi = kPi;
   double g[ 4 ][ 4 ] = {}, rr[ 4 ][ 4 ] = {};
   for ( int n1 = 0; n1 < n_; ++n1 ) {
      const double* v1 = &mesh.v1[ 3 * n1 ];
      const double* v2 = &mesh.v2[ 3 * n1 ];
      const double* v3 = &mesh.v3[ 3 * n1 ];
      // Incenter (tgetgr.f).
      const double v21x = v2[ 0 ] - v1[ 0 ];
      const double v21y = v2[ 1 ] - v1[ 1 ];
      const double v21z = v2[ 2 ] - v1[ 2 ];
      const double v31x = v3[ 0 ] - v1[ 0 ];
      const double v31y = v3[ 1 ] - v1[ 1 ];
      const double v31z = v3[ 2 ] - v1[ 2 ];
      const double v32x = v3[ 0 ] - v2[ 0 ];
      const double v32y = v3[ 1 ] - v2[ 1 ];
      const double v32z = v3[ 2 ] - v2[ 2 ];
      const double a = std::sqrt( v21x * v21x + v21y * v21y + v21z * v21z );
      const double b = std::sqrt( v31x * v31x + v31y * v31y + v31z * v31z );
      const double c = std::sqrt( v32x * v32x + v32y * v32y + v32z * v32z );
      const double cosC = ( v21x * v31x + v21y * v31y + v21z * v31z ) / a / b;
      const double temp = ( a + b - c ) / ( 2. * ( 1. + cosC ) );
      x_[ 3 * n1 ] = v1[ 0 ] + temp * ( v21x / a + v31x / b );
      x_[ 3 * n1 + 1 ] = v1[ 1 ] + temp * ( v21y / a + v31y / b );
      x_[ 3 * n1 + 2 ] = v1[ 2 ] + temp * ( v21z / a + v31z / b );
      // Normal.
      double xnorm = v31y * v21z - v31z * v21y;
      double ynorm = v31z * v21x - v31x * v21z;
      double znorm = v31x * v21y - v31y * v21x;
      const double area2 = std::sqrt( xnorm * xnorm + ynorm * ynorm + znorm * znorm ) / 2.;
      xnorm = .5 * xnorm / area2;
      ynorm = .5 * ynorm / area2;
      znorm = .5 * znorm / area2;
      // Exact self integrals.
      const double cc = std::acos( cosC );
      const double r = temp * std::sin( cc );
      const double cosB = -( v21x * v32x + v21y * v32y + v21z * v32z ) / a / c;
      const double bb = std::acos( cosB );
      const double aa = pi - bb - cc;
      const double temp2 = 1. / ( 2. * std::cos( 0.5 * cc ) );
      const double dc = ( a + b - c ) * temp2;
      const double db = std::sqrt( a * a + ( b - c ) * ( b - c ) - 2. * a * ( b - c ) * cosC ) * temp2;
      const double da = std::sqrt( b * b + ( a - c ) * ( a - c ) - 2. * b * ( a - c ) * cosC ) * temp2;
      const double ha = r / da;
      const double hb = r / db;
      const double hc = r / dc;
      const double dpa =
          ( pi - aa ) * ( da - r ) - 2.0 * ( da * std::acos( ha ) - r * std::log( ( 1. + std::sqrt( 1. - ha * ha ) ) / ha ) );
      const double dpb =
          ( pi - bb ) * ( db - r ) - 2.0 * ( db * std::acos( hb ) - r * std::log( ( 1. + std::sqrt( 1. - hb * hb ) ) / hb ) );
      const double dpc =
          ( pi - cc ) * ( dc - r ) - 2.0 * ( dc * std::acos( hc ) - r * std::log( ( 1. + std::sqrt( 1. - hc * hc ) ) / hc ) );
      const double sca = 0.5 * r *
                         ( ( std::cos( 2.0 * cc ) + std::cos( 2.0 * bb ) ) * si( ha ) -
                           ( std::sin( 2.0 * cc ) + std::sin( 2.0 * bb ) ) * csi( ha ) );
      const double scb = 0.5 * r * ( ( 1. + std::cos( 2.0 * bb ) ) * si( hb ) + std::sin( 2.0 * bb ) * csi( hb ) );
      const double scc = 0.5 * r * ( ( 1. + std::cos( 2.0 * cc ) ) * si( hc ) + std::sin( 2.0 * cc ) * csi( hc ) );
      const double cca = -std::sin( bb - cc ) * r * ( std::cos( aa ) * si( ha ) + std::sin( aa ) * csi( ha ) );
      const double ccb = std::sin( bb ) * r * ( std::cos( bb ) * si( hb ) + std::sin( bb ) * csi( hb ) );
      const double ccc = -std::sin( cc ) * r * ( std::cos( cc ) * si( hc ) + std::sin( cc ) * csi( hc ) );
      const double cor = std::sqrt( afac );
      const double cg = cg_;
      g[ 1 ][ 1 ] =
          cor * ( ( 2.0 + cg ) * pi * r + ( 1.0 + 0.5 * cg ) * cg * ( dpa + dpb + dpc ) + cg * 0.25 * ( sca + scb + scc ) );
      g[ 1 ][ 2 ] = 0.5 * cor * ( cca + ccb + ccc ) * cg;
      g[ 2 ][ 1 ] = g[ 1 ][ 2 ];
      g[ 2 ][ 2 ] = g[ 1 ][ 1 ] - 0.5 * cor * cg * ( sca + scb + scc );
      g[ 3 ][ 3 ] = cor * ( 2.0 * pi * r + ( dpa + dpb + dpc ) );
      g[ 1 ][ 3 ] = g[ 2 ][ 3 ] = g[ 3 ][ 1 ] = g[ 3 ][ 2 ] = 0.0;
      // Rotate to the lab frame.
      rr[ 1 ][ 1 ] = v21x / a;
      rr[ 1 ][ 2 ] = v21y / a;
      rr[ 1 ][ 3 ] = v21z / a;
      rr[ 3 ][ 1 ] = xnorm;
      rr[ 3 ][ 2 ] = ynorm;
      rr[ 3 ][ 3 ] = znorm;
      rr[ 2 ][ 1 ] = ( ynorm * v21z - znorm * v21y ) / a;
      rr[ 2 ][ 2 ] = ( znorm * v21x - xnorm * v21z ) / a;
      rr[ 2 ][ 3 ] = ( xnorm * v21y - ynorm * v21x ) / a;
      for ( int i = 1; i <= 3; ++i ) {
         for ( int j = 1; j <= 3; ++j ) {
            double t = 0.0;
            for ( int k = 1; k <= 3; ++k ) {
               t = t + rr[ k ][ i ] * ( g[ k ][ 1 ] * rr[ 1 ][ j ] + g[ k ][ 2 ] * rr[ 2 ][ j ] + g[ k ][ 3 ] * rr[ 3 ][ j ] );
            }
            self_[ 9 * size_t( n1 ) + 3 * ( i - 1 ) + ( j - 1 ) ] = t;
         }
      }
   }
}

void Kernel::pair( int n1, int n2, double* g ) const {
   const double* p1 = &mesh_.v1[ 3 * size_t( n2 ) ];
   const double* p2 = &mesh_.v2[ 3 * size_t( n2 ) ];
   const double* p3 = &mesh_.v3[ 3 * size_t( n2 ) ];
   const double* xc = &x_[ 3 * size_t( n1 ) ];
   double v21[ 3 ], v31[ 3 ], v32[ 3 ], v1c[ 3 ], v2c[ 3 ], v3c[ 3 ];
   for ( int c = 0; c < 3; ++c ) {
      v21[ c ] = p2[ c ] - p1[ c ];
      v31[ c ] = p3[ c ] - p1[ c ];
      v32[ c ] = p3[ c ] - p2[ c ];
      v1c[ c ] = p1[ c ] - xc[ c ];
      v2c[ c ] = p2[ c ] - xc[ c ];
      v3c[ c ] = p3[ c ] - xc[ c ];
   }
   const double d2121 = dot( v21, v21 );
   const double d2131 = dot( v21, v31 );
   const double d2132 = dot( v21, v32 );
   const double d211c = dot( v21, v1c );
   const double d212c = dot( v21, v2c );
   const double d3131 = dot( v31, v31 );
   const double d3232 = dot( v32, v32 );
   const double d3231 = dot( v32, v31 );
   const double d322c = dot( v32, v2c );
   const double d1c1c = dot( v1c, v1c );
   const double d1c31 = dot( v1c, v31 );
   const double d2c2c = dot( v2c, v2c );
   const double d2c1c = dot( v2c, v1c );
   const double d3c31 = dot( v3c, v31 );
   const double d3c32 = dot( v3c, v32 );
   const double d3c3c = dot( v3c, v3c );

   // Square roots that the Fortran functions h, hy, hxy compute with identical expressions.
   const double s21 = std::sqrt( d2121 );
   double t0 = 0.0, s0 = 0.0, sx = 0.0, sy = 0.0, sxy = 0.0, syy = 0.0, sxx = 0.0;
   for ( int j = 0; j < 21; ++j ) {
      const double y = kY[ j ];
      const double w = kW[ j ];
      const double q1 = std::sqrt( d1c1c + 2.0 * y * d1c31 + y * y * d3131 );
      const double q2 = std::sqrt( d2c2c + 2.0 * y * d322c + y * y * d3232 );
      // h(y)
      const double h = std::log( ( d212c + d2132 * y + s21 * q2 ) /
                                 ( d211c + d2131 * y + s21 * std::sqrt( d1c1c + 2 * d1c31 * y + d3131 * y * y ) ) ) /
                       s21;
      // hy(y)
      double fa = ( d211c + y * d2131 ) / q1;
      double fb = ( y * d2132 + d212c ) / q2;
      double fc = ( d211c * d211c - d1c1c * d2121 - 2.0 * y * d1c31 * d2121 + 2.0 * y * d211c * d2131 +
                    y * y * d2131 * d2131 - y * y * d2121 * d3131 );
      const double hy = ( fa - fb ) / fc;
      // hxy(y)
      fa = q1;
      fb = ( y * d211c - d2c1c - 2.0 * y * d1c31 - y * d2131 - y * y * d3231 ) / q2;
      fc = ( -d211c * d211c + d1c1c * d2121 + 2.0 * y * d1c31 * d2121 - 2.0 * y * d211c * d2131 -
             y * y * ( d2131 * d2131 ) + y * y * d2121 * d3131 );
      const double hxy = ( fa + fb ) / fc;
      // hx(y)
      fa = ( -d1c31 - y * d2131 ) / std::sqrt( d1c1c + 2.0 * y * d211c + y * y * d2121 );
      fb = ( -y * d3231 + d3c31 ) / std::sqrt( y * y * d3232 - 2.0 * y * d3c32 + d3c3c );
      fc = ( -d1c31 * d1c31 - 2.0 * y * d1c31 * d2131 - y * y * d2131 * d2131 + d1c1c * d3131 +
             2.0 * y * d211c * d3131 + y * y * d2121 * d3131 );
      const double hx = ( fa + fb ) / fc;

      t0 = t0 + w * h;
      s0 = s0 + w * hy;
      sx = sx + w * hxy;
      sy = sy + w * y * hy;
      sxy = sxy + w * y * hxy;
      syy = syy + w * y * y * hy;
      sxx = sxx + w * y * y * hx;
   }
   t0 = 0.5 * t0;
   s0 = 0.5 * s0;
   sx = 0.5 * sx;
   sy = 0.5 * sy;
   sxy = 0.5 * sxy;
   syy = 0.5 * syy;
   sxx = 0.5 * sxx;

   auto t = [&]( int a, int b ) {
      return v1c[ a ] * v1c[ b ] * s0 + ( v1c[ a ] * v21[ b ] + v21[ a ] * v1c[ b ] ) * sx +
             ( v31[ a ] * v1c[ b ] + v1c[ a ] * v31[ b ] ) * sy + ( v21[ a ] * v31[ b ] + v31[ a ] * v21[ b ] ) * sxy +
             v21[ a ] * v21[ b ] * sxx + v31[ a ] * v31[ b ] * syy;
   };
   const double t11 = t( 0, 0 ), t12 = t( 0, 1 ), t13 = t( 0, 2 ), t22 = t( 1, 1 ), t23 = t( 1, 2 ), t33 = t( 2, 2 );
   const double ar = mesh_.area[ n2 ];
   const double cg = cg_;
   const double g11 = 2.0 * ar * ( t0 + cg * t11 );
   const double g21 = 2.0 * ar * t12 * cg;
   const double g22 = 2.0 * ar * ( t0 + cg * t22 );
   const double g31 = 2.0 * ar * t13 * cg;
   const double g32 = 2.0 * ar * t23 * cg;
   const double g33 = 2.0 * ar * ( t0 + cg * t33 );
   g[ 0 ] = g11;
   g[ 1 ] = g21;
   g[ 2 ] = g31;
   g[ 3 ] = g21;
   g[ 4 ] = g22;
   g[ 5 ] = g32;
   g[ 6 ] = g31;
   g[ 7 ] = g32;
   g[ 8 ] = g33;
}

void Kernel::block( int n1, int n2, double* g ) const {
   if ( n1 == n2 ) {
      const double* s = &self_[ 9 * size_t( n1 ) ];
      for ( int i = 0; i < 9; ++i ) {
         g[ i ] = s[ i ];
      }
   } else {
      pair( n1, n2, g );
   }
}

void Kernel::fill( const int* rt, int m, const int* ct, int n, double* A, int lda ) const {
   double g[ 9 ];
   for ( int j = 0; j < n; ++j ) {
      for ( int i = 0; i < m; ++i ) {
         block( rt[ i ], ct[ j ], g );
         for ( int b = 0; b < 3; ++b ) {
            double* col = A + size_t( 3 * j + b ) * lda + 3 * i;
            col[ 0 ] = g[ b ];
            col[ 1 ] = g[ 3 + b ];
            col[ 2 ] = g[ 6 + b ];
         }
      }
   }
}


double Kernel::diag_mean() const {
   // best.f: DIAGM = DIAGM + GR(3N-2,3N-2) + GR(3N-1,3N-1) + GR(3N,3N), then DIAGM / NAR3.
   double s = 0.0;
   for ( int i = 0; i < n_; ++i ) {
      const double* g = &self_[ 9 * size_t( i ) ];
      s = s + g[ 0 ] + g[ 4 ] + g[ 8 ];
   }
   return s / ( 3.0 * n_ );
}

}  // namespace bh
