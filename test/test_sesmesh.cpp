// sesmesh unit tests.
//
// 1. The distance function: against the closed form for one sphere, and for random clusters
//    against a brute-force distance to a dense sampling of the exposed boundary of U.
// 2. Meshes of one, two and three spheres against closed-form area and volume, at several grid
//    spacings and extrapolated in h^2, and a two-sphere case whose surface pinches into two
//    pieces at a cusp.
// 3. A shell of spheres around a cavity: filled, and kept with -k.
// 4. A random cluster: the mesh must be a closed manifold.
// 5. The clean-up keeps every mesh a closed manifold of the same topology, hardly changes its
//    area and volume, and leaves no angle under 10 degrees.
// 6. The writers: what they write reads back.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "clean.h"
#include "io.h"
#include "mesh.h"
#include "ses.h"

using namespace sesmesh;

namespace {

const double PI = 3.14159265358979323846;
int failures = 0;

#define CHECK( cond, ... )                                      \
   do {                                                         \
      if ( !( cond ) ) {                                        \
         ++failures;                                            \
         std::printf( "FAIL %s:%d: ", __FILE__, __LINE__ );     \
         std::printf( __VA_ARGS__ );                            \
         std::printf( "\n" );                                   \
      }                                                         \
   } while ( 0 )

// Uniform numbers on [a, b) from mt19937, whose output the standard fixes;
// std::uniform_real_distribution differs between standard libraries.
struct Uniform {
   double a, b;
   Uniform( double lo, double hi ) : a( lo ), b( hi ) {}
   double operator()( std::mt19937& g ) const { return a + ( b - a ) * ( g() / 4294967296.0 ); }
   // x, y, z drawn in that order (the order of a call's arguments is unspecified)
   Vec3 vec( std::mt19937& g ) const {
      double x = ( *this )( g );
      double y = ( *this )( g );
      double z = ( *this )( g );
      return Vec3( x, y, z );
   }
};

struct Ref {
   double area, volume;
   int components;
};

// One sphere: the SES is the atom's own sphere.
Ref ref_one( double r ) { return { 4 * PI * r * r, 4 * PI * r * r * r / 3, 1 }; }

// Two spheres a distance d apart, by the meridian profile: atom arc, probe arc, atom arc,
// revolved about the axis (Pappus).  When the probe circle crosses the axis (rho < p) the
// surface pinches into two pieces at cusps.
Ref ref_two( double r1, double r2, double d, double p ) {
   double R1 = r1 + p, R2 = r2 + p;
   double a = ( d * d + R1 * R1 - R2 * R2 ) / ( 2 * d );
   double rho = std::sqrt( R1 * R1 - a * a );
   double x1 = a * r1 / R1, x2 = a + p * ( d - a ) / R2;   // contact points on the axis coordinate
   double u1 = x1 - a, u2 = x2 - a;
   double th1 = std::atan2( -rho, -a ), th2 = std::atan2( -rho, d - a );
   // pi * integral of (rho - sqrt(p^2 - u^2))^2 du
   auto I = [&]( double u ) {
      double s = std::sqrt( std::max( 0.0, p * p - u * u ) );
      double q = std::max( -1.0, std::min( 1.0, u / p ) );
      return ( rho * rho + p * p ) * u - u * u * u / 3 - rho * ( u * s + p * p * std::asin( q ) );
   };
   auto torus = [&]( double ta, double tb ) { return 2 * PI * p * ( rho * ( tb - ta ) - p * ( std::cos( tb ) - std::cos( ta ) ) ); };
   double y2 = x2 - d;
   double vcap = PI * ( r1 * r1 * x1 - x1 * x1 * x1 / 3 + 2 * r1 * r1 * r1 / 3 ) +
                 PI * ( 2 * r2 * r2 * r2 / 3 - r2 * r2 * y2 + y2 * y2 * y2 / 3 );
   double acap = 2 * PI * r1 * ( x1 + r1 ) + 2 * PI * r2 * ( d + r2 - x2 );
   if ( rho > p ) {
      return { acap + torus( th1, th2 ), vcap + PI * ( I( u2 ) - I( u1 ) ), 1 };
   }
   double w = std::sqrt( p * p - rho * rho );
   double tc1 = std::atan2( -rho, -w ), tc2 = std::atan2( -rho, w );
   return { acap + torus( th1, tc1 ) + torus( tc2, th2 ), vcap + PI * ( I( -w ) - I( u1 ) + I( u2 ) - I( w ) ), 2 };
}

// Three equal spheres (radius r) on an equilateral triangle of side s, the probe resting on them
// above and below (requires rho > p and H > p).  Area: contact patches (sphere minus two
// overlapping caps, the lens by Gauss-Bonnet), toroidal patches over the exposed arc, and
// reentrant geodesic triangles (spherical excess).  Volume: (1/3) of the integral of x.n over
// the same patches, with the vector area of each spherical patch from 1/2 of the loop integral
// of x cross dx along its boundary.
Ref ref_three( double r, double s, double p ) {
   double R = r + p, L = s / std::sqrt( 3.0 ), H = std::sqrt( R * R - L * L ), rho = std::sqrt( R * R - s * s / 4 );
   Vec3 c[3];
   for ( int k = 0; k < 3; ++k ) {
      c[k] = Vec3( L * std::cos( 2 * PI * k / 3 ), L * std::sin( 2 * PI * k / 3 ), 0 );
   }
   // reentrant patch over the probe at (0, 0, H)
   Vec3 P( 0, 0, H );
   Vec3 m[3];
   for ( int k = 0; k < 3; ++k ) {
      m[k] = unit( c[k] - P );
   }
   if ( dot( m[0], cross( m[1], m[2] ) ) < 0 ) {
      std::swap( m[1], m[2] );
   }
   double E = 2 * std::atan2( std::fabs( dot( m[0], cross( m[1], m[2] ) ) ),
                              1 + dot( m[0], m[1] ) + dot( m[1], m[2] ) + dot( m[2], m[0] ) );
   double are = p * p * E;
   Vec3 W;
   for ( int k = 0; k < 3; ++k ) {
      Vec3 cr = cross( m[k], m[( k + 1 ) % 3] );
      W += unit( cr ) * std::atan2( norm( cr ), dot( m[k], m[( k + 1 ) % 3] ) );
   }
   W = W * 0.5;
   double re_xn = H * p * p * norm( W ) - p * are;
   // toroidal patch of one pair
   double phi0 = std::acos( -( L / 2 ) / rho );
   double t1 = std::atan2( -rho, -s / 2 ), t2 = std::atan2( -rho, s / 2 );
   double S0 = t2 - t1, S1 = std::cos( t1 ) - std::cos( t2 );
   double S2 = ( t2 - t1 ) / 2 - ( std::sin( 2 * t2 ) - std::sin( 2 * t1 ) ) / 4;
   double ator = 2 * phi0 * p * ( rho * S0 + p * S1 );
   double tor_xn = -p * ( 2 * phi0 * ( p * rho * S0 + p * p * S1 ) + ( 2 * phi0 * rho + L * std::sin( phi0 ) ) * ( rho * S1 + p * S2 ) );
   // contact patch of atom 0
   double ca = ( s / 2 ) / R, sa = std::sqrt( 1 - ca * ca ), gam = PI / 3;
   double b0 = std::acos( ( std::cos( gam ) - ca * ca ) / ( sa * sa ) );
   double b1 = std::acos( ca * ( 1 - std::cos( gam ) ) / ( std::sin( gam ) * sa ) );
   double lens = 2 * ( PI - b0 - 2 * b1 * ca );
   double acon = r * r * ( 4 * PI - 4 * PI * ( 1 - ca ) + lens );
   Vec3 ej = unit( c[1] - c[0] ), ek = unit( c[2] - c[0] );
   Vec3 dj = unit( ek - ej * dot( ek, ej ) ), dk = unit( ej - ek * dot( ej, ek ) );
   Vec3 Wc = ( ej + ek ) * ( -PI * sa * sa ) + ( ( ej + ek ) * ( b1 * sa ) - ( dj + dk ) * ( ca * std::sin( b1 ) ) ) * sa;
   double con_xn = r * acon + r * r * dot( c[0], Wc );
   return { 3 * acon + 3 * ator + 2 * are, ( 3 * con_xn + 3 * tor_xn + 2 * re_xn ) / 3, 1 };
}

struct Fit {
   double area, volume;
};

// Cleans a copy of mesh and checks it against the raw mesh's check c.
void check_clean( const std::string& name, const SesModel& model, double h, const Mesh& mesh, const MeshCheck& c,
                  double ref_area, double ref_volume ) {
   Mesh m = mesh;
   CleanInfo ci;
   clean_mesh( model, h, CleanOptions(), m, ci );
   MeshCheck k = check_mesh( m );
   std::printf( "    cleaned: %ld triangles (%ld collapses, %ld flips, %ld moves), smallest angle %.1f deg (raw %.2f), "
                "%.2f%% under 20 deg (raw %.1f%%), area %+.2e, volume %+.2e\n",
                k.triangles, ci.collapses, ci.flips, ci.smoothed, k.min_angle, c.min_angle, 100 * k.below_20,
                100 * c.below_20, k.area / ref_area - 1, k.volume / ref_volume - 1 );
   CHECK( k.closed_manifold(), "%s h %g: cleaned mesh not a closed manifold", name.c_str(), h );
   CHECK( k.euler == c.euler && k.components.size() == c.components.size(),
          "%s h %g: clean-up changed the topology (Euler %ld -> %ld)", name.c_str(), h, c.euler, k.euler );
   CHECK( k.degenerate == 0, "%s h %g: cleaned mesh has %ld zero-area triangles", name.c_str(), h, k.degenerate );
   // No further from the reference than the raw mesh (for the clusters the raw mesh is the reference).
   CHECK( std::fabs( k.area - ref_area ) <= std::fabs( c.area - ref_area ) + 5e-4 * ref_area,
          "%s h %g: clean-up moved the area by %.2e", name.c_str(), h, k.area / c.area - 1 );
   CHECK( std::fabs( k.volume - ref_volume ) <= std::fabs( c.volume - ref_volume ) + 5e-4 * ref_volume,
          "%s h %g: clean-up moved the volume by %.2e", name.c_str(), h, k.volume / c.volume - 1 );
   CHECK( k.below_20 < 0.01 && k.below_20 < 0.1 * c.below_20, "%s h %g: %.2f%% of triangles still under 20 deg",
          name.c_str(), h, 100 * k.below_20 );
   CHECK( k.min_angle >= 10, "%s h %g: cleaned mesh has an angle of %.2f deg", name.c_str(), h, k.min_angle );
}

// Meshes at each spacing; checks topology; returns the h^2 extrapolation.
Fit mesh_series( const std::string& name, const std::vector<Ball>& atoms, double p, const Ref& ref,
                 const std::vector<double>& hs ) {
   SesModel model( atoms, p );
   std::vector<double> x, A, V;
   std::printf( "%s: reference area %.6f, volume %.6f\n", name.c_str(), ref.area, ref.volume );
   for ( double h : hs ) {
      MeshOptions opt;
      opt.h = h;
      Mesh mesh;
      BuildInfo info;
      build_mesh( model, opt, mesh, info );
      MeshCheck c = check_mesh( mesh );
      std::printf( "  h %.3f: %7ld triangles, %d component(s), Euler %3ld, area %.6f (%+.2e), volume %.6f (%+.2e)\n", h,
                   c.triangles, int( c.components.size() ), c.euler, c.area, c.area / ref.area - 1, c.volume,
                   c.volume / ref.volume - 1 );
      CHECK( c.closed_manifold(), "%s h %g: not a closed manifold (open %ld, non-manifold %ld/%ld, misoriented %ld)",
             name.c_str(), h, c.open_edges, c.nonmanifold_edges, c.nonmanifold_vertices, c.misoriented_edges );
      CHECK( int( c.components.size() ) == ref.components, "%s h %g: %d components, expected %d", name.c_str(), h,
             int( c.components.size() ), ref.components );
      CHECK( c.euler == 2 * ref.components, "%s h %g: Euler %ld", name.c_str(), h, c.euler );
      CHECK( c.degenerate == 0, "%s h %g: %ld zero-area triangles", name.c_str(), h, c.degenerate );
      for ( const ComponentInfo& ci : c.components ) {
         CHECK( ci.volume > 0, "%s h %g: a component has volume %g", name.c_str(), h, ci.volume );
      }
      if ( h >= 0.13 ) {
         check_clean( name, model, h, mesh, c, ref.area, ref.volume );
      }
      x.push_back( h * h );
      A.push_back( c.area );
      V.push_back( c.volume );
   }
   // least squares a + b h^2
   auto fit = [&]( const std::vector<double>& y ) {
      double n = double( x.size() ), sx = 0, sy = 0, sxx = 0, sxy = 0;
      for ( size_t i = 0; i < x.size(); ++i ) {
         sx += x[i];
         sy += y[i];
         sxx += x[i] * x[i];
         sxy += x[i] * y[i];
      }
      double b = ( n * sxy - sx * sy ) / ( n * sxx - sx * sx );
      return ( sy - b * sx ) / n;
   };
   Fit f = { fit( A ), fit( V ) };
   std::printf( "  extrapolated: area %.6f (%+.2e), volume %.6f (%+.2e)\n", f.area, f.area / ref.area - 1, f.volume,
                f.volume / ref.volume - 1 );
   return f;
}

// The area tolerance is wider where the surface has a conical cusp: near its tip the mesh error
// goes as h^2 log h, which the h^2 extrapolation does not remove.
void test_series( const std::string& name, const std::vector<Ball>& atoms, double p, const Ref& ref,
                  double area_tol = 2e-4 ) {
   const std::vector<double> hs = { 0.2, 0.16, 0.13, 0.1 };
   Fit f = mesh_series( name, atoms, p, ref, hs );
   CHECK( std::fabs( f.area / ref.area - 1 ) < area_tol, "%s: extrapolated area off by %.2e", name.c_str(),
          f.area / ref.area - 1 );
   CHECK( std::fabs( f.volume / ref.volume - 1 ) < 2e-4, "%s: extrapolated volume off by %.2e", name.c_str(),
          f.volume / ref.volume - 1 );
}

// --- distance function ---------------------------------------------------------------------

void test_one_sphere_function() {
   const Vec3 c( 0.3, -0.2, 0.1 );
   const double r = 1.7, p = 1.4, R = r + p;
   SesModel model( { { c, r } }, p );
   std::mt19937 rng( 1 );
   Uniform u( -5, 5 );
   double worst = 0;
   for ( int n = 0; n < 20000; ++n ) {
      Vec3 x = u.vec( rng );
      double d = norm( x - c );
      double expect = d < R ? d - r : p + ( d - R );
      worst = std::max( worst, std::fabs( model.eval( x, 100 ) - expect ) );
   }
   std::printf( "one sphere, f against |x - c| - r: worst error %.2e\n", worst );
   CHECK( worst < 1e-12, "one sphere: f off by %g", worst );
}

// Exposed boundary of U sampled densely, by itself, for a brute-force distance.
std::vector<Vec3> sample_boundary( const std::vector<Ball>& atoms, double p, int per_sphere, int per_circle, double& spacing ) {
   const size_t n = atoms.size();
   auto exposed = [&]( const Vec3& q, size_t skip1, size_t skip2, size_t skip3 ) {
      for ( size_t k = 0; k < n; ++k ) {
         if ( k == skip1 || k == skip2 || k == skip3 ) {
            continue;
         }
         double R = atoms[k].r + p;
         if ( norm2( q - atoms[k].c ) < R * R ) {
            return false;
         }
      }
      return true;
   };
   std::vector<Vec3> pts;
   spacing = 0;
   const double ga = PI * ( 3 - std::sqrt( 5.0 ) );
   for ( size_t i = 0; i < n; ++i ) {
      double R = atoms[i].r + p;
      spacing = std::max( spacing, std::sqrt( 4 * PI * R * R / per_sphere ) );
      for ( int k = 0; k < per_sphere; ++k ) {
         double z = 1 - ( 2.0 * k + 1 ) / per_sphere, rr = std::sqrt( 1 - z * z ), ph = ga * k;
         Vec3 q = atoms[i].c + Vec3( rr * std::cos( ph ), rr * std::sin( ph ), z ) * R;
         if ( exposed( q, i, n, n ) ) {
            pts.push_back( q );
         }
      }
   }
   for ( size_t i = 0; i < n; ++i ) {
      for ( size_t j = i + 1; j < n; ++j ) {
         double Ri = atoms[i].r + p, Rj = atoms[j].r + p;
         Vec3 dv = atoms[j].c - atoms[i].c;
         double d = norm( dv );
         if ( d >= Ri + Rj || d <= std::fabs( Ri - Rj ) ) {
            continue;
         }
         Vec3 e = dv / d, u, v;
         basis( e, u, v );
         double a = ( d * d + Ri * Ri - Rj * Rj ) / ( 2 * d ), rho = std::sqrt( Ri * Ri - a * a );
         for ( int k = 0; k < per_circle; ++k ) {
            double th = 2 * PI * k / per_circle;
            Vec3 q = atoms[i].c + e * a + ( u * std::cos( th ) + v * std::sin( th ) ) * rho;
            if ( exposed( q, i, j, n ) ) {
               pts.push_back( q );
            }
         }
         // triple points with every third sphere
         for ( size_t k = j + 1; k < n; ++k ) {
            double Rk = atoms[k].r + p;
            Vec3 w = atoms[k].c - ( atoms[i].c + e * a );
            double hh = dot( w, e ), wu = dot( w, u ), wv = dot( w, v );
            double A = std::hypot( wu, wv );
            if ( A < 1e-12 ) {
               continue;
            }
            double sv = ( rho * rho + norm2( w ) - Rk * Rk ) / ( 2 * rho );
            (void) hh;
            double cs = sv / A;
            if ( cs <= -1 || cs >= 1 ) {
               continue;
            }
            double ph = std::atan2( wv, wu ), half = std::acos( cs );
            for ( double th : { ph - half, ph + half } ) {
               Vec3 q = atoms[i].c + e * a + ( u * std::cos( th ) + v * std::sin( th ) ) * rho;
               if ( exposed( q, i, j, k ) ) {
                  pts.push_back( q );
               }
            }
         }
      }
   }
   return pts;
}

void test_cluster_function( unsigned seed, int natoms, double box, double p ) {
   std::mt19937 rng( seed );
   Uniform pos( 0, box ), rad( 1.0, 2.0 );
   std::vector<Ball> atoms;
   for ( int n = 0; n < natoms; ++n ) {
      atoms.push_back( { pos.vec( rng ), rad( rng ) } );
   }
   SesModel model( atoms, p );
   double spacing;
   std::vector<Vec3> pts = sample_boundary( atoms, p, 40000, 4000, spacing );
   Uniform q( -3, box + 3 );
   int inside = 0;
   double worst_low = 0, worst_high = 0;
   for ( int n = 0; n < 3000; ++n ) {
      Vec3 x = q.vec( rng );
      EvalInfo info;
      double f = model.eval( x, 100, &info );
      if ( !info.inside ) {
         double dout = 1e300;
         for ( const Ball& b : atoms ) {
            dout = std::min( dout, norm( x - b.c ) - b.r - p );
         }
         worst_high = std::max( worst_high, std::fabs( f - ( p + dout ) ) );
         continue;
      }
      ++inside;
      double d = p - f, dbf = 1e300;
      for ( const Vec3& s : pts ) {
         dbf = std::min( dbf, norm2( x - s ) );
      }
      dbf = std::sqrt( dbf );
      // the exact distance can never exceed the distance to a boundary point
      worst_low = std::max( worst_low, d - dbf );
      worst_high = std::max( worst_high, dbf - d );
      CHECK( std::fabs( norm( info.q - x ) - d ) < 1e-9, "cluster %u: nearest point not at the distance", seed );
   }
   std::printf( "cluster %u (%d spheres): %d interior points; exact - sampled distance: at most %.2e above, "
                "%.2e below (sampling %.3f A)\n",
                seed, natoms, inside, worst_low, worst_high, spacing );
   CHECK( worst_low < 1e-9, "cluster %u: distance %.3g above a sampled boundary point", seed, worst_low );
   CHECK( worst_high < spacing, "cluster %u: distance %.3g below the sampled boundary", seed, worst_high );
}

// --- meshes ---------------------------------------------------------------------------------

std::vector<Ball> shell( double radius, double r ) {
   const double g = ( 1 + std::sqrt( 5.0 ) ) / 2;
   std::vector<Vec3> v;
   for ( int s1 : { -1, 1 } ) {
      for ( int s2 : { -1, 1 } ) {
         v.push_back( Vec3( 0, s1, s2 * g ) );
         v.push_back( Vec3( s1, s2 * g, 0 ) );
         v.push_back( Vec3( s2 * g, 0, s1 ) );
         v.push_back( Vec3( 0, s1 * g, s2 / g ) );   // the dual dodecahedron
         v.push_back( Vec3( s1 / g, 0, s2 * g ) );
         v.push_back( Vec3( s1 * g, s2 / g, 0 ) );
         for ( int s3 : { -1, 1 } ) {
            v.push_back( Vec3( s1, s2, s3 ) );
         }
      }
   }
   std::vector<Ball> atoms;
   for ( const Vec3& x : v ) {
      atoms.push_back( { unit( x ) * radius, r } );
   }
   return atoms;
}

void test_cavity() {
   std::vector<Ball> atoms = shell( 5.0, 1.9 );
   SesModel model( atoms, 1.0 );
   MeshOptions opt;
   opt.h = 0.15;
   Mesh filled, kept;
   BuildInfo fi, ki;
   build_mesh( model, opt, filled, fi );
   opt.fill_cavities = false;
   build_mesh( model, opt, kept, ki );
   MeshCheck cf = check_mesh( filled ), ck = check_mesh( kept );
   check_clean( "cavity shell", model, 0.15, filled, cf, cf.area, cf.volume );
   check_clean( "cavity shell, kept", model, 0.15, kept, ck, ck.area, ck.volume );
   std::printf( "cavity shell (%zu spheres): filled: %d component(s), volume %.4f, %ld cavity nodes; kept: %d "
                "components, volumes",
                atoms.size(), int( cf.components.size() ), cf.volume, fi.cavity_nodes, int( ck.components.size() ) );
   for ( const ComponentInfo& c : ck.components ) {
      std::printf( " %.4f", c.volume );
   }
   std::printf( "\n" );
   CHECK( cf.closed_manifold() && ck.closed_manifold(), "cavity: not closed manifolds" );
   CHECK( cf.components.size() == 1 && cf.euler == 2, "cavity filled: %d components", int( cf.components.size() ) );
   CHECK( fi.cavity_nodes > 0, "cavity: no cavity found" );
   CHECK( ck.components.size() == 2, "cavity kept: %d components", int( ck.components.size() ) );
   if ( ck.components.size() == 2 ) {
      double outer = std::max( ck.components[0].volume, ck.components[1].volume );
      double inner = std::min( ck.components[0].volume, ck.components[1].volume );
      CHECK( inner < 0, "cavity kept: the cavity's surface faces the wrong way" );
      CHECK( std::fabs( outer - cf.volume ) < 1e-9 * outer, "cavity: outer surface differs when filled" );
   }
}

void test_cluster_mesh( unsigned seed, int natoms, double box, double p, double h ) {
   std::mt19937 rng( seed );
   Uniform pos( 0, box ), rad( 1.0, 2.0 );
   std::vector<Ball> atoms;
   for ( int n = 0; n < natoms; ++n ) {
      atoms.push_back( { pos.vec( rng ), rad( rng ) } );
   }
   SesModel model( atoms, p );
   MeshOptions opt;
   opt.h = h;
   Mesh mesh;
   BuildInfo info;
   build_mesh( model, opt, mesh, info );
   MeshCheck c = check_mesh( mesh );
   std::printf( "cluster %u mesh (%d spheres, h %.2f): %ld triangles, %d component(s), Euler %ld, %s\n", seed, natoms, h,
                c.triangles, int( c.components.size() ), c.euler, c.closed_manifold() ? "closed manifold" : "NOT closed" );
   CHECK( c.closed_manifold(), "cluster %u: not a closed manifold", seed );
   for ( const ComponentInfo& ci : c.components ) {
      CHECK( ci.volume > 0, "cluster %u: a component has volume %g after filling", seed, ci.volume );
   }
   check_clean( "cluster " + std::to_string( seed ), model, h, mesh, c, c.area, c.volume );
}

// Reads the whitespace-separated numbers of a line.
std::vector<double> numbers( const std::string& line ) {
   std::istringstream in( line );
   std::vector<double> v;
   std::string tok;
   while ( in >> tok ) {
      char* end = nullptr;
      double x = std::strtod( tok.c_str(), &end );
      if ( *end ) {
         break;
      }
      v.push_back( x );
   }
   return v;
}

void test_writers() {
   SesModel model( { { Vec3( 0.1, 0.2, 0.3 ), 1.7 }, { Vec3( 2.1, 0.2, 0.3 ), 1.5 } }, 1.4 );
   MeshOptions opt;
   opt.h = 0.25;
   Mesh m;
   BuildInfo info;
   build_mesh( model, opt, m, info );
   CleanInfo ci;
   clean_mesh( model, opt.h, CleanOptions(), m, ci );
   MeshCheck c = check_mesh( m );
   const std::string base = "sesmesh_test_io";
   write_msms( base, "test", m, 2, 1.4, opt.h );
   write_c3p( base, "test", m, 2, 1.4, 12.5, 34.5 );
   write_best( base, m, 12.5, 34.5 );
   write_bod( base, m, 14305, 20, 1.002 );
   const long nv = long( m.x.size() ), nt = long( m.tri.size() );
   std::string line;
   {
      std::ifstream f( base + ".vert" );
      std::getline( f, line );
      std::getline( f, line );
      std::getline( f, line );
      std::vector<double> h = numbers( line );
      long rows = 0;
      while ( std::getline( f, line ) ) {
         rows += numbers( line ).size() == 9;
      }
      CHECK( !h.empty() && long( h[0] ) == nv && rows == nv, "writers: .vert has %ld of %ld vertices", rows, nv );
   }
   {
      std::ifstream f( base + ".face" );
      for ( int k = 0; k < 3; ++k ) {
         std::getline( f, line );
      }
      long rows = 0;
      while ( std::getline( f, line ) ) {
         rows += numbers( line ).size() == 5;
      }
      CHECK( rows == nt, "writers: .face has %ld of %ld triangles", rows, nt );
   }
   {
      std::ifstream f( base + ".c3p" );
      std::getline( f, line );
      std::vector<double> h = numbers( line );
      CHECK( h.size() == 3 && long( h[0] ) == nv && long( h[1] ) == 3 * nt / 2 && long( h[2] ) == nt,
             "writers: .c3p header" );
      std::ifstream g( base + ".c3v" );
      for ( int k = 0; k < 13; ++k ) {
         std::getline( g, line );
      }
      std::vector<double> r = numbers( line );   // component, centre x y z, volume, area, accessible area
      CHECK( r.size() == 7 && r[4] == 34.5 && r[5] == 12.5, "writers: .c3v line 13 (%s)", line.c_str() );
   }
   {
      std::ifstream f( base + ".best" );
      std::getline( f, line );
      std::vector<double> h = numbers( line );
      CHECK( h.size() == 3 && long( h[0] ) == nt && h[1] == 12.5 && h[2] == 34.5, "writers: .best header" );
      double area = 0;
      long rows = 0;
      while ( std::getline( f, line ) ) {
         std::vector<double> v = numbers( line );
         if ( v.size() == 9 ) {
            ++rows;
            Vec3 a( v[0], v[1], v[2] ), b( v[3], v[4], v[5] ), cc( v[6], v[7], v[8] );
            area += 0.5 * norm( cross( b - a, cc - a ) );
         }
      }
      CHECK( rows == nt && std::fabs( area - c.area ) < 1e-6 * c.area, "writers: .best has %ld triangles, area %g", rows,
             area );
   }
   {
      std::ifstream f( base + ".bod" );
      long rows = 0, mass = 0;
      while ( std::getline( f, line ) ) {
         rows += line.compare( 0, 9, "TRIANGLE " ) == 0;
         mass += line == "MASS 14305 Da";
      }
      CHECK( rows == nt && mass == 1, "writers: .bod has %ld triangles", rows );
   }
   for ( const char* ext : { ".vert", ".face", ".c3p", ".c3v", ".best", ".bod" } ) {
      std::remove( ( base + ext ).c_str() );
   }
   std::printf( "writers: msms, c3p, best and bod read back (%ld vertices, %ld triangles)\n", nv, nt );
}

}  // namespace

int main() {
   test_one_sphere_function();
   test_cluster_function( 7, 12, 6, 1.4 );
   test_cluster_function( 11, 30, 9, 1.4 );

   test_series( "one sphere", { { Vec3( 0.11, 0.07, -0.05 ), 1.7 } }, 1.5, ref_one( 1.7 ) );
   {
      double r1 = 1.5, r2 = 1.8, d = 2.5, p = 1.4;
      test_series( "two overlapping spheres", { { Vec3( 0, 0, 0 ), r1 }, { unit( Vec3( 1, 0.3, 0.2 ) ) * d, r2 } }, p,
                   ref_two( r1, r2, d, p ) );
   }
   {
      double r = 1.5, d = 3.6, p = 1.4;
      test_series( "two separate spheres", { { Vec3( 0.05, 0, 0 ), r }, { Vec3( 0.05, d, 0 ), r } }, p,
                   ref_two( r, r, d, p ) );
   }
   {
      double r = 1.0, d = 4.2, p = 1.5;
      test_series( "two spheres, cusp", { { Vec3( 0, 0, 0 ), r }, { Vec3( 0, 0, d ), r } }, p, ref_two( r, r, d, p ),
                   1.5e-3 );
   }
   {
      double r = 1.6, s = 3.0, p = 1.4, L = s / std::sqrt( 3.0 );
      std::vector<Ball> atoms;
      for ( int k = 0; k < 3; ++k ) {
         atoms.push_back( { Vec3( L * std::cos( 2 * PI * k / 3 + 0.3 ), L * std::sin( 2 * PI * k / 3 + 0.3 ), 0.02 ), r } );
      }
      test_series( "three spheres", atoms, p, ref_three( r, s, p ) );
   }
   test_cavity();
   test_cluster_mesh( 7, 12, 6, 1.4, 0.15 );
   test_cluster_mesh( 11, 30, 9, 1.4, 0.12 );
   test_cluster_mesh( 13, 60, 12, 1.4, 0.2 );
   test_writers();

   if ( failures ) {
      std::printf( "%d check(s) FAILED\n", failures );
      return 1;
   }
   std::printf( "all sesmesh tests passed\n" );
   return 0;
}
