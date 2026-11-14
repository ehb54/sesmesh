// Unit tests of the Coalesce port (src/coalesce.h).  The byte-for-byte comparison with the
// Fortran rcoal is in test/test_rcoal.sh.
//
// 1. heap_index sorts, and gives a permutation when lengths are equal.
// 2. Coalesce on cleaned sesmesh meshes of random clusters, at several edge fractions: the
//    output is a closed, consistently oriented manifold of the input's topology, with no
//    triangle of zero area, and fewer triangles at larger edge fractions.
// 3. find_edgfrac reaches its target within 1%.
// 4. c3p_surface gives exactly the surface rcoal reads from the files write_c3p writes.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "clean.h"
#include "coalesce.h"
#include "io.h"
#include "mesh.h"
#include "rcoal_io.h"
#include "ses.h"

using namespace sesmesh;

namespace {

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

void test_heap_index() {
   std::mt19937 rng( 3 );
   for ( int n : { 1, 2, 3, 10, 1000 } ) {
      std::vector<double> a( n );
      for ( double& v : a ) {
         v = double( rng() % 50 );   // many ties
      }
      std::vector<int> ix;
      rcoal::heap_index( a, ix );
      std::vector<int> sorted = ix;
      std::sort( sorted.begin(), sorted.end() );
      bool perm = true;
      for ( int j = 0; j < n; ++j ) {
         perm = perm && sorted[j] == j;
      }
      bool asc = true;
      for ( int j = 1; j < n; ++j ) {
         asc = asc && a[ix[j - 1]] <= a[ix[j]];
      }
      CHECK( perm && asc, "heap_index, n = %d: %s", n, perm ? "not ascending" : "not a permutation" );
   }
}

// A sesmesh mesh as rcoal reads one: 1-based vertices, the edges in order of first
// appearance, all triangles in component 1.
rcoal::Surface surface( const Mesh& m ) {
   rcoal::Surface s;
   s.x.resize( m.x.size() + 1 );
   for ( size_t i = 0; i < m.x.size(); ++i ) {
      s.x[i + 1] = { m.x[i].x, m.x[i].y, m.x[i].z };
   }
   std::set<std::pair<int, int>> seen;
   for ( const auto& t : m.tri ) {
      rcoal::T3 u = { t[0] + 1, t[1] + 1, t[2] + 1 };
      s.tri.push_back( u );
      for ( int k = 0; k < 3; ++k ) {
         int a = u[k], b = u[( k + 1 ) % 3];
         if ( seen.insert( { std::min( a, b ), std::max( a, b ) } ).second ) {
            s.edges.push_back( { a, b } );
         }
      }
   }
   s.comp.assign( s.tri.size(), 1 );
   return s;
}

void test_cluster( unsigned seed, int natoms, double box, double p, double h ) {
   std::mt19937 rng( seed );
   auto uni = [&]( double a, double b ) { return a + ( b - a ) * ( rng() / 4294967296.0 ); };
   std::vector<Ball> atoms;
   for ( int n = 0; n < natoms; ++n ) {
      double x = uni( 0, box ), y = uni( 0, box ), z = uni( 0, box );
      atoms.push_back( { Vec3( x, y, z ), uni( 1.0, 2.0 ) } );
   }
   SesModel model( atoms, p );
   MeshOptions opt;
   opt.h = h;
   Mesh mesh;
   BuildInfo info;
   build_mesh( model, opt, mesh, info );
   CleanInfo ci;
   clean_mesh( model, h, CleanOptions(), mesh, ci );
   MeshCheck c0 = check_mesh( mesh );
   CHECK( c0.closed_manifold(), "cluster %u: input not a closed manifold", seed );
   rcoal::Surface s = surface( mesh );

   // c3p_surface is what rcoal reads from write_c3p's files.
   {
      const std::string base = "coalesce_test_c3p_" + std::to_string( seed );
      write_c3p( base, "test", mesh, natoms, p, c0.area, c0.volume );
      rcoal::Surface f = rcoal::read_c3p( base + ".c3p" );
      rcoal::Surface m = c3p_surface( mesh, c0.area, c0.volume );
      std::remove( ( base + ".c3p" ).c_str() );
      std::remove( ( base + ".c3v" ).c_str() );
      bool same = f.x == m.x && f.edges == m.edges && f.tri == m.tri && f.comp == m.comp &&
                  f.begin_area == m.begin_area && f.true_area == m.true_area && f.true_vol == m.true_vol;
      std::printf( "cluster %u: c3p_surface %s what rcoal reads from write_c3p's files\n", seed,
                   same ? "is" : "is NOT" );
      CHECK( same, "cluster %u: c3p_surface differs from read_c3p of write_c3p", seed );
   }

   int last = int( s.tri.size() ) + 1;
   for ( double e : { 0.5, 0.8, 1.2, 2.0 } ) {
      rcoal::Result R = rcoal::coalesce( s, e );
      Mesh out = mesh;
      out.tri.clear();
      for ( const auto& t : R.tri ) {
         out.tri.push_back( { t[0] - 1, t[1] - 1, t[2] - 1 } );
      }
      MeshCheck c = check_mesh( out );
      std::printf( "cluster %u, %ld triangles, e %.1f: %ld triangles, Euler %ld, %s, %ld of zero area\n", seed,
                   c0.triangles, e, c.triangles, c.euler, c.closed_manifold() ? "closed manifold" : "NOT closed",
                   c.degenerate );
      CHECK( c.closed_manifold(), "cluster %u, e %.1f: not a closed manifold", seed, e );
      CHECK( c.euler == c0.euler, "cluster %u, e %.1f: Euler %ld, input %ld", seed, e, c.euler, c0.euler );
      CHECK( c.degenerate == 0, "cluster %u, e %.1f: %ld triangles of zero area", seed, e, c.degenerate );
      CHECK( R.st.npoly == int( R.tri.size() ), "cluster %u, e %.1f: npoly %d, %zu triangles", seed, e, R.st.npoly,
             R.tri.size() );
      CHECK( int( c.triangles ) < last, "cluster %u, e %.1f: %ld triangles, not fewer than at the smaller e", seed, e,
             c.triangles );
      last = int( c.triangles );
   }

   int target = int( s.tri.size() ) / 4;
   int n = 0;
   double e = rcoal::find_edgfrac( s, target, &n );
   std::printf( "cluster %u: find_edgfrac(%d) = %.6f, %d triangles\n", seed, target, e, n );
   CHECK( std::abs( n - target ) <= target / 100, "cluster %u: find_edgfrac reached %d for target %d", seed, n, target );
   CHECK( rcoal::coalesce( s, e ).st.npoly == n, "cluster %u: find_edgfrac's count does not repeat", seed );
}

}  // namespace

int main() {
   test_heap_index();
   test_cluster( 7, 12, 6, 1.4, 0.15 );
   test_cluster( 11, 30, 9, 1.4, 0.2 );
   if ( failures ) {
      std::printf( "%d failure(s)\n", failures );
      return 1;
   }
   std::printf( "all coalesce tests passed\n" );
   return 0;
}
