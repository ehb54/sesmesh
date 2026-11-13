// sesmesh: mesh clean-up; see clean.h.
#include "clean.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <utility>
#include <vector>

namespace sesmesh {

namespace {

const double PI = 3.14159265358979323846;

double now() {
   return std::chrono::duration<double>( std::chrono::steady_clock::now().time_since_epoch() ).count();
}

// Smallest angle of a triangle (radians); 0 if degenerate.
double min_angle( const Vec3& A, const Vec3& B, const Vec3& C ) {
   Vec3 ab = B - A, ac = C - A, bc = C - B;
   double cr = norm( cross( ab, ac ) );
   if ( !( cr > 0 ) ) {
      return 0;
   }
   double a1 = std::atan2( cr, dot( ab, ac ) );
   double a2 = std::atan2( cr, -dot( ab, bc ) );
   return std::min( a1, std::min( a2, PI - a1 - a2 ) );
}

// Moves x onto the surface: for x inside U at distance d from its nearest boundary point q,
// x + (q - x)(d - p)/d is at distance p from q, on the surface if q stays the nearest point.
bool project( const SesModel& model, double dmax, Vec3& x, Vec3& n, int& atom, int& feature ) {
   const double p = model.probe();
   if ( !( p > 0 ) ) {
      return false;
   }
   for ( int it = 0; it < 8; ++it ) {
      EvalInfo inf;
      double f = model.eval( x, dmax, &inf );
      if ( !inf.inside || inf.feature == F_NONE ) {
         return false;
      }
      Vec3 dir = inf.q - x;
      double d = norm( dir );
      if ( !( d > 0 ) ) {
         return false;
      }
      if ( std::fabs( f ) < 1e-11 ) {
         n = dir / d;
         atom = inf.atom;
         feature = inf.feature;
         return true;
      }
      x = x + dir * ( ( d - p ) / d );
   }
   return false;
}

struct Dyn {
   Mesh& m;
   std::vector<char> talive, valive;
   std::vector<std::vector<int>> vt;   // alive triangles at each vertex (pruned lazily)

   explicit Dyn( Mesh& mesh ) : m( mesh ) {
      talive.assign( m.tri.size(), 1 );
      valive.assign( m.x.size(), 1 );
      vt.resize( m.x.size() );
      for ( size_t t = 0; t < m.tri.size(); ++t ) {
         for ( int k = 0; k < 3; ++k ) {
            vt[m.tri[t][k]].push_back( int( t ) );
         }
      }
   }
   void prune( int v ) {
      std::vector<int>& l = vt[v];
      l.erase( std::remove_if( l.begin(), l.end(), [&]( int t ) { return !talive[t]; } ), l.end() );
   }
   static bool has( const std::array<int, 3>& t, int v ) { return t[0] == v || t[1] == v || t[2] == v; }
   void neighbours( int v, std::vector<int>& out ) {
      prune( v );
      out.clear();
      for ( int t : vt[v] ) {
         for ( int k = 0; k < 3; ++k ) {
            if ( m.tri[t][k] != v ) {
               out.push_back( m.tri[t][k] );
            }
         }
      }
      std::sort( out.begin(), out.end() );
      out.erase( std::unique( out.begin(), out.end() ), out.end() );
   }
   // The two triangles on edge (a, b); false unless there are exactly two.
   bool edge_tris( int a, int b, int& t1, int& t2 ) {
      prune( a );
      t1 = t2 = -1;
      int n = 0;
      for ( int t : vt[a] ) {
         if ( has( m.tri[t], b ) ) {
            ( n == 0 ? t1 : t2 ) = t;
            ++n;
         }
      }
      return n == 2;
   }
   static int third( const std::array<int, 3>& t, int a, int b ) {
      for ( int k = 0; k < 3; ++k ) {
         if ( t[k] != a && t[k] != b ) {
            return t[k];
         }
      }
      return -1;
   }
};

}  // namespace

void clean_mesh( const SesModel& model, double h, const CleanOptions& opt, Mesh& m, CleanInfo& info ) {
   double t0 = now();
   info = CleanInfo();
   Dyn d( m );
   const double dmax = model.probe() + 2 * h;
   const double lmax = opt.collapse * h;
   const double tiny = 1e-3 * h;
   const double cos_flip = std::cos( opt.flip_dihedral * PI / 180 );
   const double good = opt.good_angle * PI / 180;
   std::vector<int> na, nb, nc, around;

   // Smallest angle of triangle t with vertex `from` moved to P (from < 0: as it is).
   auto angle_with = [&]( int t, int from, const Vec3& P ) {
      const std::array<int, 3>& tr = m.tri[t];
      Vec3 X[3];
      for ( int k = 0; k < 3; ++k ) {
         X[k] = tr[k] == from ? P : m.x[tr[k]];
      }
      return min_angle( X[0], X[1], X[2] );
   };
   // True if moving `from` to P (with surface normal nP) keeps triangle t facing the same way.
   auto keeps_facing = [&]( int t, int from, const Vec3& P, const Vec3& nP ) {
      const std::array<int, 3>& tr = m.tri[t];
      Vec3 X[3], Y[3];
      for ( int k = 0; k < 3; ++k ) {
         X[k] = m.x[tr[k]];
         Y[k] = tr[k] == from ? P : X[k];
      }
      Vec3 n0 = cross( X[1] - X[0], X[2] - X[0] ), n1 = cross( Y[1] - Y[0], Y[2] - Y[0] );
      double l0 = norm( n0 ), l1 = norm( n1 );
      if ( !( l1 > 0 ) || dot( n1, nP ) <= 0 ) {
         return false;
      }
      return l0 < 1e-12 * l1 || dot( n0, n1 ) > 0.5 * l0 * l1;
   };

   for ( int pass = 0; pass < opt.passes; ++pass ) {
      long collapsed = 0, flipped = 0;

      // Edge collapses, shortest first.
      std::vector<std::pair<double, std::pair<int, int>>> edges;
      for ( size_t t = 0; t < m.tri.size(); ++t ) {
         if ( !d.talive[t] ) {
            continue;
         }
         for ( int k = 0; k < 3; ++k ) {
            int a = m.tri[t][k], b = m.tri[t][( k + 1 ) % 3];
            if ( a < b ) {
               double l = norm( m.x[a] - m.x[b] );
               if ( l < lmax ) {
                  edges.push_back( { l, { a, b } } );
               }
            }
         }
      }
      std::sort( edges.begin(), edges.end() );
      for ( const auto& e : edges ) {
         int a = e.second.first, b = e.second.second;
         if ( !d.valive[a] || !d.valive[b] ) {
            continue;
         }
         double len = norm( m.x[a] - m.x[b] );
         int t1, t2;
         if ( len >= lmax || !d.edge_tris( a, b, t1, t2 ) ) {
            continue;
         }
         int c = Dyn::third( m.tri[t1], a, b ), dd = Dyn::third( m.tri[t2], a, b );
         d.neighbours( a, na );
         d.neighbours( b, nb );
         std::vector<int> common;
         std::set_intersection( na.begin(), na.end(), nb.begin(), nb.end(), std::back_inserter( common ) );
         d.neighbours( c, nc );
         size_t vc = nc.size();
         d.neighbours( dd, nc );
         size_t vd = nc.size();
         if ( common.size() != 2 || c == dd || vc <= 3 || vd <= 3 || na.size() + nb.size() > 24 ) {
            ++info.rejected_link;
            continue;
         }
         // Triangles that survive, and the smallest angle around the edge now.
         around.clear();
         double worst_old = PI;
         for ( int v : { a, b } ) {
            for ( int t : d.vt[v] ) {
               worst_old = std::min( worst_old, angle_with( t, -1, Vec3() ) );
               if ( t != t1 && t != t2 ) {
                  around.push_back( t );
               }
            }
         }
         // Candidate positions: either end, or the projected midpoint.
         struct Cand {
            Vec3 x, n;
            int atom, feature;
         };
         Cand cand[3] = { { m.x[a], m.n[a], m.atom[a], m.feature[a] }, { m.x[b], m.n[b], m.atom[b], m.feature[b] },
                          { ( m.x[a] + m.x[b] ) * 0.5, Vec3(), -1, F_NONE } };
         int ncand = 2;
         if ( project( model, dmax, cand[2].x, cand[2].n, cand[2].atom, cand[2].feature ) &&
              norm( cand[2].x - ( m.x[a] + m.x[b] ) * 0.5 ) < 0.25 * h ) {
            ncand = 3;
         }
         int best = -1;
         double best_worst = -1;
         for ( int ci = 0; ci < ncand; ++ci ) {
            bool ok = true;
            double worst = PI;
            for ( int t : around ) {
               int from = Dyn::has( m.tri[t], a ) ? a : b;
               if ( !keeps_facing( t, from, cand[ci].x, cand[ci].n ) ) {
                  ok = false;
                  break;
               }
               worst = std::min( worst, angle_with( t, from, cand[ci].x ) );
            }
            if ( ok && worst > best_worst ) {
               best_worst = worst;
               best = ci;
            }
         }
         if ( best < 0 || !( best_worst > worst_old || best_worst > good || len < tiny ) ) {
            ++info.rejected_shape;
            continue;
         }
         // Collapse b into a.
         for ( int t : d.vt[b] ) {
            if ( !d.talive[t] ) {
               continue;
            }
            if ( Dyn::has( m.tri[t], a ) ) {
               d.talive[t] = 0;
            } else {
               for ( int k = 0; k < 3; ++k ) {
                  if ( m.tri[t][k] == b ) {
                     m.tri[t][k] = a;
                  }
               }
               d.vt[a].push_back( t );
            }
         }
         d.vt[b].clear();
         d.valive[b] = 0;
         m.x[a] = cand[best].x;
         m.n[a] = cand[best].n;
         m.atom[a] = cand[best].atom;
         m.feature[a] = cand[best].feature;
         d.prune( a );
         d.prune( c );
         d.prune( dd );
         ++collapsed;
      }

      // Edge flips where they raise the smallest angle of two nearly coplanar triangles.
      for ( size_t t = 0; t < m.tri.size(); ++t ) {
         if ( !d.talive[t] ) {
            continue;
         }
         for ( int k = 0; k < 3; ++k ) {
            int a = m.tri[t][k], b = m.tri[t][( k + 1 ) % 3];
            int t1, t2;
            if ( !d.edge_tris( a, b, t1, t2 ) ) {
               continue;
            }
            if ( t2 == int( t ) ) {
               std::swap( t1, t2 );
            }
            // t1 = (a, b, c) and t2 = (b, a, dd), counter-clockwise; the quadrilateral is b, c, a, dd.
            int c = Dyn::third( m.tri[t1], a, b ), dd = Dyn::third( m.tri[t2], a, b );
            const Vec3 &A = m.x[a], &B = m.x[b], &C = m.x[c], &D = m.x[dd];
            Vec3 n1 = cross( B - A, C - A ), n2 = cross( A - B, D - B );
            double l1 = norm( n1 ), l2 = norm( n2 );
            if ( !( l1 > 0 && l2 > 0 ) || dot( n1, n2 ) < cos_flip * l1 * l2 ) {
               continue;
            }
            double old_min = std::min( min_angle( A, B, C ), min_angle( B, A, D ) );
            double new_min = std::min( min_angle( C, A, D ), min_angle( D, B, C ) );
            if ( !( new_min > old_min + 1e-3 ) ) {
               continue;
            }
            Vec3 m1 = cross( A - C, D - C ), m2 = cross( B - D, C - D );
            Vec3 avg = n1 / l1 + n2 / l2;
            if ( dot( m1, avg ) <= 0 || dot( m2, avg ) <= 0 ) {
               continue;
            }
            d.neighbours( a, na );
            d.neighbours( b, nb );
            d.neighbours( c, nc );
            if ( na.size() <= 3 || nb.size() <= 3 || std::binary_search( nc.begin(), nc.end(), dd ) ) {
               continue;
            }
            m.tri[t1] = { c, a, dd };
            m.tri[t2] = { dd, b, c };
            auto drop = [&]( int v, int tt ) {
               std::vector<int>& l = d.vt[v];
               l.erase( std::remove( l.begin(), l.end(), tt ), l.end() );
            };
            drop( a, t2 );
            drop( b, t1 );
            d.vt[c].push_back( t2 );
            d.vt[dd].push_back( t1 );
            ++flipped;
            break;   // this triangle has changed
         }
      }
      // Tangential smoothing: each vertex moves halfway to its neighbours' centroid in the
      // tangent plane, then back onto the surface, if that raises the smallest angle around it.
      long moved = 0;
      if ( pass < opt.smooth ) {
         for ( size_t v = 0; v < m.x.size(); ++v ) {
            if ( !d.valive[v] ) {
               continue;
            }
            d.neighbours( int( v ), na );
            if ( na.size() < 3 ) {
               continue;
            }
            Vec3 c;
            for ( int u : na ) {
               c += m.x[u];
            }
            c = c / double( na.size() );
            Vec3 del = c - m.x[v];
            del = del - m.n[v] * dot( del, m.n[v] );
            Vec3 P = m.x[v] + del * 0.5, nP;
            int at = -1, fe = F_NONE;
            if ( !project( model, dmax, P, nP, at, fe ) || norm( P - m.x[v] ) > 0.5 * h ) {
               continue;
            }
            double worst_old = PI, worst_new = PI;
            bool ok = true;
            for ( int t : d.vt[v] ) {
               worst_old = std::min( worst_old, angle_with( t, -1, Vec3() ) );
               if ( !keeps_facing( t, int( v ), P, nP ) ) {
                  ok = false;
                  break;
               }
               worst_new = std::min( worst_new, angle_with( t, int( v ), P ) );
            }
            if ( !ok || !( worst_new > worst_old ) ) {
               continue;
            }
            m.x[v] = P;
            m.n[v] = nP;
            m.atom[v] = at;
            m.feature[v] = fe;
            ++moved;
         }
      }
      info.collapses += collapsed;
      info.flips += flipped;
      info.smoothed += moved;
      if ( collapsed == 0 && flipped == 0 && moved == 0 ) {
         break;
      }
   }

   // Compact.
   std::vector<int> newid( m.x.size(), -1 );
   Mesh out;
   for ( size_t v = 0; v < m.x.size(); ++v ) {
      if ( d.valive[v] ) {
         newid[v] = int( out.x.size() );
         out.x.push_back( m.x[v] );
         out.n.push_back( m.n[v] );
         out.atom.push_back( m.atom[v] );
         out.feature.push_back( m.feature[v] );
      }
   }
   for ( size_t t = 0; t < m.tri.size(); ++t ) {
      if ( d.talive[t] ) {
         out.tri.push_back( { newid[m.tri[t][0]], newid[m.tri[t][1]], newid[m.tri[t][2]] } );
      }
   }
   m = std::move( out );
   info.seconds = now() - t0;
}

}  // namespace sesmesh
