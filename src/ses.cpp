// sesmesh: the exact solvent-excluded surface function; see ses.h.
#include "ses.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace sesmesh {

namespace {

const double PI = 3.14159265358979323846;
const double TWO_PI = 2 * PI;

}  // namespace

void PointHash::build( const std::vector<Vec3>& pts, double cell ) {
   cell_ = cell;
   idx_.clear();
   start_.clear();
   if ( pts.empty() ) {
      return;
   }
   lo_ = pts[0];
   Vec3 hi = pts[0];
   for ( const Vec3& p : pts ) {
      for ( int a = 0; a < 3; ++a ) {
         lo_[a] = std::min( lo_[a], p[a] );
         hi[a] = std::max( hi[a], p[a] );
      }
   }
   for ( int a = 0; a < 3; ++a ) {
      dim_[a] = int( ( hi[a] - lo_[a] ) / cell_ ) + 1;
   }
   long ncell = long( dim_[0] ) * dim_[1] * dim_[2];
   start_.assign( ncell + 1, 0 );
   std::vector<long> cellof( pts.size() );
   for ( size_t n = 0; n < pts.size(); ++n ) {
      int c[3];
      for ( int a = 0; a < 3; ++a ) {
         c[a] = std::min( dim_[a] - 1, int( ( pts[n][a] - lo_[a] ) / cell_ ) );
      }
      cellof[n] = ( long( c[2] ) * dim_[1] + c[1] ) * dim_[0] + c[0];
      ++start_[cellof[n] + 1];
   }
   for ( long c = 0; c < ncell; ++c ) {
      start_[c + 1] += start_[c];
   }
   idx_.resize( pts.size() );
   std::vector<long> fill( start_.begin(), start_.end() - 1 );
   for ( size_t n = 0; n < pts.size(); ++n ) {
      idx_[fill[cellof[n]]++] = int( n );
   }
}

SesModel::SesModel( const std::vector<Ball>& atoms, double probe ) : p_( probe ), natoms_( int( atoms.size() ) ) {
   if ( !( probe >= 0 ) ) {
      throw std::invalid_argument( "the probe radius must be zero or positive" );
   }
   if ( atoms.empty() ) {
      throw std::invalid_argument( "no atoms" );
   }
   for ( size_t a = 0; a < atoms.size(); ++a ) {
      if ( !( atoms[a].r >= 0 ) ) {
         throw std::invalid_argument( "atom " + std::to_string( a + 1 ) + " has a negative radius" );
      }
   }

   // Drop every sphere that lies inside another one: it changes neither U nor its boundary.
   // Of two identical spheres the first is kept.
   std::vector<Vec3> c( atoms.size() );
   double Rmax = 0;
   for ( size_t a = 0; a < atoms.size(); ++a ) {
      c[a] = atoms[a].c;
      Rmax = std::max( Rmax, atoms[a].r + p_ );
   }
   PointHash all;
   all.build( c, std::max( Rmax, 0.5 ) );
   const double eps = 1e-10;
   std::vector<int> keep;
   for ( size_t a = 0; a < atoms.size(); ++a ) {
      double Ra = atoms[a].r + p_;
      bool inside = false;
      all.visit( c[a], 2 * Rmax, [&]( int b ) {
         if ( inside || b == int( a ) ) {
            return;
         }
         double Rb = atoms[b].r + p_;
         double d = norm( c[a] - c[b] );
         if ( d + Ra < Rb - eps ) {
            inside = true;
         } else if ( std::fabs( d + Ra - Rb ) <= eps && ( Ra < Rb || ( Ra == Rb && b < int( a ) ) ) ) {
            inside = true;
         }
      } );
      if ( !inside ) {
         keep.push_back( int( a ) );
      }
   }

   s_.resize( keep.size() );
   std::vector<Vec3> kc( keep.size() );
   for ( size_t n = 0; n < keep.size(); ++n ) {
      Sph& s = s_[n];
      s.c = atoms[keep[n]].c;
      s.R = atoms[keep[n]].r + p_;
      s.R2 = s.R * s.R;
      s.atom = keep[n];
      s.buried = false;
      kc[n] = s.c;
      rmax_ = std::max( rmax_, s.R );
      for ( int a = 0; a < 3; ++a ) {
         if ( n == 0 || s.c[a] - s.R < lo_[a] ) {
            lo_[a] = s.c[a] - s.R;
         }
         if ( n == 0 || s.c[a] + s.R > hi_[a] ) {
            hi_[a] = s.c[a] + s.R;
         }
      }
   }
   hash_.build( kc, std::max( rmax_, 0.5 ) );

   // Neighbours: spheres that overlap.
   for ( size_t i = 0; i < s_.size(); ++i ) {
      s_[i].nb_begin = int( nbr_.size() );
      std::vector<int> nb;
      hash_.visit( s_[i].c, 2 * rmax_, [&]( int j ) {
         if ( j != int( i ) ) {
            double rr = s_[i].R + s_[j].R;
            if ( norm2( s_[i].c - s_[j].c ) < rr * rr ) {
               nb.push_back( j );
            }
         }
      } );
      std::sort( nb.begin(), nb.end() );
      nbr_.insert( nbr_.end(), nb.begin(), nb.end() );
      s_[i].nb_end = int( nbr_.size() );
   }

   // Circles where two spheres meet, and the arcs of each that lie inside no third sphere.
   std::vector<char> has_exposed( s_.size(), 0 );
   std::vector<std::pair<double, double>> buried, merged;
   for ( size_t i = 0; i < s_.size(); ++i ) {
      s_[i].cir_begin = int( circ_.size() );
      for ( int nn = s_[i].nb_begin; nn < s_[i].nb_end; ++nn ) {
         int j = nbr_[nn];
         if ( j <= int( i ) ) {
            continue;
         }
         const Sph& si = s_[i];
         const Sph& sj = s_[j];
         Vec3 dij = sj.c - si.c;
         double d = norm( dij );
         if ( d <= std::fabs( si.R - sj.R ) || d >= si.R + sj.R ) {
            continue;
         }
         Circle C;
         C.n = dij / d;
         double a = ( d * d + si.R2 - sj.R2 ) / ( 2 * d );
         double rho2 = si.R2 - a * a;
         if ( rho2 <= 0 ) {
            continue;
         }
         C.rho = std::sqrt( rho2 );
         C.m = si.c + C.n * a;
         basis( C.n, C.u, C.v );
         C.i = int( i );
         C.j = j;

         // Buried intervals of the angle th, where the circle point
         // m + rho (cos th u + sin th v) lies inside sphere k.
         bool all_buried = false;
         buried.clear();
         for ( int kk = si.nb_begin; kk < si.nb_end && !all_buried; ++kk ) {
            int k = nbr_[kk];
            if ( k == j ) {
               continue;
            }
            Vec3 w = s_[k].c - C.m;
            double du = dot( w, C.u ), dv = dot( w, C.v );
            double A = std::hypot( du, dv );
            double sv = ( rho2 + norm2( w ) - s_[k].R2 ) / ( 2 * C.rho );
            // inside k  <=>  du cos th + dv sin th > sv
            if ( A <= 1e-300 ) {
               if ( sv < 0 ) {
                  all_buried = true;
               }
               continue;
            }
            double cs = sv / A;
            if ( cs >= 1 ) {
               continue;
            }
            if ( cs <= -1 ) {
               all_buried = true;
               continue;
            }
            double half = std::acos( cs );
            double ph = std::atan2( dv, du );
            double lo = std::fmod( ph - half, TWO_PI );
            if ( lo < 0 ) {
               lo += TWO_PI;
            }
            double hi = lo + 2 * half;
            if ( hi > TWO_PI ) {
               buried.push_back( { lo, TWO_PI } );
               buried.push_back( { 0, hi - TWO_PI } );
            } else {
               buried.push_back( { lo, hi } );
            }
         }
         if ( all_buried ) {
            continue;
         }
         C.full = buried.empty();
         C.arc_begin = C.arc_end = int( arcs_.size() );
         C.tp_begin = C.tp_end = int( tp_.size() );
         if ( !C.full ) {
            std::sort( buried.begin(), buried.end() );
            merged.clear();
            for ( const auto& b : buried ) {
               if ( !merged.empty() && b.first <= merged.back().second ) {
                  merged.back().second = std::max( merged.back().second, b.second );
               } else {
                  merged.push_back( b );
               }
            }
            // Exposed arcs are the gaps, the last one wrapping through 0.
            const double tiny = 1e-13;
            for ( size_t g = 0; g < merged.size(); ++g ) {
               double lo = merged[g].second;
               double hi = g + 1 < merged.size() ? merged[g + 1].first : merged[0].first + TWO_PI;
               if ( hi - lo <= tiny ) {
                  continue;
               }
               if ( lo >= TWO_PI ) {
                  lo -= TWO_PI;
                  hi -= TWO_PI;
               }
               arcs_.push_back( { lo, hi } );
               tp_.push_back( C.m + ( C.u * std::cos( lo ) + C.v * std::sin( lo ) ) * C.rho );
               tp_.push_back( C.m + ( C.u * std::cos( hi ) + C.v * std::sin( hi ) ) * C.rho );
            }
            C.arc_end = int( arcs_.size() );
            C.tp_end = int( tp_.size() );
            if ( C.arc_end == C.arc_begin ) {
               continue;
            }
         }
         has_exposed[i] = has_exposed[j] = 1;
         circ_.push_back( C );
      }
      s_[i].cir_end = int( circ_.size() );
   }

   // A sphere with no exposed circle is either wholly exposed or wholly buried; one test point
   // tells which.
   const Vec3 dir = unit( Vec3( 0.26726, 0.53452, 0.80178 ) );
   for ( size_t i = 0; i < s_.size(); ++i ) {
      if ( !has_exposed[i] ) {
         s_[i].buried = !exposed_on_sphere( int( i ), s_[i].c + dir * s_[i].R );
      }
   }
}

int SesModel::spheres_buried() const {
   int n = 0;
   for ( const Sph& s : s_ ) {
      n += s.buried;
   }
   return n;
}

bool SesModel::exposed_on_sphere( int i, const Vec3& q ) const {
   for ( int n = s_[i].nb_begin; n < s_[i].nb_end; ++n ) {
      const Sph& o = s_[nbr_[n]];
      if ( norm2( q - o.c ) < o.R2 ) {
         return false;
      }
   }
   return true;
}

bool SesModel::in_arcs( const Circle& c, double th ) const {
   for ( int a = c.arc_begin; a < c.arc_end; ++a ) {
      const Arc& r = arcs_[a];
      if ( ( th >= r.lo && th <= r.hi ) || ( th + TWO_PI >= r.lo && th + TWO_PI <= r.hi ) ) {
         return true;
      }
   }
   return false;
}

void SesModel::gather( const Vec3& x, double radius, std::vector<int>& out ) const {
   out.clear();
   hash_.visit( x, radius + rmax_, [&]( int i ) {
      double r = s_[i].R + radius;
      if ( norm2( x - s_[i].c ) < r * r ) {
         out.push_back( i );
      }
   } );
}

double SesModel::eval( const Vec3& x, double dmax, EvalInfo* info ) const {
   std::vector<int> cand;
   gather( x, dmax, cand );
   return eval( x, cand, dmax, info );
}

double SesModel::eval( const Vec3& x, const std::vector<int>& cand, double dmax, EvalInfo* info ) const {
   bool inside = false;
   for ( int i : cand ) {
      if ( norm2( x - s_[i].c ) < s_[i].R2 ) {
         inside = true;
         break;
      }
   }

   if ( !inside ) {
      double best = dmax;
      int bi = -1;
      for ( int i : cand ) {
         double d = norm( x - s_[i].c ) - s_[i].R;
         if ( d < best ) {
            best = d;
            bi = i;
         }
      }
      if ( info ) {
         info->inside = false;
         info->feature = bi >= 0 ? F_SPHERE : F_NONE;
         info->atom = bi >= 0 ? s_[bi].atom : -1;
         info->q = bi >= 0 ? s_[bi].c + unit( x - s_[bi].c ) * s_[bi].R : x;
      }
      return p_ + best;
   }

   double best = dmax, best2 = dmax * dmax;
   int feat = F_NONE, fsph = -1;
   Vec3 q;

   // Exposed sphere patches: the nearest point of a sphere containing x is its radial projection.
   for ( int i : cand ) {
      const Sph& s = s_[i];
      if ( s.buried ) {
         continue;
      }
      Vec3 w = x - s.c;
      double d2 = norm2( w );
      if ( d2 >= s.R2 ) {
         continue;
      }
      double d = std::sqrt( d2 );
      double dist = s.R - d;
      if ( dist >= best || d <= 0 ) {
         continue;
      }
      Vec3 qq = s.c + w * ( s.R / d );
      if ( !exposed_on_sphere( i, qq ) ) {
         continue;
      }
      best = dist;
      best2 = dist * dist;
      feat = F_SPHERE;
      fsph = i;
      q = qq;
   }

   // Exposed circles and their exposed end points (the triple points).
   for ( int i : cand ) {
      const Sph& s = s_[i];
      for ( int ci = s.cir_begin; ci < s.cir_end; ++ci ) {
         const Circle& C = circ_[ci];
         Vec3 w = x - C.m;
         double lim = C.rho + best;
         if ( norm2( w ) >= lim * lim ) {
            continue;
         }
         double hh = dot( w, C.n );
         Vec3 wp = w - C.n * hh;
         double t = norm( wp );
         double dist2 = hh * hh + ( t - C.rho ) * ( t - C.rho );
         if ( dist2 < best2 ) {
            bool ok = true;
            Vec3 qq;
            if ( t > 1e-12 * C.rho ) {
               if ( !C.full ) {
                  double th = std::atan2( dot( wp, C.v ), dot( wp, C.u ) );
                  if ( th < 0 ) {
                     th += TWO_PI;
                  }
                  ok = in_arcs( C, th );
               }
               qq = C.m + wp * ( C.rho / t );
            } else {
               // x on the axis: every point of the circle is equally near.
               double th = C.full ? 0 : arcs_[C.arc_begin].lo;
               qq = C.m + ( C.u * std::cos( th ) + C.v * std::sin( th ) ) * C.rho;
            }
            if ( ok ) {
               best2 = dist2;
               best = std::sqrt( dist2 );
               feat = F_CIRCLE;
               fsph = C.i;
               q = qq;
            }
         }
         for ( int k = C.tp_begin; k < C.tp_end; ++k ) {
            double d2 = norm2( x - tp_[k] );
            if ( d2 < best2 ) {
               best2 = d2;
               best = std::sqrt( d2 );
               feat = F_TRIPLE;
               fsph = C.i;
               q = tp_[k];
            }
         }
      }
   }

   if ( info ) {
      info->inside = true;
      info->feature = feat;
      info->atom = fsph >= 0 ? s_[fsph].atom : -1;
      info->q = q;
   }
   return p_ - best;
}

}  // namespace sesmesh
