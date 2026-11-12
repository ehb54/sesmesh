// sesmesh: triangulate the solvent-excluded surface of a set of spheres.
//
//   sesmesh [options] input.xyzr
//
// See README.md at the top of the repository.
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mesh.h"
#include "ses.h"

using namespace sesmesh;

namespace {

const char* USAGE =
   "usage: sesmesh [options] input.xyzr\n"
   "\n"
   "Triangulates the solvent-excluded surface (SES) of the spheres in input.xyzr, one\n"
   "sphere per line as x y z r (further columns are ignored; # starts a comment).\n"
   "\n"
   "  -p probe     probe radius in A (default 1.5)\n"
   "  -g h         grid spacing in A (default 0.2); vertices lie on the exact surface\n"
   "  -x h2,h3,..  also mesh at these spacings, report area and volume for each, and\n"
   "               extrapolate both to h = 0 by a least-squares fit a + b h^2\n"
   "  -o base      output base name (default: the input's name without .xyzr)\n"
   "  -f list      output formats, a comma list of msms (base.vert, base.face) and\n"
   "               c3p (base.c3p, base.c3v as msroll writes them, for rcoal), or none\n"
   "               (default msms,c3p)\n"
   "  -k           keep internal cavities (default: fill them, keeping the outer surface)\n"
   "  -t n         threads (default: all)\n"
   "  -v           report the time each stage takes\n"
   "\n"
   "base.c3v gets the extrapolated area and volume when -x is given, else the mesh's own.\n";

[[noreturn]] void fail( const std::string& msg ) {
   std::fprintf( stderr, "sesmesh: %s\n", msg.c_str() );
   std::exit( 1 );
}

double to_double( const std::string& s, const char* what ) {
   char* end = nullptr;
   errno = 0;
   double v = std::strtod( s.c_str(), &end );
   if ( s.empty() || *end || errno ) {
      fail( std::string( "bad " ) + what + ": " + s );
   }
   return v;
}

std::vector<std::string> split( const std::string& s, char sep ) {
   std::vector<std::string> out;
   std::string cur;
   std::istringstream in( s );
   while ( std::getline( in, cur, sep ) ) {
      out.push_back( cur );
   }
   return out;
}

std::vector<Ball> read_xyzr( const std::string& path ) {
   std::ifstream in( path );
   if ( !in ) {
      fail( "cannot read " + path );
   }
   std::vector<Ball> atoms;
   std::string line;
   long lineno = 0;
   while ( std::getline( in, line ) ) {
      ++lineno;
      size_t hash = line.find( '#' );
      if ( hash != std::string::npos ) {
         line.erase( hash );
      }
      std::istringstream ls( line );
      std::string tok[4];
      int n = 0;
      while ( n < 4 && ls >> tok[n] ) {
         ++n;
      }
      if ( n == 0 ) {
         continue;
      }
      if ( n < 4 ) {
         fail( path + " line " + std::to_string( lineno ) + ": expected x y z r" );
      }
      double v[4];
      for ( int k = 0; k < 4; ++k ) {
         char* end = nullptr;
         v[k] = std::strtod( tok[k].c_str(), &end );
         if ( *end || !std::isfinite( v[k] ) ) {
            fail( path + " line " + std::to_string( lineno ) + ": not a number: " + tok[k] );
         }
      }
      atoms.push_back( { Vec3( v[0], v[1], v[2] ), v[3] } );
   }
   if ( atoms.empty() ) {
      fail( "no spheres in " + path );
   }
   return atoms;
}

// MSMS face types: 1 contact, 2 reentrant, 3 toroidal.
int msms_type( int feature ) { return feature == F_SPHERE ? 1 : ( feature == F_TRIPLE ? 2 : 3 ); }

FILE* open_out( const std::string& path ) {
   FILE* f = std::fopen( path.c_str(), "w" );
   if ( !f ) {
      fail( "cannot write " + path );
   }
   return f;
}

void close_out( FILE* f, const std::string& path ) {
   if ( std::ferror( f ) || std::fclose( f ) != 0 ) {
      fail( "error writing " + path );
   }
}

void write_msms( const std::string& base, const std::string& name, const Mesh& m, int natoms, double probe, double h ) {
   std::string pv = base + ".vert", pf = base + ".face";
   FILE* f = open_out( pv );
   std::fprintf( f, "# MSMS solvent excluded surface vertices for %s (sesmesh)\n", name.c_str() );
   std::fprintf( f, "#vertex #sphere density probe_r\n" );
   std::fprintf( f, "%7ld %7d %6.2f %5.2f\n", long( m.x.size() ), natoms, 1 / ( h * h ), probe );
   for ( size_t v = 0; v < m.x.size(); ++v ) {
      std::fprintf( f, "%14.8f %14.8f %14.8f %10.6f %10.6f %10.6f %7d %7d %2d\n", m.x[v].x, m.x[v].y, m.x[v].z,
                    m.n[v].x, m.n[v].y, m.n[v].z, 0, m.atom[v] + 1, msms_type( m.feature[v] ) );
   }
   close_out( f, pv );
   f = open_out( pf );
   std::fprintf( f, "# MSMS solvent excluded surface faces for %s (sesmesh)\n", name.c_str() );
   std::fprintf( f, "#faces  #sphere density probe_r\n" );
   std::fprintf( f, "%7ld %7d %6.2f %5.2f\n", long( m.tri.size() ), natoms, 1 / ( h * h ), probe );
   for ( const auto& t : m.tri ) {
      std::fprintf( f, "%7d %7d %7d %2d %6d\n", t[0] + 1, t[1] + 1, t[2] + 1, msms_type( m.feature[t[0]] ), 0 );
   }
   close_out( f, pf );
}

// msroll's .c3p (vertices with normals, edges, triangles with their edges) and .c3v (areas and
// volume; rcoal reads the area and volume of component 1 from its line 13).
void write_c3p( const std::string& base, const std::string& name, const Mesh& m, int natoms, double probe,
                double area, double volume ) {
   const long nt = long( m.tri.size() );
   std::vector<std::pair<uint64_t, long>> key( 3 * nt );
   for ( long t = 0; t < nt; ++t ) {
      for ( int k = 0; k < 3; ++k ) {
         uint64_t a = uint64_t( m.tri[t][k] ), b = uint64_t( m.tri[t][( k + 1 ) % 3] );
         key[3 * t + k] = { ( std::min( a, b ) << 32 ) | std::max( a, b ), 3 * t + k };
      }
   }
   std::sort( key.begin(), key.end() );
   std::vector<long> eid( 3 * nt );
   std::vector<uint64_t> edges;
   for ( size_t i = 0; i < key.size(); ++i ) {
      if ( i == 0 || key[i].first != key[i - 1].first ) {
         edges.push_back( key[i].first );
      }
      eid[key[i].second] = long( edges.size() );
   }

   std::string pc = base + ".c3p", pv = base + ".c3v";
   FILE* f = open_out( pc );
   std::fprintf( f, "%ld %ld %ld\n", long( m.x.size() ), long( edges.size() ), nt );
   for ( size_t v = 0; v < m.x.size(); ++v ) {
      std::fprintf( f, "%15.8f %15.8f %15.8f %9.6f %9.6f %9.6f   0.000000   0.000000   0.000000   1 %6d %2d\n",
                    m.x[v].x, m.x[v].y, m.x[v].z, m.n[v].x, m.n[v].y, m.n[v].z, m.atom[v] + 1,
                    msms_type( m.feature[v] ) );
   }
   for ( uint64_t e : edges ) {
      std::fprintf( f, "%7ld %7ld   1     0   0\n", long( e >> 32 ) + 1, long( e & 0xffffffffu ) + 1 );
   }
   for ( long t = 0; t < nt; ++t ) {
      std::fprintf( f, "%8ld %8ld %8ld %7d %7d %7d   1     1   0\n", eid[3 * t], eid[3 * t + 1], eid[3 * t + 2],
                    m.tri[t][0] + 1, m.tri[t][1] + 1, m.tri[t][2] + 1 );
   }
   close_out( f, pc );

   Vec3 c;
   for ( const Vec3& x : m.x ) {
      c += x;
   }
   if ( !m.x.empty() ) {
      c = c / double( m.x.size() );
   }
   f = open_out( pv );
   std::fprintf( f, "molecule name           = %s\n", name.c_str() );
   std::fprintf( f, "number of atoms         = %10d\n", natoms );
   std::fprintf( f, "probe radius            = %14.3f\n", probe );
   std::fprintf( f, "contact area            = %14.3f\n", 0.0 );
   std::fprintf( f, "reentrant area          = %14.3f\n", 0.0 );
   std::fprintf( f, "molecular area          = %14.3f\n", area );
   std::fprintf( f, "accessible area         = %14.3f\n", 0.0 );
   std::fprintf( f, "solvent-excluded volume = %14.3f\n", volume );
   std::fprintf( f, "\n" );
   std::fprintf( f, "component           center              volume    molecular  accessible\n" );
   std::fprintf( f, " number      x        y        z                    area       area\n" );
   std::fprintf( f, "\n" );
   std::fprintf( f, "    1  %10.3f %8.3f %8.3f %11.3f %11.3f %11.3f\n", c.x, c.y, c.z, volume, area, 0.0 );
   close_out( f, pv );
}

// Least-squares fit y = a + b x; returns a (and b).
double fit_line( const std::vector<double>& x, const std::vector<double>& y, double& b ) {
   const size_t n = x.size();
   double sx = 0, sy = 0, sxx = 0, sxy = 0;
   for ( size_t i = 0; i < n; ++i ) {
      sx += x[i];
      sy += y[i];
      sxx += x[i] * x[i];
      sxy += x[i] * y[i];
   }
   double d = n * sxx - sx * sx;
   b = d != 0 ? ( n * sxy - sx * sy ) / d : 0;
   return ( sy - b * sx ) / n;
}

}  // namespace

int main( int argc, char** argv ) {
   double probe = 1.5, h = 0.2;
   std::vector<double> extra;
   std::string input, base;
   bool want_msms = true, want_c3p = true, fill = true, verbose = false;
   int threads = 0;

   for ( int a = 1; a < argc; ++a ) {
      std::string o = argv[a];
      auto value = [&]() -> std::string {
         if ( a + 1 >= argc ) {
            fail( "option " + o + " needs a value" );
         }
         return argv[++a];
      };
      if ( o == "-h" || o == "--help" ) {
         std::fputs( USAGE, stdout );
         return 0;
      } else if ( o == "-p" ) {
         probe = to_double( value(), "probe radius" );
      } else if ( o == "-g" ) {
         h = to_double( value(), "grid spacing" );
      } else if ( o == "-x" ) {
         for ( const std::string& s : split( value(), ',' ) ) {
            extra.push_back( to_double( s, "grid spacing" ) );
         }
      } else if ( o == "-o" ) {
         base = value();
      } else if ( o == "-f" ) {
         want_msms = want_c3p = false;
         for ( const std::string& s : split( value(), ',' ) ) {
            if ( s == "msms" ) {
               want_msms = true;
            } else if ( s == "c3p" ) {
               want_c3p = true;
            } else if ( s != "none" ) {
               fail( "unknown format " + s );
            }
         }
      } else if ( o == "-k" ) {
         fill = false;
      } else if ( o == "-t" ) {
         threads = int( to_double( value(), "thread count" ) );
      } else if ( o == "-v" ) {
         verbose = true;
      } else if ( !o.empty() && o[0] == '-' ) {
         fail( "unknown option " + o + "\n" + USAGE );
      } else if ( input.empty() ) {
         input = o;
      } else {
         fail( "more than one input file\n" + std::string( USAGE ) );
      }
   }
   if ( input.empty() ) {
      std::fputs( USAGE, stderr );
      return 1;
   }
   if ( !( probe >= 0 ) ) {
      fail( "the probe radius must be zero or positive" );
   }
   if ( !( h > 0 ) ) {
      fail( "the grid spacing must be positive" );
   }
   for ( double x : extra ) {
      if ( !( x > 0 ) ) {
         fail( "grid spacings must be positive" );
      }
   }
   std::string name = input.substr( input.find_last_of( '/' ) == std::string::npos ? 0 : input.find_last_of( '/' ) + 1 );
   if ( name.size() > 5 && name.compare( name.size() - 5, 5, ".xyzr" ) == 0 ) {
      name.erase( name.size() - 5 );
   }
   if ( base.empty() ) {
      base = name;
   }
#ifdef _OPENMP
   if ( threads > 0 ) {
      omp_set_num_threads( threads );
   }
   threads = omp_get_max_threads();
#else
   threads = 1;
#endif

   std::vector<Ball> atoms = read_xyzr( input );
   try {
      SesModel model( atoms, probe );
      std::printf( "sesmesh: %s, %d spheres (%d inside others, %d buried), probe %.3f A; %d exposed circles, "
                   "%d exposed arcs; %d threads\n",
                   input.c_str(), model.atoms(), model.atoms() - model.spheres_kept(), model.spheres_buried(), probe,
                   model.circles_exposed(), model.arcs_exposed(), threads );

      std::vector<double> hs( extra );
      hs.push_back( h );   // the mesh written out is the last one
      std::vector<double> h2, areas, vols;
      Mesh mesh;
      MeshCheck chk;
      for ( size_t k = 0; k < hs.size(); ++k ) {
         MeshOptions opt;
         opt.h = hs[k];
         opt.fill_cavities = fill;
         BuildInfo info;
         build_mesh( model, opt, mesh, info );
         chk = check_mesh( mesh );
         std::printf( "h %.4f: grid %dx%dx%d, %ld band bricks; %ld vertices, %ld triangles, %zu component%s, %s, "
                      "Euler %ld; area %.4f, volume %.4f; smallest angle %.2f deg, %.1f%% under 10 deg\n",
                      hs[k], info.nodes_per_axis[0], info.nodes_per_axis[1], info.nodes_per_axis[2], info.band_bricks,
                      chk.vertices, chk.triangles, chk.components.size(), chk.components.size() == 1 ? "" : "s",
                      chk.closed_manifold() ? "closed manifold" : "NOT a closed manifold", chk.euler, chk.area,
                      chk.volume, chk.min_angle, 100 * chk.below_10 );
         if ( !chk.closed_manifold() || chk.degenerate ) {
            std::printf( "  open edges %ld, non-manifold edges %ld, misoriented edges %ld, non-manifold vertices %ld, "
                         "zero-area triangles %ld\n",
                         chk.open_edges, chk.nonmanifold_edges, chk.misoriented_edges, chk.nonmanifold_vertices,
                         chk.degenerate );
         }
         if ( chk.components.size() > 1 ) {
            for ( size_t c = 0; c < chk.components.size() && c < 8; ++c ) {
               std::printf( "  component %zu: %ld triangles, area %.4f, volume %.4f\n", c + 1,
                            chk.components[c].triangles, chk.components[c].area, chk.components[c].volume );
            }
         }
         if ( info.cavity_nodes ) {
            std::printf( "  %ld grid nodes in internal cavities %s\n", info.cavity_nodes, fill ? "filled" : "kept" );
         }
         if ( verbose ) {
            std::printf( "  %ld node and %ld root evaluations; seconds: bricks %.2f, nodes %.2f, fill %.2f, "
                         "vertices %.2f, triangles %.2f\n",
                         info.node_evals, info.root_evals, info.seconds_classify, info.seconds_nodes, info.seconds_fill,
                         info.seconds_vertices, info.seconds_triangles );
         }
         h2.push_back( hs[k] * hs[k] );
         areas.push_back( chk.area );
         vols.push_back( chk.volume );
      }

      double area = chk.area, volume = chk.volume;
      if ( hs.size() > 1 ) {
         double ba, bv;
         area = fit_line( h2, areas, ba );
         volume = fit_line( h2, vols, bv );
         std::printf( "extrapolated to h = 0 (a + b h^2, %zu spacings): area %.4f (b %.4f), volume %.4f (b %.4f)\n",
                      hs.size(), area, ba, volume, bv );
      }
      if ( want_msms ) {
         write_msms( base, name, mesh, model.atoms(), probe, h );
         std::printf( "wrote %s.vert, %s.face\n", base.c_str(), base.c_str() );
      }
      if ( want_c3p ) {
         write_c3p( base, name, mesh, model.atoms(), probe, area, volume );
         std::printf( "wrote %s.c3p, %s.c3v\n", base.c_str(), base.c_str() );
      }
      if ( !chk.closed_manifold() ) {
         return 2;
      }
   } catch ( const std::exception& e ) {
      fail( e.what() );
   }
   return 0;
}
