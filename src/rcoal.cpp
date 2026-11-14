// rcoal: makes a series of BEST input meshes with Coalesce.  A C++ drop-in for rcoal 2.0
// (rcoal.f) and rcoal 3.0 (rcoal3.f) in ehb54/best-hydro, as fixed there in PRs #1 and #2: the
// same command line, the same output files, byte for byte, and the same stdout, including
// the "Output files written:" lines that SOMO reads.
//
// Copyright (c) 1998-2025 Sergio R. Aragon (Coalesce and rcoal, in ehb54/best-hydro)
// Copyright (c) 2026 Emre H. Brookes (this port, and the fixes of best-hydro PRs #1 and #2)
// MIT licence: see LICENSE.
//
//   rcoal -f name.c3p -nmax N -nmin N -n K   K meshes from nmin up to nmax triangles
//   rcoal -f name.c3p -e edgfrac             one mesh
//   -s        MSMS input: give name.face; name.vert and name.msms.log are read too
//   -search   find e for each target by bisection, instead of rcoal's two-run fit
//   -v3       rcoal 3.0's messages (the default when run as rcoal3, and with -s)
//
// Differences from the Fortran:
//  - Options are matched exactly.  The Fortran matches -n, -e and -s as prefixes, so a -n
//    given before -nmax or -nmin was overwritten by their values.
//  - No limit on the length of file names or on the number of meshes.
//  - -search is new.
#include "coalesce.h"
#include "fortio.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sesmesh::rcoal;

namespace {

bool v3 = false;                 // rcoal 3.0's messages
std::string prog = "rcoal";      // the name the program was run as, for messages
std::vector<std::string> outf;   // the files written, in order

void out( const ListRecord& r ) { std::cout << r.str() << '\n'; }

// The output base name: the input name up to its extension (.c3p, or .fac of .face).
std::string base_name( const std::string& in, bool sanner ) {
   std::string ext = sanner ? ".fac" : ".c3p";
   size_t p = in.find( ext );
   if ( p == std::string::npos ) {
      throw std::runtime_error( "input file name must end in " + std::string( sanner ? ".face" : ".c3p" ) + ": " + in );
   }
   return in.substr( 0, p );
}

// msroll's .c3p and, from its .c3v, the volume and molecular area of component 1.
Surface read_c3p( const std::string& name ) {
   Surface s;
   s.input_name = name;
   ListReader r( name );
   std::vector<double> v;
   r.read( v, 3 );
   int nvert = int( v[0] ), nedge = int( v[1] ), ntriang = int( v[2] );
   s.x.resize( nvert + 1 );
   for ( int n = 1; n <= nvert; ++n ) {
      r.read( v, 12 );
      s.x[n] = { v[0], v[1], v[2] };
   }
   s.edges.resize( nedge );
   for ( int n = 0; n < nedge; ++n ) {
      r.read( v, 5 );
      s.edges[n] = { int( v[0] ), int( v[1] ) };
      for ( int k = 0; k < 2; ++k ) {
         if ( s.edges[n][k] < 1 || s.edges[n][k] > nvert ) {
            throw std::runtime_error( "edge with a vertex number out of range in " + name );
         }
      }
   }
   s.tri.resize( ntriang );
   s.comp.resize( ntriang );
   double totarea = 0.0;
   for ( int n = 0; n < ntriang; ++n ) {
      r.read( v, 9 );
      s.tri[n] = { int( v[3] ), int( v[4] ), int( v[5] ) };
      s.comp[n] = int( v[6] );
      for ( int k = 0; k < 3; ++k ) {
         if ( s.tri[n][k] < 1 || s.tri[n][k] > nvert ) {
            throw std::runtime_error( "triangle with a vertex number out of range in " + name );
         }
      }
      totarea = totarea + tri_area( s.x[s.tri[n][0]], s.x[s.tri[n][1]], s.x[s.tri[n][2]] );
   }
   s.begin_area = totarea;
   // .c3v, line 13: component, centre (3), volume, molecular area, accessible area
   std::string c3v = name.substr( 0, name.size() - 4 ) + ".c3v";
   ListReader rv( c3v );
   for ( int n = 0; n < 12; ++n ) {
      rv.skip_line();
   }
   rv.read( v, 7 );
   s.true_vol = v[4];
   s.true_area = v[5];
   return s;
}

// MSMS's .vert and .face, and from its .msms.log the analytical area and numerical volume.
Surface read_msms( const std::string& face ) {
   Surface s;
   s.input_name = face;
   s.sanner = true;
   std::string base = base_name( face, true );
   std::vector<double> v;
   ListReader rv( base + ".vert" );
   rv.skip_line();
   rv.skip_line();
   rv.read( v, 4 );
   int nvert = int( v[0] );
   s.x.resize( nvert + 1 );
   for ( int n = 1; n <= nvert; ++n ) {
      rv.read( v, 8 );
      s.x[n] = { v[0], v[1], v[2] };
   }
   ListReader rf( face );
   rf.skip_line();
   rf.skip_line();
   rf.read( v, 4 );
   int ntriang = int( v[0] );
   s.tri.resize( ntriang );
   s.comp.assign( ntriang, 1 );
   double totarea = 0.0;
   for ( int n = 0; n < ntriang; ++n ) {
      rf.read( v, 5 );
      s.tri[n] = { int( v[0] ), int( v[1] ), int( v[2] ) };
      for ( int k = 0; k < 3; ++k ) {
         if ( s.tri[n][k] < 1 || s.tri[n][k] > nvert ) {
            throw std::runtime_error( "face with a vertex number out of range in " + face );
         }
      }
      totarea = totarea + tri_area( s.x[s.tri[n][0]], s.x[s.tri[n][1]], s.x[s.tri[n][2]] );
   }
   s.begin_area = totarea;
   // The unique edges, in order of first appearance, as the Fortran builds them.
   std::vector<std::vector<int>> adj( nvert + 1 );
   for ( const T3& t : s.tri ) {
      for ( int m = 0; m < 3; ++m ) {
         int a = t[m], b = t[( m + 1 ) % 3];
         int lo = std::min( a, b ), hi = std::max( a, b );
         bool dup = false;
         for ( int w : adj[lo] ) {
            if ( w == hi ) {
               dup = true;
               break;
            }
         }
         if ( !dup ) {
            adj[lo].push_back( hi );
            s.edges.push_back( { a, b } );
         }
      }
   }
   // name.msms.log: the analytical SES area, and the numerical volume and area
   std::ifstream lg( base + ".msms.log" );
   std::string line;
   auto toks = []( const std::string& l ) {
      std::vector<std::string> t;
      size_t i = 0;
      while ( i < l.size() ) {
         while ( i < l.size() && ( l[i] == ' ' || l[i] == '\t' || l[i] == ',' ) ) {
            ++i;
         }
         size_t j = i;
         while ( j < l.size() && l[j] != ' ' && l[j] != '\t' && l[j] != ',' ) {
            ++j;
         }
         if ( j > i ) {
            t.push_back( l.substr( i, j - i ) );
         }
         i = j;
      }
      return t;
   };
   while ( std::getline( lg, line ) ) {
      auto t = toks( line );
      if ( t.size() >= 2 && t[0].compare( 0, 10, "ANALYTICAL" ) == 0 && t[1].compare( 0, 7, "SURFACE" ) == 0 ) {
         std::getline( lg, line );
         std::getline( lg, line );
         auto u = toks( line );
         if ( u.size() >= 6 ) {
            s.true_area = std::strtod( u[5].c_str(), nullptr );
         }
      } else if ( !t.empty() && t[0].compare( 0, 9, "NUMERICAL" ) == 0 ) {
         std::getline( lg, line );
         std::getline( lg, line );
         auto u = toks( line );
         if ( u.size() >= 4 ) {
            s.poly_vol = std::strtod( u[2].c_str(), nullptr );
            s.poly_area = std::strtod( u[3].c_str(), nullptr );
         }
      }
   }
   s.true_vol = s.poly_vol;
   return s;
}

// The area histogram: subroutines bin and dstrb, 1-based as in the Fortran.
void bin( double e, double s, std::vector<double>& bins, int numbin ) {
   bins.assign( numbin + 1, 0.0 );
   for ( int j = 1; j <= ( numbin - 1 ) / 2; ++j ) {
      bins[j] = e - double( ( numbin + 1 ) / 2 - j ) * s * e / 2.5;
      bins[j + ( numbin + 1 ) / 2] = e + double( j ) * s * e / 2.5;
   }
   bins[( numbin + 1 ) / 2] = e;
}

void dstrb( double z, std::vector<int>& nd, const std::vector<double>& bins, int numbin ) {
   for ( int j = 1; j <= ( numbin - 2 ) / 2; ++j ) {
      if ( z <= bins[-j + 1 + ( numbin / 2 )] && z > bins[( numbin / 2 ) - j] ) {
         nd[( numbin / 2 ) + 1 - j]++;
         return;
      }
      if ( z > bins[( numbin / 2 ) - 1 + j] && z <= bins[( numbin / 2 ) + j] ) {
         nd[( numbin / 2 ) - 1 + j]++;
         return;
      }
   }
   if ( z <= bins[1] ) {
      nd[1]++;
   } else if ( z > bins[numbin] ) {
      nd[numbin]++;
   } else if ( z <= bins[numbin] && z > bins[numbin - 1] ) {
      nd[numbin - 1]++;
   }
}

// The summary that scoalesce prints after each triangulation (his: as written to the .his
// file, whose labels differ slightly).
void summary( std::ostream& os, const Surface& s, const Stats& st, bool his ) {
   const double ver2 = 3.1, ver3 = 3.5;
   std::vector<ListRecord> L;
   auto add = [&]( const ListRecord& r ) { L.push_back( r ); };
   add( ListRecord().c( "COALESCE PROCESS RESULTS -Version= " ).r( v3 ? ver3 : ver2 ) );
   add( ListRecord().c( " " ) );
   add( ListRecord().c( s.sanner ? "Input File (Sanner): " : "Input File (Connolly): " ).c80( s.input_name ) );
   add( ListRecord().c( " Edge fraction = " ).r( st.edgfrac ) );
   add( ListRecord().c( " " ) );
   add( ListRecord().c( " Beginning polygon Area =  " ).r( s.begin_area ) );
   add( ListRecord().c( his ? "Initial number of triangle records= " : "Number of triangle records, pass 1= " ).i( st.ntriang ) );
   add( ListRecord().c( "Number of degenerate triangles, pass 1= " ).i( st.idegen ) );
   add( ListRecord().c( "Number of small triangles deleted, pass 1= " ).i( st.idelete1 ) );
   add( ListRecord().c( "Number of polygons generated, pass 1= " ).i( st.npolygon1 ) );
   add( ListRecord().c( his ? "total area, pass 1 = " : "Total area, pass 1 = " ).r( st.totarea1 ) );
   add( ListRecord().c( " " ) );
   add( ListRecord().c( "Number of small triangles deleted, pass 2= " ).i( st.idelete ) );
   add( ListRecord().c( "Number of polygons generated, pass 2= " ).i( st.npoly ) );
   add( ListRecord().c( his ? "total area, pass 2 = " : "Total area, pass 2 = " ).r( st.totarea ) );
   add( ListRecord().c( "Average Triangle Area (A^2)= " ).r( st.avarea ) );
   add( ListRecord().c( "Average edgelength (A)= " ).r( st.avelength / st.nedge ) );
   add( ListRecord().c( "Cut-off length (A) = " ).r( st.cutofflength ) );
   add( ListRecord().c( "Found" ).i( st.nobtuse ).c( " slender obtuse triangles, processed" ).i( st.icon ) );
   if ( !his ) {
      add( ListRecord().c( "Collapses skipped by the topology check: " ).i( st.nlinkskip ) );
      add( ListRecord().c( "Flips skipped by the topology check: " ).i( st.nflipskip ) );
   } else {
      add( ListRecord().c( " " ) );
   }
   for ( auto& r : L ) {
      os << r.str() << '\n';
   }
}

// Writes one triangulation and its .his file, unless a triangle has two vertices at the same
// position as written; nzero receives the number of such triangles.
void write_result( const Surface& s, const Result& R, double edgfrac, int& nzero ) {
   const Stats& st = R.st;
   std::string base = base_name( s.input_name, s.sanner );
   // the coordinates as written
   std::vector<std::string> cv( 3 * st.npoly );
   nzero = 0;
   for ( int m = 0; m < st.npoly; ++m ) {
      for ( int k = 0; k < 3; ++k ) {
         const V3& p = s.x[R.tri[m][k]];
         cv[3 * m + k] = fort_e( p[0], 17, 8 ) + fort_e( p[1], 17, 8 ) + fort_e( p[2], 17, 8 );
      }
      if ( cv[3 * m] == cv[3 * m + 1] || cv[3 * m + 1] == cv[3 * m + 2] || cv[3 * m] == cv[3 * m + 2] ) {
         ++nzero;
      }
   }
   if ( nzero > 0 ) {
      out( ListRecord().c( "Triangles with two vertices at the same position: " ).i( nzero ).c( "; not written" ) );
      return;
   }
   char num[16];
   std::snprintf( num, sizeof num, "%05d", st.npoly % 100000 );
   std::string output = base + "_" + num;
   std::string outputhis = output + ".his";

   {
      std::ofstream h( outputhis );
      summary( h, s, st, true );
      h << ListRecord().c( "Histogram bins and count follows:" ).str() << '\n';
      std::vector<double> bins, bins2;
      std::vector<int> nd1( 51, 0 ), nd2( 51, 0 );
      int numbin = 12;
      bin( 0.1 * st.avarea, 0.5, bins, numbin + 1 );
      for ( int m = 0; m < st.npoly; ++m ) {
         if ( R.area[m] <= 0.2 * st.avarea ) {
            dstrb( R.area[m], nd1, bins, numbin );
         }
      }
      h << ListRecord().c( "Area distribution around" ).r( .1 * st.avarea ).c( " A^2" ).str() << '\n';
      int nsmallarea = 0;
      for ( int n = 2; n <= numbin; ++n ) {
         h << ListRecord().r( bins[n] ).i( nd1[n] ).str() << '\n';
         nsmallarea += nd1[n];
      }
      if ( nsmallarea > 0 ) {
         out( ListRecord().c( " " ) );
         out( ListRecord().c( "************************" ) );
         out( ListRecord().c( "N =  " ).i( st.npoly ) );
         out( ListRecord().c( "ALERT! " ).i( nsmallarea ).c( " Triangles with Area< 20% average area Remain!" ) );
         out( ListRecord().c( "Review histogram." ) );
         out( ListRecord().c( "Edgefrac = " ).r( edgfrac ) );
         out( ListRecord().c( "************************" ) );
      }
      numbin = 22;
      bin( st.avarea, 0.2, bins2, numbin + 1 );
      for ( int m = 0; m < st.npoly; ++m ) {
         if ( R.area[m] > 0.2 * st.avarea ) {
            dstrb( R.area[m], nd2, bins2, numbin );
         }
      }
      h << ListRecord().c( "Area distribution around" ).r( st.avarea ).c( "A^2" ).str() << '\n';
      for ( int n = 2; n <= numbin; ++n ) {
         h << ListRecord().r( bins2[n] ).i( nd2[n] ).str() << '\n';
      }
      if ( !h ) {
         throw std::runtime_error( "cannot write " + outputhis );
      }
   }
   {
      std::ofstream o( output );
      o << ListRecord().i( st.npoly ).r( s.true_area ).r( s.true_vol ).str() << '\n';
      for ( int m = 0; m < st.npoly; ++m ) {
         o << cv[3 * m] << cv[3 * m + 1] << cv[3 * m + 2] << '\n';
      }
      if ( !o ) {
         throw std::runtime_error( "cannot write " + output );
      }
   }
   outf.push_back( output );
   outf.push_back( outputhis );
}

// One call of subroutine scoalesce: coalesce with edge fraction edgfrac, print what the
// Fortran prints, and, if write, write the result.  Returns the number of triangles.
int scoalesce( const Surface& s, double edgfrac, bool write, int& nzero ) {
   nzero = 0;
   out( ListRecord().c( " " ) );
   out( ListRecord().c( "Begin Coalesce Processing" ) );
   out( ListRecord().c( " " ) );
   std::string base = base_name( s.input_name, s.sanner );
   if ( v3 ) {
      out( ListRecord().c( "basename = " ).c80( base ) );
   }
   if ( s.sanner ) {
      // scoalesce3 reads the MSMS files again on every call, and says so.
      out( ListRecord().c( "Number triangles = " ).i( long( s.tri.size() ) ) );
      out( ListRecord().i( long( s.edges.size() ) ).c( " unique edges found" ) );
      out( ListRecord().c( "Input file = " ).c80( base + ".msms.log" ) );
      out( ListRecord().c( "Areas: True, Poly = " ).r( s.true_area ).r( s.poly_area ) );
      out( ListRecord().c( "Volumes: True, Poly =" ).r( s.true_vol ).r( s.poly_vol ) );
   }
   std::vector<std::string> log;
   Result R = coalesce( s, edgfrac, &log );
   for ( auto& l : log ) {
      std::cout << l << '\n';
   }
   out( ListRecord().c( "  " ) );
   summary( std::cout, s, R.st, false );
   if ( write ) {
      write_result( s, R, edgfrac, nzero );
   }
   return R.st.npoly;
}

void list_files() {
   for ( size_t i = 0; i + 1 < outf.size(); i += 2 ) {
      out( ListRecord().c( " " ) );
      out( ListRecord().c( " Output files written:" ) );
      out( ListRecord().c( "                      " ).c80( outf[i] ) );
      out( ListRecord().c( "                      " ).c80( outf[i + 1] ) );
      out( ListRecord().c( " " ) );
   }
}

// The target triangle counts, as rcoal computes and prints them.
std::vector<int> targets( int nmax, int nmin, int& n ) {
   int ndelta = 0;
   if ( nmax == nmin ) {
      n = 1;
   } else if ( n > 1 ) {
      ndelta = ( nmax - nmin ) / ( n - 1 );
   }
   std::vector<int> nt( n );
   out( ListRecord().c( " " ) );
   out( ListRecord().c( "Target values of n" ) );
   for ( int j = 0; j < n; ++j ) {
      nt[j] = nmin + j * ndelta;
      out( ListRecord().c( "             " ).i( nt[j] ) );
   }
   out( ListRecord().c( " " ) );
   return nt;
}

void usage() {
   if ( v3 ) {
      out( ListRecord().c( "Usage: rcoal -f filename -nmax val -nmin val " ) );
      out( ListRecord().c( "-n numb [-s]" ) );
      out( ListRecord().c( "Usage: rcoal -f filename -e val [-s]" ) );
   } else {
      out( ListRecord().c( "Usage: rcoal -f filename -nmax val -nmin val " ) );
      out( ListRecord().c( "-n numb " ) );
      out( ListRecord().c( "Usage: rcoal -f filename -e val " ) );
   }
   out( ListRecord().c( "Also: -search (find e for each target by bisection), -v3 (rcoal 3.0's messages)" ) );
}

}  // namespace

int main( int argc, char** argv ) {
   prog = argv[0];
   size_t slash = prog.find_last_of( "/\\" );
   if ( slash != std::string::npos ) {
      prog = prog.substr( slash + 1 );
   }
   if ( prog.compare( 0, 6, "rcoal3" ) == 0 ) {
      v3 = true;
   }
   if ( argc < 3 ) {
      usage();
      return 0;
   }
   std::string file;
   int nmax = -1, nmin = -1, n = -1;
   double edgefrac = -1.0;
   bool sanner = false, search = false;
   for ( int k = 1; k < argc; ++k ) {
      std::string a = argv[k];
      std::string next = k + 1 < argc ? argv[k + 1] : "";
      if ( a == "-f" ) {
         file = next;
      } else if ( a == "-nmax" ) {
         nmax = std::atoi( next.c_str() );
      } else if ( a == "-nmin" ) {
         nmin = std::atoi( next.c_str() );
      } else if ( a == "-n" ) {
         n = std::atoi( next.c_str() );
      } else if ( a == "-e" ) {
         edgefrac = std::atof( next.c_str() );
      } else if ( a == "-s" ) {
         sanner = true;
      } else if ( a == "-v3" ) {
         v3 = true;
      } else if ( a == "-search" ) {
         search = true;
      }
   }
   if ( sanner ) {
      v3 = true;
   }
   out( ListRecord().c( v3 ? "Rcoal Ver. 3.0" : "Rcoal Ver. 2.0" ) );
   out( ListRecord().c( "SFSU, Aragon Lab." ) );
   int nzero = 0;
   try {
      Surface s = sanner ? read_msms( file ) : read_c3p( file );
      int ntriang = int( s.tri.size() );
      if ( edgefrac > 0.0 ) {
         int n3 = scoalesce( s, edgefrac, true, nzero );
         if ( nzero == 0 ) {
            out( ListRecord().c( "Actual coalesce n" ).i( n3 ) );
         }
      } else if ( nmax == -1 || nmin <= 3 || n == 0 || n == -1 ) {
         out( ListRecord().c( "Check input parameters: n=0?, nmin <3,n,nmax missing?" ) );
      } else if ( nmax < nmin ) {
         out( ListRecord().c( "nmax is less than nmin" ) );
      } else if ( ntriang < nmax ) {
         out( ListRecord().c( "Nmax must be less than" ).i( ntriang + 1 ).c( "in this instance." ) );
      } else {
         std::vector<int> nt = targets( nmax, nmin, n );
         double p = 0, a = 0;
         if ( !search ) {
            // N(e) = ntriang / (1 + a e^p), fitted to runs at e = 0.7 and 1.0
            int n1 = scoalesce( s, 0.7, false, nzero );
            if ( v3 ) {
               out( ListRecord().c( "Just finished first scoalesce call" ) );
            }
            int n2 = scoalesce( s, 1.0, false, nzero );
            out( ListRecord().c( "initial n " ).i( n1 ).i( n2 ) );
            double y2 = std::log( double( ntriang ) / double( n2 ) - 1.0 );
            double y1 = std::log( double( ntriang ) / double( n1 ) - 1.0 );
            double x2 = std::log( 1.0 );
            double x1 = std::log( 0.7 );
            p = ( y2 - y1 ) / ( x2 - x1 );
            a = std::exp( y2 - p * x2 );
            out( ListRecord().c( " p = " ).r( p ) );
            out( ListRecord().c( " a = " ).r( a ) );
         }
         for ( int j = 0; j < n; ++j ) {
            double e = search ? find_edgfrac( s, nt[j] ) : 1.000 * std::pow( ( double( ntriang ) / double( nt[j] ) - 1.0 ) / a, 1 / p );
            out( ListRecord().c( " e = " ).r( e ) );
            int nn = scoalesce( s, e, true, nzero );
            if ( nzero > 0 ) {
               break;
            }
            out( ListRecord().c( "Target n         " ).i( nt[j] ) );
            out( ListRecord().c( "Actual coalesce n" ).i( nn ) );
            if ( !search ) {
               out( ListRecord().c( " " ) );
            }
         }
      }
   } catch ( const std::exception& e ) {
      std::cerr << prog << ": " << e.what() << '\n';
      return 2;
   }
   if ( nzero > 0 ) {
      for ( const std::string& f : outf ) {
         std::remove( f.c_str() );
      }
      std::cerr << " rcoal: a triangulation has triangles of zero area\n"
                << " (two vertices at the same position), which BEST\n"
                << " cannot use.  No triangulation written;\n"
                << " try a coarser input triangulation.\n"
                << "STOP 1\n";
      return 1;
   }
   list_files();
   return 0;
}
