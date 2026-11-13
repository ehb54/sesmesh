# sesmesh

`sesmesh` triangulates the solvent-excluded surface (SES, the molecular surface of Richards and
Connolly) of a set of spheres. It makes closed surface meshes of macromolecules for
boundary-element hydrodynamics ([BEST](https://github.com/ehb54/best-hydro), through its mesh
coarsener rcoal) and for ZENO.

- It writes msroll's `.c3p`/`.c3v` format, which rcoal reads, and MSMS's `.vert`/`.face`.
- It is a free replacement for msroll, which is commercial, and MSMS, whose binary is free only
  for academic use.
- It is written from the published methods in C++17, needs no libraries, uses OpenMP when it
  is available, and is MIT-licensed.

## Building

```bash
make                   # bin/sesmesh, with OpenMP
make OPENMP=           # without OpenMP (one thread); needed with Apple's clang
make test              # unit tests, a few seconds
```

## Running

```
sesmesh [options] input.xyzr
  -p probe     probe radius in A (default 1.5)
  -g h         grid spacing in A (default 0.2)
  -n N         choose the grid spacing for a mesh of about N triangles (instead of -g)
  -x h2,h3,..  also mesh at these spacings and extrapolate area and volume to h = 0 in h^2
  -r           raw marching-tetrahedra mesh: skip the clean-up
  -c frac      clean-up: collapse edges shorter than frac times h (default 0.3)
  -o base      output base name (default: the input's name without .xyzr)
  -f list      msms, c3p, best, bod or none (default msms,c3p)
  -M mass      mass in Da for the .bod file
  -T temp      temperature in C for the .bod file (default 20)
  -V visc      solvent viscosity in cP for the .bod file (default 1.002)
  -k           keep internal cavities (default: fill them)
  -t n         threads
  -v           time each stage
```

The input has one sphere per line, `x y z r`. Further columns are ignored, and `#` starts a
comment. To reproduce msroll's surface, use the radii msroll would use.

Examples:
- `sesmesh -p 1.5 -n 25000 -x 0.3,0.25,0.2 protein.xyzr` writes `protein.c3p`/`.c3v` and
  `protein.vert`/`.face` with about 25,000 triangles, for rcoal. The area and volume in
  `protein.c3v` are extrapolated from four spacings.
- `sesmesh -p 1.5 -g 0.2 -f bod -M 14305 protein.xyzr` writes a fine mesh for ZENO.

With `-n`, sesmesh also meshes at a finer spacing and warns if the topology differs: a coarse
grid can open a narrow channel or close a thin neck.

## Output

- `base.vert`, `base.face`: MSMS format (rcoal 3.0 reads it with `-s`). In the `.vert` file,
  each vertex has its outward normal, the input number of the nearest atom, and the type of
  patch the vertex is on: 1 contact, 2 reentrant, 3 toroidal.
- `base.c3p`, `base.c3v`: msroll's format, which rcoal reads.
  - `.c3p` holds the vertices with normals, the edges, and the triangles with their edges, all
    in component 1.
  - Line 13 of `.c3v` holds the area and volume that rcoal passes on to BEST. These are the
    h² extrapolated values when `-x` is given, otherwise the mesh's own.

- `base.best`: BEST's triangle format. The first line is the number of triangles, then the
  area and volume (as for `.c3v`), then one line with the nine coordinates of each triangle.
- `base.bod`: ZENO's format: `UNITS A`, `TEMP`, `VISCOSITY`, `MASS` (with `-M`), then one
  `TRIANGLE` line per triangle.

Triangles run counter-clockwise seen from outside. Coordinates are written with 8 decimals.

## Method

1. **Exact surface function.** Let U be the union of the balls of radius r + p (the
   solvent-accessible solid). Then f(x) = p − dist(x, ∂U) inside U and p + dist(x, U)
   outside.
   - f is negative inside the SES and positive outside it, and its zero set is exactly the
     SES, cusps included.
   - The point of ∂U nearest to x lies on one of three features, each not inside any other
     ball:
     - a sphere patch (the radial projection of x),
     - a circle where two spheres meet (the nearest point of that circle),
     - a point where three spheres meet.
   - The exposed arcs of every circle are computed once; their ends are the exposed triple
     points. A spatial hash finds the spheres near a point.
2. **Narrow band.** The grid is cut into bricks of 8×8×8 nodes. f is 1-Lipschitz, so one
   value at a brick's centre shows whether the whole brick is clearly outside, clearly inside,
   or near the surface. Only bricks near the surface are evaluated node by node.
3. **Cavities.** A flood fill from the edge of the grid marks the outside that solvent can
   reach. It moves along the same 14 neighbour directions as the tetrahedron edges. Outside
   regions it does not reach are internal cavities, and they are counted as inside.
4. **Marching tetrahedra.** Each grid cube is split into the six Freudenthal tetrahedra
   around its main diagonal. This split matches between neighbouring cubes, so the cut is a
   closed, oriented 2-manifold. No triangle leaves its tetrahedron, so the mesh cannot
   intersect itself.
   - Each vertex is a root of the exact f on its edge (false position with the Illinois
     step), so it lies on the true surface.
   - A vertex is kept at least 10⁻⁴ of its edge from either end, so that a surface through a
     grid node gives no zero-area triangles.
   - The grid is shifted by irrational fractions of h, so that its nodes do not line up with
     the sphere centres.
5. **Clean-up** (skipped with `-r`). Marching tetrahedra leaves about a quarter of its
   triangles with an angle under 20°: needles where the surface passes close to a grid node,
   and caps. Three operations remove them, repeated a few times, and none moves a vertex off
   the exact surface:
   - **Edge collapse.** An edge shorter than 0.3 h becomes one vertex, at either end or at the
     midpoint projected onto the surface, whichever leaves the best angles.
     - The projection is one Newton step on f, exact while the nearest feature stays the same.
     - The link condition keeps the mesh a manifold of the same topology.
     - A collapse is made only if no triangle turns over and the smallest angle around the
       edge improves or stays above 20°.
   - **Edge flip.** The shared edge of two nearly coplanar triangles is swapped for the other
     diagonal when that raises their smallest angle.
   - **Tangential smoothing.** A vertex moves halfway to its neighbours' centroid within the
     tangent plane and is projected back onto the surface. The move is kept only if it raises
     the smallest angle around the vertex.
6. **Area and volume** come from the mesh. Their error falls as h², so meshes at a few
   spacings extrapolate to the smooth surface's values (`-x`). The extrapolation uses the
   meshes before clean-up, whose error follows h² most closely.

## Accuracy

`make test` runs the unit tests:
- f against the closed form for one sphere, and against a brute-force distance to a dense
  sampling of ∂U for random clusters;
- meshes of one, two and three spheres against the closed-form area and volume (h = 0.2 to
  0.1 Å, extrapolated in h²). One case pinches to a cusp;
- cavity filling;
- closed-manifold checks on random clusters;
- the clean-up: every cleaned mesh keeps its topology, its area and volume are no further from
  the closed form than the raw mesh's, and it keeps fewer than 1% of triangles under 20°;
- the writers: what they write reads back.

The extrapolated area and volume of the analytic cases agree with the closed forms to about
10⁻⁵. In the cusp case the area agrees to 7×10⁻⁴: near a conical tip the error goes as
h² log h.

Lysozyme (PDB 1AKI) with msroll's 1001 spheres and a 1.5 Å probe, compared with msroll's
analytic SES (area 5633.206 Å², volume 23368.748 Å³). Every mesh from h = 0.5 to 0.05 Å is
one closed manifold with Euler characteristic 2. Before clean-up:

| h (Å) | triangles | area | volume | wall time (64 threads) |
|---|---|---|---|---|
| 0.5 | 204,532 | −0.136% | −0.026% | 0.5 s |
| 0.25 | 821,140 | −0.035% | −0.007% | 0.7 s |
| 0.2 | 1,283,832 | −0.024% | −0.004% | 0.9 s |
| 0.1 | 5,139,668 | −0.006% | −0.001% | 2.5 s |
| 0.05 | 20,562,584 | −0.002% | −0.0003% | 9.3 s |
| h² fit, 0.2/0.15/0.1 | | −0.0003% | −0.00002% | |

The clean-up removes about a third of the triangles. Triangles with an angle under 20° fall
from 25% to under 0.05%, and area and volume change by about 10⁻⁴ or less:

| h (Å) | triangles | under 20° | area | volume | clean-up time (one thread) |
|---|---|---|---|---|---|
| 0.5 | 141,472 | 37 | −0.148% | −0.020% | 3.6 s |
| 0.25 | 566,178 | 73 | −0.038% | −0.005% | 14 s |
| 0.2 | 884,082 | 76 | −0.026% | −0.003% | 22 s |

ZENO 5.3.1 (10⁷ walks and interior samples) was run on the full meshes of three programs,
without coarsening, and each program's results were extrapolated in 1/N:

| mesher | meshes | triangles | D (10⁻⁶ cm²/s) | [η] (cm³/g) |
|---|---|---|---|---|
| sesmesh, before clean-up | 13 | 35k–3.3M | 1.11515 ± 0.00009 | 3.1777 ± 0.0005 |
| NanoShaper | 6 | 38k–1.08M | 1.11528 ± 0.00012 | 3.1768 ± 0.0007 |
| msroll | 12 | 17k–43k | 1.11544 ± 0.00028 | 3.1757 ± 0.0033 |

They agree within 0.03% in D and 0.06% in [η]. ZENO on the cleaned 0.5–0.2 Å meshes gives D
1.11495–1.11532 and [η] 3.1775–3.1815, within its statistical error (0.02% and 0.055%) of
these limits. msroll's error is larger because it wrote no lysozyme mesh finer than 43k
triangles.

## Limitations

These are planned next:
- **Speed.** The clean-up runs on one thread, so for a fine mesh it takes far longer than the
  meshing itself.
- **A few slivers remain.** Fewer than 0.1% of triangles are under 20°, but the smallest angle
  can still be about 1°.
- **rcoal.**
  - rcoal 2.0 as distributed reads cleaned meshes. However, some of its outputs have
    non-manifold or misoriented edges.
  - With the fixes in ehb54/best-hydro#1 and ehb54/best-hydro#2 its outputs are closed. It
    coarsens a cleaned mesh by at most about 12 times, so its outputs follow the input. For
    lysozyme with SOMO's targets (2,152–3,588 triangles):

    | `-n` | rcoal outputs |
    |---|---|
    | 12,000 | 1,944–2,748 |
    | 25,000 | 2,674–3,566 |
    | 40,000 | 3,448–4,404 |

  - BEST through rcoal, extrapolated linearly in 1/N:
    - At a 30,000-triangle input, sesmesh agrees with msroll at the same input size within
      0.03% in D and [η].
    - The extrapolation over rcoal's narrow range of N carries a systematic of a few tenths of
      a percent in D, and up to about 1% in [η]. It depends on the input mesh, for msroll's
      meshes as well. ZENO on a full mesh does not have it.

## References

- M. L. Connolly, J. Appl. Cryst. 16, 548–558 (1983): the solvent-excluded surface.
- W. E. Lorensen and H. E. Cline, Computer Graphics 21(4), 163–169 (1987): marching cubes. This
  code uses marching tetrahedra on the Freudenthal (Kuhn) subdivision instead.
- D. Xu and Y. Zhang, PLoS ONE 4, e8140 (2009): computing the SES as a distance field on a
  grid.

## Licence

MIT; see `LICENSE`.
