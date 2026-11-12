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
  -x h2,h3,..  also mesh at these spacings and extrapolate area and volume to h = 0 in h^2
  -o base      output base name (default: the input's name without .xyzr)
  -f list      msms, c3p or none (default msms,c3p)
  -k           keep internal cavities (default: fill them)
  -t n         threads
  -v           time each stage
```

The input has one sphere per line, `x y z r`. Further columns are ignored, and `#` starts a
comment. To reproduce msroll's surface, use the radii msroll would use.

For example, `sesmesh -p 1.5 -g 0.5 -x 0.3,0.25,0.2 protein.xyzr` writes `protein.c3p`/`.c3v`
and `protein.vert`/`.face` from a 0.5 Å grid. The area and volume in `protein.c3v` are
extrapolated from four spacings.

## Output

- `base.vert`, `base.face`: MSMS format (rcoal 3.0 reads it with `-s`). In the `.vert` file,
  each vertex has its outward normal, the input number of the nearest atom, and the type of
  patch the vertex is on: 1 contact, 2 reentrant, 3 toroidal.
- `base.c3p`, `base.c3v`: msroll's format, which rcoal reads.
  - `.c3p` holds the vertices with normals, the edges, and the triangles with their edges, all
    in component 1.
  - Line 13 of `.c3v` holds the area and volume that rcoal passes on to BEST. These are the
    h² extrapolated values when `-x` is given, otherwise the mesh's own.

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
5. **Area and volume** come from the mesh. Their error falls as h², so meshes at a few
   spacings extrapolate to the smooth surface's values (`-x`).

## Accuracy

`make test` runs the unit tests:
- f against the closed form for one sphere, and against a brute-force distance to a dense
  sampling of ∂U for random clusters;
- meshes of one, two and three spheres against the closed-form area and volume (h = 0.2 to
  0.1 Å, extrapolated in h²). One case pinches to a cusp;
- cavity filling;
- closed-manifold checks on random clusters.

The extrapolated area and volume of the analytic cases agree with the closed forms to about
10⁻⁵. In the cusp case the area agrees to 7×10⁻⁴: near a conical tip the error goes as
h² log h.

Lysozyme (PDB 1AKI) with msroll's 1001 spheres and a 1.5 Å probe, compared with msroll's
analytic SES (area 5633.206 Å², volume 23368.748 Å³). Every mesh from h = 0.5 to 0.05 Å is
one closed manifold with Euler characteristic 2.

| h (Å) | triangles | area | volume | wall time (64 threads) |
|---|---|---|---|---|
| 0.5 | 204,532 | −0.136% | −0.026% | 0.5 s |
| 0.25 | 821,140 | −0.035% | −0.007% | 0.7 s |
| 0.2 | 1,283,832 | −0.024% | −0.004% | 0.9 s |
| 0.1 | 5,139,668 | −0.006% | −0.001% | 2.5 s |
| 0.05 | 20,562,584 | −0.002% | −0.0003% | 9.3 s |
| h² fit, 0.2/0.15/0.1 | | −0.0003% | −0.00002% | |

ZENO 5.3.1 (10⁷ walks) on the 0.25, 0.2 and 0.15 Å meshes gives a translational diffusion
coefficient within 0.021% and an intrinsic viscosity within 0.084% of ZENO on a 0.25 Å
NanoShaper mesh. That is within about 1σ of the combined statistical error.

## Limitations

These are planned next:
- **No clean-up yet.** Vertices are not yet projected along ∇f, and slivers are not removed.
  About 11% of triangles have an angle under 10°.
- **Vertices near grid nodes** can sit about 2×10⁻⁵ Å apart. A tool that matches vertices by
  coordinates rounded to 6 significant digits sees those tiny triangles as degenerate.
- **rcoal.**
  - rcoal 2.0 as distributed stops on these meshes at its limit of 200 slender obtuse
    triangles.
  - With the fixes in ehb54/best-hydro#1 and ehb54/best-hydro#2, it reads them and writes
    closed meshes. It then stops coarsening at about 6,000 triangles for lysozyme.
- **Output formats.** There is no BEST-format or ZENO `.bod` output yet.

## References

- M. L. Connolly, J. Appl. Cryst. 16, 548–558 (1983): the solvent-excluded surface.
- W. E. Lorensen and H. E. Cline, Computer Graphics 21(4), 163–169 (1987): marching cubes. This
  code uses marching tetrahedra on the Freudenthal (Kuhn) subdivision instead.
- D. Xu and Y. Zhang, PLoS ONE 4, e8140 (2009): computing the SES as a distance field on a
  grid.

## Licence

MIT; see `LICENSE`.
