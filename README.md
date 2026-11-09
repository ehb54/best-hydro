# BEST and rcoal

Boundary-element hydrodynamics for rigid macromolecules, by Sergio R. Aragon
(San Francisco State University).

- **BEST** computes the transport properties of a rigid body of any shape from a
  triangulation of its surface. It gives the translational and rotational diffusion
  tensors and their coupling, the centers of resistance and diffusion, and the
  intrinsic viscosity. It solves the Stokes equations with stick boundary conditions
  by the boundary element method of Youngren and Acrivos (J. Fluid Mech. 69, 377
  (1975)).
- **rcoal** prepares the triangulations. Its Coalesce routine removes the slender and
  very small triangles that molecular-surface programs produce, and coarsens a
  triangulation using only its original vertices. rcoal runs it to make a series of
  triangulations over a chosen range of triangle counts, so that BEST's results can be
  extrapolated linearly in 1/N to an infinite number of triangles.

The sources here are as the author sent them. Fixes and improvements come as pull
requests.

## Contents

| Directory | |
|---|---|
| `best/` | BEST 4.0 (April 2025), with the author's PGI make file and `best.contents` |
| `rcoal/` | rcoal 2.0 (2007), for MSROLL triangulations |
| `rcoal3/` | rcoal 3.0 (2012), which also reads MSMS triangulations |
| `doc/` | the BEST Hydrodynamics Package manual, version 3.9 (2014) |
| `examples/` | a sample BEST input, the rcoal 2.0 test case, and an rcoal 3.0 example (crambin, MSMS) |

`index.f` is a sort routine used only by Coalesce. It is in each directory, as
distributed.

## Building

You need gfortran and LAPACK/BLAS.

```bash
make                                  # Linux, reference LAPACK/BLAS
make LAPACK="-lopenblas"              # or another LAPACK/BLAS
make LAPACK="-framework Accelerate"   # macOS
```

This builds `bin/best`, `bin/rcoal` and `bin/rcoal3`. `make test` (with the same
`LAPACK=...`) runs the examples and checks the results; it needs Perl. The author compiled with PGI and
`-r8`, so that every REAL is 8 bytes; the Makefile uses gfortran's equivalent,
`-fdefault-real-8 -fdefault-double-8`. He recommends optimization no higher than `-O2`.

BEST stores a dense 3N x 3N matrix for N triangles: about 0.067 N² GB, N in
thousands. Its run time grows as N³.

## Running

The input is a triangulated molecular surface. Connolly's MSROLL is not included: it
is a commercial program. rcoal 3.0 also reads triangulations from Sanner's MSMS. The
manual in `doc/` describes each step and its options; in short:

```bash
msroll -m protein.pdb -p 1.5 -r radii -t protein.c3p -v protein.c3v -f 0.7
rcoal -f protein.c3p -nmax 7000 -nmin 2500 -n 4    # reads protein.c3v too
best -f protein_02500 -v -mw 14300
```

rcoal writes each triangulation to `<name>_<N>` (with a histogram in `<name>_<N>.his`),
with the molecule's exact area and volume on the first line. BEST writes
`<input><flags>.be` and a `.log`. By default BEST rescales the triangle areas to the
exact surface area. This area correction also regularizes the integral equation, so the
manual advises against turning it off (`-na`).

To try the examples, copy them to a scratch directory first, since the programs write
their output next to their input:

- `examples/best/2eu7_02562`: `best -f 2eu7_02562 -v -mw 30000` takes a minute or less.
- `examples/rcoal`: the `.his` files give the edge fraction of each output, so for
  example `rcoal -f test.c3p -e 0.8` reproduces `test_05770`.
- `examples/rcoal3`: `rcoal3 -f 1crn.face -nmax 2500 -nmin 1500 -n 3 -s`.

## References

- S. R. Aragon, "A precise boundary element method for macromolecular transport
  properties", J. Comput. Chem. 25, 1191-1205 (2004).
- S. R. Aragon and D. K. Hahn, "Precise boundary element computation of protein
  transport properties: diffusion tensors, specific volume, and hydration",
  Biophys. J. 91, 1591-1603 (2006).
- D. K. Hahn and S. R. Aragon, J. Chem. Theory Comput. 2, 1416-1428 (2006), for the
  intrinsic viscosity.

## Credits

Written by Sergio R. Aragon. David Hahn programmed the viscosity calculation and many
enhancements. Tilman Rosales contributed to the programs and to the hydration strategy
for proteins, and Martin Perez to the one for DNA. Stuart Allison gave the initial push. Jerome Fung noticed the
translation-rotation coupling error fixed in version 4.0. The work was supported by the
NIH MBRS SCORE program.

## License

Released by the author under the MIT License; see `LICENSE`.
