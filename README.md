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

The tag [`v4.0-aragon`](https://github.com/ehb54/best-hydro/tree/v4.0-aragon) holds the
programs as the author released them. The main branch carries the fixes and improvements
made since; see [Versions](#versions).

## Contents

| Directory | |
|---|---|
| `best/` | BEST 4.0 (April 2025), with the author's PGI make file and `best.contents` |
| `hmatrix/` | the H-matrix solver BEST uses above 10,000 triangles (C++17); see [Large meshes](#large-meshes) |
| `rcoal/` | rcoal 2.0 (2007), for MSROLL triangulations |
| `rcoal3/` | rcoal 3.0 (2012), which also reads MSMS triangulations |
| `doc/` | the BEST Hydrodynamics Package manual, version 3.9 (2014) |
| `examples/` | two sample BEST inputs, the rcoal 2.0 test case, and an rcoal 3.0 example (crambin, MSMS) |

`index.f` is a sort routine used only by Coalesce. It is in each directory, as
distributed.

## Building

You need gfortran, g++ (C++17; GCC 8 or later) and LAPACK/BLAS.

```bash
make                                                # Linux, reference LAPACK/BLAS
make LAPACK="-lopenblas"                            # or another LAPACK/BLAS
make LAPACK="-framework Accelerate" CXXLIBS=-lc++   # macOS, where g++ is clang
make HMATRIX=                                       # BEST without the H-matrix solver (no C++)
```

This builds `bin/best`, `bin/rcoal` and `bin/rcoal3`. `make test` (with the same
`LAPACK=...`) runs the examples and checks the results; it needs Perl. The author compiled with PGI and
`-r8`, so that every REAL is 8 bytes; the Makefile uses gfortran's equivalent,
`-fdefault-real-8 -fdefault-double-8`. He recommends optimization no higher than `-O2`.
The Makefile also turns off fused multiply-adds (`-ffp-contract=off`), so the results do not
depend on whether the compiler targets a processor that has them.

Up to 10,000 triangles BEST stores a dense 3N x 3N matrix for N triangles: about
0.067 N² GB, N in thousands, and its run time grows as N³. Above that it uses an
H-matrix solver, whose memory and time grow about as N log N; see
[Large meshes](#large-meshes).

## Running

The input is a triangulated molecular surface. Connolly's MSROLL is not included: it
is a commercial program. rcoal 3.0 also reads triangulations from Sanner's MSMS.
[sesmesh](https://github.com/ehb54/sesmesh), a free (MIT) mesher, writes both formats
from the same spheres and probe radius. rcoal coarsens its meshes into closed ones. For
lysozyme, with rcoal starting from meshes of the same size (30,000 triangles), BEST gives
the same D and [η] from sesmesh's meshes as from MSROLL's, within 0.03%; the sesmesh
README has the comparisons. The manual in `doc/` describes each step and its options; in
short:

```bash
msroll -m protein.pdb -p 1.5 -r radii -t protein.c3p -v protein.c3v -f 0.7
rcoal -f protein.c3p -nmax 7000 -nmin 2500 -n 4    # reads protein.c3v too
best -f protein_02500 -v -mw 14300
```

rcoal writes each triangulation to `<name>_<N>` (with a histogram in `<name>_<N>.his`),
with the molecule's exact area and volume on the first line. BEST writes
`<input><flags>.be` and a `.log`. By default BEST rescales the triangle areas to the
exact surface area, and the manual advises against turning this off (`-na`). The
rescaling also kept BEST's matrix away from a near-singularity. BEST now removes that
near-singularity itself, so `-na`, and meshes whose own area is exact (such as sesmesh's),
give sound results too.

A triangulation may hold several separate closed surfaces, for example a molecule and a
detached atom. Each surface brings its own near-singularity, and BEST removes each one. It
finds the surfaces from the vertices the triangles share, and the `.log` says how many there
are. If some edges do not belong to exactly two triangles, the surface is not closed. The
`.log` then gives a warning, and BEST treats all the triangles as one surface, as before.

If a result is not a finite number (NaN or Inf), BEST stops with an error and exit code 1,
and writes nothing to the `.be` file.

To try the examples, copy them to a scratch directory first, since the programs write
their output next to their input:

- `examples/best/2eu7_02562`: `best -f 2eu7_02562 -v -mw 30000` takes a minute or less.
- `examples/best/1aki49_atom_01806`: two separate closed surfaces, made by sesmesh at probe 0.
  One is 49 atoms of lysozyme, the other one more atom 12 Å away. `make test` uses it.
- `examples/rcoal`: the `.his` files give the edge fraction of each output, so for
  example `rcoal -f test.c3p -e 0.8` reproduces `test_05770`.
- `examples/rcoal3`: `rcoal3 -f 1crn.face -nmax 2500 -nmin 1500 -n 3 -s`.

## Large meshes

Above 10,000 triangles BEST uses an H-matrix solver (`hmatrix/`) instead of the dense
matrix. It computes the same matrix entries, bit for bit. The blocks between
well-separated groups of triangles are stored in low-rank form, accurate to a relative
eps. It then solves by GMRES, preconditioned by a coarse H-LU factorization of the same
matrix. The results go through BEST's own code, so the `.be` and `.log` files are the same
kind as before; the `.log` also records the solver's memory, ranks, iterations and
residuals.

| Option | |
|---|---|
| `-hm` | the H-matrix solver at any number of triangles |
| `-dn` | the dense solver at any number of triangles |
| `-he eps` | the H-matrix accuracy (default `1e-4`) |
| `-hlu` | a direct H-LU factorization at eps instead of GMRES |
| `-hs` | the H-matrix and its factors in single precision: half the memory |

**Accuracy.** At the default eps of 1e-4 the H-matrix solver matches the dense solver
closely. These are the largest relative differences in D_t, D_r, CHI and [η]:

- **GMRES:**
  - 2.7e-6 on 14 meshes and flag sets of 2,000 to 15,000 triangles;
  - 6e-7 on a 30,050-triangle mesh.
- **`-hlu`:** 2.3e-5.

Both are far below the error of the 1/N extrapolation. `-he 1e-6` changes the GMRES results
by 6e-7 or less, on meshes of up to 254,292 triangles. At 30,050 triangles it is within
3e-9 of dense, and it takes about 1.5 times the memory and the time.

A larger eps saves little. With `-he 1e-3` the low-rank blocks shrink by a third, but the
peak memory only by 16 to 21%: the dense blocks between neighboring triangles and the
preconditioner stay the same. The results move by 2e-5 to 3e-5, about 50 times as far as at
1e-4, and all in the same direction, so the 1/N extrapolation does not average the error
out. To save memory, use `-hs` instead (see Single precision, below).

**Size.** These are lysozyme (1AKI) meshes from sesmesh, run with `-vc -mw 14305` and the
defaults (GMRES, eps 1e-4). The node is a 2 × EPYC 7713 with 64 threads. "Dense" is the
72 N² bytes that the dense matrix alone would take.

| Triangles | H-matrix | Preconditioner | GMRES iterations | Time | Peak memory | Dense |
|---|---|---|---|---|---|---|
| 30,050 | 4.0 GB | 0.7 GB | 7 | 10 s | 5 GB | 65 GB (315 s) |
| 60,278 | 8.7 GB | 1.6 GB | 8 | 23 s | 12 GB | 262 GB |
| 120,976 | 18.8 GB | 3.5 GB | 9 | 49 s | 26 GB | 1.1 TB |
| 254,292 | 42.5 GB | 7.8 GB | 9 | 1.7 min | 58 GB | 4.7 TB |
| 508,720 | 90.3 GB | 16.9 GB | 11 | 3.8 min | 124 GB | 19 TB |
| 1,019,462 | 190.7 GB | 36.1 GB | 11 | 7.9 min | 263 GB | 75 TB |

The memory per triangle grows slowly, from 155 KB at 120,000 triangles to 187 KB at
1 million, so a node with 256 GB stops just short of 1 million triangles at eps 1e-4.

**Single precision.** `-hs` stores the H-matrix and its factors in single precision. The
arithmetic stays in double precision: the products are summed in double, and the
compression and the factorization work on double copies. It takes half the memory and 8 to
11% less time, and it changes the results very little. These are the largest relative
differences in D_t, D_r, CHI and [η] from the same runs without `-hs`:

- **GMRES:** 1.1e-7 on the 14 meshes of 2,000 to 15,000 triangles; 6e-11 at 30,050 to
  1,019,462 triangles, and 2e-11 with `-he 1e-6` (30,050 to 254,292).
- **`-hlu`:** 2.2e-7 on the 14 small meshes; 3e-11 at 120,976 to 508,720 triangles.

The same meshes with `-hs`, and one of 2 million triangles:

| Triangles | H-matrix | Preconditioner | GMRES iterations | Time | Peak memory |
|---|---|---|---|---|---|
| 30,050 | 2.0 GB | 0.35 GB | 7 | 9 s | 3 GB |
| 60,278 | 4.3 GB | 0.8 GB | 8 | 20 s | 6 GB |
| 120,976 | 9.4 GB | 1.7 GB | 9 | 42 s | 13 GB |
| 254,292 | 21.2 GB | 3.9 GB | 9 | 1.5 min | 30 GB |
| 508,720 | 45.2 GB | 8.5 GB | 11 | 3.3 min | 63 GB |
| 1,019,462 | 95.4 GB | 18.0 GB | 11 | 7.3 min | 131 GB |
| 2,039,690 | 201.0 GB | 38.4 GB | 13 | 16 min | 274 GB |

So with `-hs` 1 million triangles fit on a node with 256 GB, and 2 million on this node
(512 GB). Without `-hs` the 2-million-triangle mesh would need about twice the memory,
more than the node has: above, `-hs` took 0.50 to 0.55 times the peak memory without it,
and 0.50 from 509k triangles up.

`-hlu` factors the whole matrix at eps instead, with no iterations. It needs more memory,
and it is slower except on small meshes. At 508,720 triangles it took 11 minutes and
170 GB (with `-hs`, 10 minutes and 81 GB), against GMRES's 3.8 minutes and 124 GB, and its
results were within 8e-6 of GMRES's. Its peak memory is about 1.7 times its factors.

The solver uses `OMP_NUM_THREADS` threads, or every hardware thread if that is not set.

**Memory checks.** The memory per triangle depends on the surface as well as on N. A rough
surface, such as the union of the atoms' spheres (probe 0, with its crevices and tunnels),
has about three times as many blocks per triangle as a smooth SES, and more than twice the
memory. With `-hs` and GMRES, a probe-0 mesh of lysozyme with 60,562 triangles peaked at
13.8 GB (228 KB per triangle), against 6.2 GB (103 KB) for its probe-3 SES of 60,004. A
probe-0 mesh of 1ADO with 667,230 triangles took 221 GB (332 KB), against 64 GB (145 KB)
for its probe-3 SES of 441,814.

- **Before the build.** BEST estimates the peak memory of the run and stops with an error if
  the node does not have that much available above the floor below. It counts the blocks
  exactly, from the block tree, and measures their ranks on a sample of 2% of them. The
  sample is kept for the matrix, so it costs no extra time. The `.log` gives the estimate and
  its parts. In 61 runs, on smooth and rough meshes of 26,000 to 667,000 triangles, in
  single and double precision and with 8 to 64 threads, it came out 1% to 12% above the
  peak with GMRES, and 8% to 14% above with `-hlu`. The estimate covers the first GMRES
  preconditioner. If BEST has to rebuild it finer (see Rough meshes, below), it checks the
  memory again first. `BEST_HMATRIX_MEMCHECK=0` turns both checks off.
- **During the run.** While the solver builds, factors and solves, a thread watches the
  memory available on the node. If it falls below a floor, the solver stops with an error and
  frees its memory, instead of running the node out of memory. On a shared node that would
  hit every job. The floor is 5% of the node's memory, and at least 1 GB.
  `BEST_HMATRIX_MEMFLOOR` sets it in GB, and 0 turns the watch off.

**Rough meshes.** Probe-0 surfaces (unions of atom spheres) and meshes with narrow tunnels
have two kinds of defect: sliver triangles folded onto their neighbors, and sheets facing each
other across gaps of hundredths of an Å. Two collocation points that close, on triangles with
opposite normals, make two rows of the matrix nearly the same. Each such pair gives a
near-null mode: a stress trapped in the fold or gap. That stress exerts no net force, so the
results hardly depend on it. On a 9,254-triangle probe-0 fragment of lysozyme, the dense
solver, GMRES and the H-matrix at eps 1e-6 agreed to every printed digit.

The coarse preconditioner breaks down on these meshes, though. Its smallest pivot ratio was
1e-29 to 1e-41, where smooth meshes give 1e-5 to 0.2, and GMRES then diverged or the
factorization hit a zero pivot. BEST now rebuilds the preconditioner at eps 1e-2 when that
happens, and at 1e-3 if needed. It also does so when GMRES does not converge, or breaks down
(its residual grows). A finer preconditioner takes more memory and time. On a
60,540-triangle probe-0 lysozyme mesh with `-hs`, it took 4.1 GB at 1e-2 and 7.5 GB at 1e-3,
against 1.5 GB at 1e-1 and 10.3 GB for the H-matrix. Both gave the same results, in 10 to 14
and 6 to 10 iterations. The `.log` then also lists the mesh's slivers and its close pairs of
triangles with opposite normals. Smooth meshes never take this path, so their results are
unchanged. A larger probe radius avoids the defects. `BEST_HMATRIX_PRE_EPS` sets the first eps
(default 1e-1).

## Versions

- **`v4.0-aragon`**: BEST 4.0, rcoal 2.0 and rcoal 3.0 as Sergio Aragon released them in
  2025. The Fortran sources are exactly as he sent them; only this README, the license and
  the Makefile were added.
- **`main`**: work since `v4.0-aragon`. Each change comes as a pull request. Releases will be
  tagged `v4.1`, `v4.2` and so on, and listed here.

Work since `v4.0-aragon`:

- rcoal: crash and memory fixes
  ([#1](https://github.com/ehb54/best-hydro/pull/1)); coarsened meshes kept closed and
  manifold, and triangulations that would give NaN in BEST refused with an error
  ([#2](https://github.com/ehb54/best-hydro/pull/2)).
- `make test` and continuous integration on GitHub
  ([#3](https://github.com/ehb54/best-hydro/pull/3)); its rcoal tests expect closed meshes
  where two of the author's examples are not closed
  ([#8](https://github.com/ehb54/best-hydro/pull/8)).
- BEST: 17 to 29 times faster with `-vc` on 2000 to 8000 triangles. It solves with the LU
  factors instead of inverting the matrix, and fills the matrix with OpenMP threads; the
  results are the same to rounding
  ([#9](https://github.com/ehb54/best-hydro/pull/9)).
- BEST: removes the near-null space of its matrix (a uniform pressure on the surface)
  before solving. Without the area correction (`-na`, or a mesh whose own area is exact)
  the torques, the centers, D_r and the viscosity were wrong; on a smooth 7948-triangle
  lysozyme mesh CHI came out as 1054 instead of 3.24. With the area correction, D_t, CHI
  and [η] change by 2.3e-6 or less
  ([#10](https://github.com/ehb54/best-hydro/pull/10)); `make test` checks a run with `-na`
  ([#12](https://github.com/ehb54/best-hydro/pull/12)).
- BEST: translation-rotation coupling. Dtr kept one of its three terms, and the
  transform to the center of resistance used the coupling tensor where its transpose
  belongs. The center of diffusion was off by 0.1 to 0.3 Å, and Drr by up to 2.4e-4. The
  results now agree with the direct inverse of the 6×6 friction tensor, and `make test`
  checks that a rotated mesh gives the same results
  ([#11](https://github.com/ehb54/best-hydro/pull/11)).
- BEST: with `-vc`, 11 solves instead of 19. The viscosity calculation reuses the
  diffusion calculation's solutions for the unit translations, and gets the strain flows
  about the center of viscosity from those about the centroid. The solves go through one
  routine, `BESOLV` in `best/besolve.f90`, which holds the matrix. The printed results are
  unchanged on 14 meshes and flag sets, and `make test` now checks `-vc`
  ([#13](https://github.com/ehb54/best-hydro/pull/13)).
- BEST: an H-matrix solver for large meshes, used above 10,000 triangles. It ran a
  lysozyme mesh of 1.02 million triangles on one node in 8 minutes and 263 GB. Its results
  agree with the dense solver's to about 1e-6 (2.7e-6 at most on 14 test meshes). With
  `-hs` it stores the matrix in single precision, in half the memory: a mesh of 2.04
  million triangles took 16 minutes and 274 GB
  ([#14](https://github.com/ehb54/best-hydro/pull/14)).
- BEST: sound solves on rough and multi-part meshes. A triangulation of several separate
  closed surfaces gets one near-singularity removed per surface. Before, D_t could come out
  negative, with exit code 0. When the H-matrix solver's coarse preconditioner breaks down
  (probe-0 and tunnel meshes), BEST rebuilds it finer, and logs the mesh's folds and gaps.
  GMRES stops as soon as it breaks down, and BEST no longer writes NaN to the `.be` file. On
  32 smooth lysozyme and pepsin meshes (probe 1.5 and 3, 29,000 to 448,000 triangles), the `.be`
  files are unchanged, byte for byte ([#16](https://github.com/ehb54/best-hydro/pull/16)).
- BEST: the H-matrix solver's memory estimate counts the matrix's blocks and samples their
  ranks, instead of scaling with the number of triangles alone. On rough surfaces the old
  estimate was half the peak: 95 GB for 221 GB on a 667,000-triangle probe-0 mesh. In 61 runs
  the new estimate came out 1% to 14% above the peak, where the old one ranged from 57% under
  to 38% over. A thread now also watches the node's available memory, and stops the solver
  with an error before the node runs out
  ([#17](https://github.com/ehb54/best-hydro/pull/17)).

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

E. H. Brookes (Department of Chemistry and Biochemistry, University of Montana) maintains
this repository for US-SOMO. The work since `v4.0-aragon` is supported by NIH grant
R01 GM120600.

## License

Released by the author under the MIT License. The changes since `v4.0-aragon` are under
the same license. See `LICENSE`.
