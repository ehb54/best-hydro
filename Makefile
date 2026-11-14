# Builds BEST, rcoal (2.0) and rcoal3 (3.0) with gfortran into bin/.
#
#   make                                 Linux with reference LAPACK/BLAS
#   make LAPACK="-lopenblas"             another LAPACK/BLAS
#   make LAPACK="-framework Accelerate"  macOS
#   make OPENMP=                         BEST without OpenMP
#
# The author built these with PGI and -r8, which makes every REAL 8 bytes; the
# gfortran equivalent is -fdefault-real-8 -fdefault-double-8.  He recommends
# optimization no higher than -O2.  Some source lines run past column 72, which
# gfortran rejects by default, hence -ffixed-line-length-none.
#
# BEST fills its matrix with OpenMP threads (OMP_NUM_THREADS); the LAPACK
# calls use the threads of the LAPACK/BLAS library (e.g. OPENBLAS_NUM_THREADS).
#
# best/besolve.f90 is a module (GR and the solves with it), so it comes first;
# its .mod file goes in build/.
#
# Above 10000 triangles BEST uses the H-matrix solver in hmatrix/ (C++17,
# std::thread; threads: OMP_NUM_THREADS, or every hardware thread).  It needs a
# C++ compiler; its objects go in build/ and BEST is linked with CXXLIBS.
#   make HMATRIX=                        BEST without it (dense solver only)
#   make CXXLIBS=-lc++                   macOS with clang++ as CXX
# After changing HMATRIX, run make clean first.

FC     = gfortran
FFLAGS = -O2
R8     = -fdefault-real-8 -fdefault-double-8
FIXED  = -ffixed-line-length-none
OPENMP = -fopenmp
LAPACK = -llapack -lblas
CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pthread -ffp-contract=off
CXXLIBS  = -lstdc++
HMATRIX  = 1

BEST_SRC   = best/besolve.f90 best/best.f best/input.f best/tgetgr.f best/triangleint.f \
             best/jacobi.f best/matin3.f best/visco.f best/difco.f best/output.f
# The H-matrix solver.  -ffp-contract=off keeps the compiler from fusing
# multiplies and adds, so that its matrix entries match TGETGR's bit for bit.
HM_SRC = hmatrix/src/api.cpp hmatrix/src/kernel.cpp hmatrix/src/cluster.cpp \
         hmatrix/src/hmatrix.cpp hmatrix/src/mat.cpp hmatrix/src/gmres.cpp \
         hmatrix/src/parallel.cpp
HM_HDR = hmatrix/best_hmatrix.h $(wildcard hmatrix/src/*.h)
HM_OBJ = $(HM_SRC:hmatrix/src/%.cpp=build/%.o)

ifeq ($(HMATRIX),)
BEST_HM      = best/hmstub.f90
BEST_HM_LIBS =
else
BEST_HM      = $(HM_OBJ)
BEST_HM_LIBS = $(CXXLIBS) -pthread
endif

RCOAL_SRC  = rcoal/rcoal.f rcoal/scoalesce.f rcoal/index.f
RCOAL3_SRC = rcoal3/rcoal3.f rcoal3/scoalesce3.f rcoal3/index.f

all: bin/best bin/rcoal bin/rcoal3

bin:
	mkdir -p bin

build:
	mkdir -p build

bin/best: $(BEST_SRC) $(BEST_HM) | bin build
	$(FC) $(FFLAGS) $(OPENMP) $(R8) $(FIXED) -Jbuild -o $@ $(BEST_SRC) $(BEST_HM) $(LAPACK) $(BEST_HM_LIBS)

build/%.o: hmatrix/src/%.cpp $(HM_HDR) | build
	$(CXX) $(CXXFLAGS) -Ihmatrix -c -o $@ $<

bin/rcoal: $(RCOAL_SRC) | bin
	$(FC) $(FFLAGS) $(R8) $(FIXED) -o $@ $(RCOAL_SRC)

bin/rcoal3: $(RCOAL3_SRC) | bin
	$(FC) $(FFLAGS) $(R8) $(FIXED) -o $@ $(RCOAL3_SRC)

# Runs the examples and checks the results (needs Perl).
test: all
	test/run_tests.sh bin

clean:
	rm -rf bin build

.PHONY: all test clean
