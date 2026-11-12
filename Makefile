# Builds BEST, rcoal (2.0) and rcoal3 (3.0) with gfortran, and the surface mesher sesmesh
# with a C++17 compiler, into bin/.
#
#   make                                 Linux with reference LAPACK/BLAS
#   make LAPACK="-lopenblas"             another LAPACK/BLAS
#   make LAPACK="-framework Accelerate"  macOS
#   make OPENMP=                         sesmesh without OpenMP (one thread)
#   make test-sesmesh                    sesmesh's unit tests
#
# The author built these with PGI and -r8, which makes every REAL 8 bytes; the
# gfortran equivalent is -fdefault-real-8 -fdefault-double-8.  He recommends
# optimization no higher than -O2.  Some source lines run past column 72, which
# gfortran rejects by default, hence -ffixed-line-length-none.

FC     = gfortran
FFLAGS = -O2
R8     = -fdefault-real-8 -fdefault-double-8
FIXED  = -ffixed-line-length-none
LAPACK = -llapack -lblas

CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra
OPENMP   = -fopenmp
OMPFLAGS = $(if $(OPENMP),$(OPENMP),-Wno-unknown-pragmas)

BEST_SRC   = best/best.f best/input.f best/tgetgr.f best/triangleint.f \
             best/jacobi.f best/matin3.f best/visco.f best/difco.f best/output.f
RCOAL_SRC  = rcoal/rcoal.f rcoal/scoalesce.f rcoal/index.f
RCOAL3_SRC = rcoal3/rcoal3.f rcoal3/scoalesce3.f rcoal3/index.f
MESH_SRC   = mesher/ses.cpp mesher/mesh.cpp
MESH_HDR   = mesher/vec3.h mesher/ses.h mesher/mesh.h

all: bin/best bin/rcoal bin/rcoal3 bin/sesmesh

bin:
	mkdir -p bin

bin/best: $(BEST_SRC) | bin
	$(FC) $(FFLAGS) $(R8) $(FIXED) -o $@ $(BEST_SRC) $(LAPACK)

bin/rcoal: $(RCOAL_SRC) | bin
	$(FC) $(FFLAGS) $(R8) $(FIXED) -o $@ $(RCOAL_SRC)

bin/rcoal3: $(RCOAL3_SRC) | bin
	$(FC) $(FFLAGS) $(R8) $(FIXED) -o $@ $(RCOAL3_SRC)

bin/sesmesh: mesher/sesmesh.cpp $(MESH_SRC) $(MESH_HDR) | bin
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ mesher/sesmesh.cpp $(MESH_SRC)

bin/test_sesmesh: mesher/test_sesmesh.cpp $(MESH_SRC) $(MESH_HDR) | bin
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ mesher/test_sesmesh.cpp $(MESH_SRC)

test-sesmesh: bin/test_sesmesh
	bin/test_sesmesh

clean:
	rm -rf bin

.PHONY: all clean test-sesmesh
