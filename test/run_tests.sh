#!/bin/bash
# Runs BEST, rcoal and rcoal3 on the examples and checks the results.
#
#   test/run_tests.sh [bin-directory]     (default: bin; "make test" runs this)
#
# Work files go in a new directory under $TEST_WORKDIR, or the system's
# temporary directory if that is not set, and are kept if a test fails.
# Needs Perl.

set -u
root=$(cd "$(dirname "$0")/.." && pwd)
bin=$(cd "${1:-$root/bin}" && pwd) || exit 1
t=$root/test
work=$(mktemp -d "${TEST_WORKDIR:-${TMPDIR:-/tmp}}/best-hydro-test.XXXXXX") || exit 1
fails=0
ok()  { echo "ok   - $*"; }
bad() { echo "FAIL - $*"; fails=$((fails+1)); }

# BEST: the author's sample input
d=$work/best; mkdir -p "$d"; cd "$d" || exit 1
cp "$root/examples/best/2eu7_02562" .
if "$bin/best" -f 2eu7_02562 -v -mw 30000 > best.out 2>&1 \
   && perl "$t/check_best.pl" 2eu7_02562vm.be "$t/expected_best_2eu7.txt"; then
   ok "best: 2eu7_02562 -v -mw 30000"
else
   bad "best: 2eu7_02562 -v -mw 30000 (see $d)"
fi

# BEST: the viscosity at the center of viscosity (-vc)
if "$bin/best" -f 2eu7_02562 -vc -mw 30000 > best_vc.out 2>&1 \
   && perl "$t/check_best.pl" 2eu7_02562vcm.be "$t/expected_best_2eu7_vc.txt"; then
   ok "best: 2eu7_02562 -vc -mw 30000"
else
   bad "best: 2eu7_02562 -vc -mw 30000 (see $d)"
fi

# BEST: the same mesh, rotated, must give the same results
perl "$t/permute_mesh.pl" 2eu7_02562 2eu7_rot
if "$bin/best" -f 2eu7_rot -v -mw 30000 > best_rot.out 2>&1 \
   && perl "$t/check_best_frame.pl" 2eu7_02562vm.be 2eu7_02562vm.log 2eu7_rotvm.be; then
   ok "best: 2eu7_02562 rotated gives the same results"
else
   bad "best: 2eu7_02562 rotated should give the same results (see $d)"
fi

# BEST: the same mesh without the area correction.  The correction keeps GR
# away from singular; without it, BEST must deflate GR's near-null space (a
# uniform pressure on the surface) or the results are wrong.
if "$bin/best" -f 2eu7_02562 -na -v -mw 30000 > best_na.out 2>&1 \
   && perl "$t/check_best.pl" 2eu7_02562vmn.be "$t/expected_best_2eu7_na.txt"; then
   ok "best: 2eu7_02562 -na -v -mw 30000"
else
   bad "best: 2eu7_02562 -na -v -mw 30000 (see $d)"
fi

# BEST with the H-matrix solver (-hm), forced on this small mesh.  With
# BEST_HMATRIX_CHECK set, BEST also fills GR densely and logs whether the
# solver's entries are TGETGR's bit for bit.  The results must agree with the
# dense ones above to 2e-5 (or the last digit printed): GMRES with -v, and the
# direct H-LU with -vc.
d=$work/best-hm; mkdir -p "$d"; cd "$d" || exit 1
cp "$root/examples/best/2eu7_02562" .
BEST_HMATRIX_CHECK=1 "$bin/best" -f 2eu7_02562 -v -mw 30000 -hm > best_hm.out 2>&1
if grep -q "built without the H-matrix solver" 2eu7_02562vm.log; then
   echo "skip - best -hm: BEST was built without the H-matrix solver (HMATRIX=)"
else
   if grep -q "GR check: 0 of [0-9]* entries and 0 of [0-9]* incenter" 2eu7_02562vm.log; then
      ok "best -hm: the H-matrix solver's GR entries are TGETGR's, bit for bit"
   else
      bad "best -hm: the H-matrix solver's GR entries should be TGETGR's (see $d)"
   fi
   if perl "$t/check_best_pair.pl" "$work/best/2eu7_02562vm.be" 2eu7_02562vm.be 2e-5; then
      ok "best -hm: 2eu7_02562 -v -mw 30000 agrees with the dense solver"
   else
      bad "best -hm: 2eu7_02562 -v -mw 30000 should agree with the dense solver (see $d)"
   fi
   if "$bin/best" -f 2eu7_02562 -vc -mw 30000 -hm -hlu > best_hlu.out 2>&1 \
      && perl "$t/check_best_pair.pl" "$work/best/2eu7_02562vcm.be" 2eu7_02562vcm.be 2e-5; then
      ok "best -hm -hlu: 2eu7_02562 -vc -mw 30000 agrees with the dense solver"
   else
      bad "best -hm -hlu: 2eu7_02562 -vc -mw 30000 should agree with the dense solver (see $d)"
   fi
fi

# rcoal: the author's test case.  Each mesh is regenerated from the edge
# fraction recorded in its .his file, and must be closed.  It must also
# reproduce the author's mesh, unless test/expected_rcoal.txt gives another
# triangle count for it: there the author's mesh is not closed, and rcoal now
# skips the collapses that made it so.
d=$work/rcoal; mkdir -p "$d"; cd "$d" || exit 1
cp "$root/examples/rcoal/test.c3p" "$root/examples/rcoal/test.c3v" .
for his in "$root"/examples/rcoal/test_*.his; do
   ref=${his%.his}
   name=$(basename "$ref")
   e=$(awk '/Edge fraction/ { print $NF; exit }' "$his")
   n=$(awk -v m="$name" '$1 == m { print $2; exit }' "$t/expected_rcoal.txt")
   if [ -z "$n" ]; then
      if "$bin/rcoal" -f test.c3p -e "$e" > "$name.out" 2>&1 && [ -f "$name" ] \
         && perl "$t/meshcheck.pl" "$name" && perl "$t/meshcmp.pl" "$ref" "$name" 1e-4; then
         ok "rcoal: -e $e reproduces $name"
      else
         bad "rcoal: -e $e should reproduce $name (see $d)"
      fi
   else
      out=$(printf '%s_%05d' "${name%_*}" "$n")
      if "$bin/rcoal" -f test.c3p -e "$e" > "$name.out" 2>&1 && [ -f "$out" ] \
         && perl "$t/meshcheck.pl" "$out"; then
         ok "rcoal: -e $e gives a closed $out (the author's $name is not closed)"
      else
         bad "rcoal: -e $e should give a closed $out (see $d)"
      fi
   fi
done

# rcoal3: crambin, triangulated by MSMS
d=$work/rcoal3; mkdir -p "$d"; cd "$d" || exit 1
cp "$root"/examples/rcoal3/1crn.* .
if "$bin/rcoal3" -f 1crn.face -nmax 2500 -nmin 1500 -n 3 -s > rcoal3.out 2>&1; then
   meshes=$(ls 1crn_[0-9][0-9][0-9][0-9][0-9] 2>/dev/null)
   n=$(echo $meshes | wc -w)
   if [ "$n" -eq 3 ] && perl "$t/meshcheck.pl" $meshes; then
      ok "rcoal3: 1crn.face -s gives 3 closed meshes"
   else
      bad "rcoal3: 1crn.face -s should give 3 closed meshes, got $n (see $d)"
   fi
else
   bad "rcoal3: 1crn.face -s failed (see $d)"
fi

if [ $fails -eq 0 ]; then
   rm -rf "$work"
   echo "all tests passed"
else
   echo "$fails test(s) failed; files kept in $work"
   exit 1
fi
