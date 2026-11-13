#!/bin/sh
# Compares bin/rcoal and bin/rcoal3 with the Fortran rcoal 2.0 and 3.0 of ehb54/best-hydro, as
# fixed there in PRs #1 and #2 (commit 347191b): every output file must be byte for byte the
# same, and so must stdout.  Run from the top of the tree after make (make test does both).
#
# Inputs, in test/rcoal:
#  - 1crn.face, 1crn.vert, 1crn.msms.log: crambin from MSMS, best-hydro's rcoal 3.0 example;
#  - crn.c3p.gz, crn.c3v: crambin from sesmesh (main e5f329b),
#    sesmesh -p 1.5 -n 6000 -f c3p, from msroll's radii; 5962 triangles.
# For each case, NAME.sha256 holds the sha256 of each file the Fortran wrote and NAME.stdout
# its stdout; it wrote nothing to stderr and exited with 0.  The Fortran was built with
# best-hydro's Makefile (gfortran 13, Linux x86-64).
#
# Then -search, which the Fortran does not have, must reach each target within 1%.
set -u
top=$(cd "$(dirname "$0")/.." && pwd)
data=$top/test/rcoal
work=$top/bin/test_rcoal
fail=0

sha() {
   if command -v sha256sum > /dev/null 2>&1; then
      sha256sum "$@"
   else
      shasum -a 256 "$@"
   fi
}

# setup NAME: a fresh directory with the inputs
setup() {
   d=$work/$1
   rm -rf "$d"
   mkdir -p "$d"
   cd "$d" || exit 1
   cp "$data/1crn.face" "$data/1crn.vert" "$data/1crn.msms.log" "$data/crn.c3v" .
   gzip -dc "$data/crn.c3p.gz" > crn.c3p
}

# compare NAME PROGRAM ARGS...: run and compare with the Fortran's results
compare() {
   name=$1
   prog=$2
   shift 2
   setup "$name"
   "$top/bin/$prog" "$@" > stdout 2> stderr
   rc=$?
   rm -f 1crn.face 1crn.vert 1crn.msms.log crn.c3p crn.c3v
   ok=1
   if [ $rc -ne 0 ]; then
      echo "FAIL $name: exit code $rc"
      ok=0
   fi
   if [ -s stderr ]; then
      echo "FAIL $name: stderr is not empty:"
      cat stderr
      ok=0
   fi
   if ! cmp -s stdout "$data/$name.stdout"; then
      echo "FAIL $name: stdout differs from the Fortran's"
      diff "$data/$name.stdout" stdout | head -20
      ok=0
   fi
   ls | grep -v -e '^stdout$' -e '^stderr$' -e '^sha256$' | LC_ALL=C sort | while read -r f; do sha "$f"; done > sha256
   if ! cmp -s sha256 "$data/$name.sha256"; then
      echo "FAIL $name: output files differ from the Fortran's (expected, then these):"
      diff "$data/$name.sha256" sha256
      ok=0
   fi
   if [ $ok -eq 1 ]; then
      echo "rcoal $name: $(wc -l < sha256 | tr -d ' ') files identical to the Fortran's"
   else
      fail=1
   fi
   cd "$top" || exit 1
}

compare crn_n rcoal -f crn.c3p -nmax 2500 -nmin 1500 -n 3
compare crn_e rcoal -f crn.c3p -e 0.8
compare crn_v3 rcoal3 -f crn.c3p -nmax 2000 -nmin 1200 -n 2
compare msms rcoal3 -f 1crn.face -nmax 2500 -nmin 1500 -n 3 -s

# -search
setup search
"$top/bin/rcoal" -f crn.c3p -nmax 2500 -nmin 1500 -n 3 -search > stdout 2> stderr
rc=$?
res=$(awk '/Target n/ { t = $3 } /Actual coalesce n/ { a = $4; d = a - t; if ( d < 0 ) d = -d;
           printf "%s %d->%d", ( 100 * d <= t ? "ok" : "BAD" ), t, a; printf " " }' stdout)
nfiles=$(ls crn_?????.his 2> /dev/null | wc -l | tr -d ' ')
if [ $rc -ne 0 ] || [ -s stderr ] || [ "$nfiles" -ne 3 ] || echo "$res" | grep -q BAD || [ -z "$res" ]; then
   echo "FAIL search: exit code $rc, $nfiles meshes, targets: $res"
   fail=1
else
   echo "rcoal -search: $res"
fi
cd "$top" || exit 1

if [ $fail -ne 0 ]; then
   echo "rcoal comparison FAILED"
   exit 1
fi
echo "all rcoal tests passed"
