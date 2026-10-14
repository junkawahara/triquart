#!/bin/sh
# tests/test_cli.sh — CLI regression for the parameterized-distance output.
#
# ((a,b),(c,d)) vs (a,b,c,d): the single quartet is resolved in T1 only, so
# A=B=D=E=0, C=1 and B + p (C + D) = p.  Trial 8 (Finding 2) found that the
# CLI printed this with one decimal, i.e. p=0.01 as 0.0 and p=0.25 as 0.2.
# The expected strings below are what %.{LDBL_DIG}Lg prints for exact p.
#
#   usage: sh tests/test_cli.sh build/soda13
set -u
bin=${1:-build/soda13}
dir=$(dirname "$0")/cli
fail=0
check() {  # check <p> <expected quartet line>
    out=$("$bin" --quartet -p "$1" "$dir/resolved4.nwk" "$dir/star4.nwk") || { echo "cli: exit status $? for p=$1"; fail=1; return; }
    if [ "$out" != "$2" ]; then
        echo "cli: p=$1"
        echo "  expected: $2"
        echo "  got:      $out"
        fail=1
    fi
}
check 0.01 "quartets: A=0 B=0 C=1 D=0 E=0 distance=1 parameterized(p=0.01)=0.01"
check 0.25 "quartets: A=0 B=0 C=1 D=0 E=0 distance=1 parameterized(p=0.25)=0.25"
check 1    "quartets: A=0 B=0 C=1 D=0 E=0 distance=1 parameterized(p=1)=1"
check 0    "quartets: A=0 B=0 C=1 D=0 E=0 distance=1 parameterized(p=0)=0"
check 0.1  "quartets: A=0 B=0 C=1 D=0 E=0 distance=1 parameterized(p=0.1)=0.1"
# Triplets of the same pair: 4 triplets, each resolved in T1 only, so the
# result is 4p; p = 0.1 and p = 0.01 exercise the product rounding.
checkt() {  # checkt <p> <expected triplet line>
    out=$("$bin" --triplet -p "$1" "$dir/resolved4.nwk" "$dir/star4.nwk") || { echo "cli: exit status $? for p=$1"; fail=1; return; }
    if [ "$out" != "$2" ]; then echo "cli: triplet p=$1"; echo "  expected: $2"; echo "  got:      $out"; fail=1; fi
}
checkt 0.5  "triplets: A=0 B=0 C=4 D=0 E=0 distance=4 parameterized(p=0.5)=2"
checkt 0.1  "triplets: A=0 B=0 C=4 D=0 E=0 distance=4 parameterized(p=0.1)=0.4"
checkt 0.01 "triplets: A=0 B=0 C=4 D=0 E=0 distance=4 parameterized(p=0.01)=0.04"
if [ $fail -eq 0 ]; then echo "test_cli: ok"; fi
exit $fail
