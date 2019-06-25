#!/bin/bash

testdir="$(dirname $(realpath $BASH_SOURCE))"

tmpdir=""
datadir="$1"
if [ -z "$datadir" ] ; then
    tmpdir=$(mktemp -d)
    echo "using temp dir: $tmpdir"
    unzip -d $tmpdir $testdir/test_beedeposource.zip || exit 1
    datadir="$tmpdir/data"
    
fi
depofiles=$(for n in  data/?/?-truthDepo.json; do echo -n '"'$n'"',; done)


wire-cell -L debug -l stdout -C "depofiles=["$depofiles"]" -c $testdir/test_beedeposource.jsonnet || exit 1
if [ -n "$tmpdir" ] ; then
    rm -rf $tmpdir
fi
