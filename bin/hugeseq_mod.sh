#!/bin/bash -e

module load hugeseq/1.2
export TMP=$1
echo $TMP
shift

$*
