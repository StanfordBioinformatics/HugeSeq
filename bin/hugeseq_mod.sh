#!/bin/bash -e

module load hugeseq/2.0
export TMP=$1
echo $TMP
shift
export LOGFILE=$1
echo $LOGFILE
shift

$*
