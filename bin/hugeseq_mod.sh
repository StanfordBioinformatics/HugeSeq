#!/bin/bash -e

MODULE_CMD=~/app/Modules/default/init/sh

if [ -z "$HUGESEQ_HOME" ]
then
	source $MODULE_CMD
	export MODULEPATH=`cd \`dirname $0\`/../modulefiles; pwd`:$MODULEPATH
fi

module load hugeseq

$*
