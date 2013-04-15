#!/bin/bash -e

MODULE_CMD=/srv/gs1/apps/Modules/default/init/sh

if [ -z "$HUGESEQ_HOME" ]
then
	source $MODULE_CMD
	export MODULEPATH=`cd \`dirname $0\`/../../modulefiles; pwd`:$MODULEPATH
	echo $MODULEPATH
fi

module load hugeseq

$*
