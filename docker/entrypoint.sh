#!/bin/sh

. /opt/conda/etc/profile.d/conda.sh
conda activate base
conda activate mc

if [ "$@" == "none" ]; then
	bash
else
	$@
fi
