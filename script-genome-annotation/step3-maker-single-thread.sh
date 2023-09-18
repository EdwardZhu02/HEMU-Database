#!/bin/bash

# Environment req: maker

export LIBDIR=/home/test3/miniconda3/envs/maker/share/RepeatMasker/Matrices
export LIBDIR=/home/test3/miniconda3/envs/maker/share/RepeatMasker/Libraries

# Check log folder
if [ ! -d "./logs" ]; then
	mkdir ./logs
fi

maker \
	-base Microstegium_rnd1 \
	-cpus 24 \
	maker_opts-firstrun-microstegium.ctl maker_bopts.ctl maker_exe.ctl