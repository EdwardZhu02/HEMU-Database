#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o ./logs/maker_round1.out
#SBATCH -e ./logs/maker_round1.err 

source activate maker

export LIBDIR='/public1/home/sc80041/miniconda3/envs/maker/share/RepeatMasker/Matrices'
export LIBDIR='/public1/home/sc80041/miniconda3/envs/maker/share/RepeatMasker/Libraries'

mpiexec -n 40 maker \
	-base Themeda_rnd1 \
	maker_opts-firstrun-themeda-parallel.ctl maker_bopts.ctl maker_exe.ctl 