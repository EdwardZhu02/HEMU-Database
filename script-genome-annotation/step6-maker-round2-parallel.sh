#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o ./logs/maker_round2.out
#SBATCH -e ./logs/maker_round2.log

source activate maker

export LIBDIR='/public1/home/sc80041/miniconda3/envs/maker/share/RepeatMasker/Matrices'
export LIBDIR='/public1/home/sc80041/miniconda3/envs/maker/share/RepeatMasker/Libraries'

mpiexec -n 40 maker \
	-base Chrysopogon_rnd2 \
	maker_opts-secondrun-chrysopogon-parallel.ctl maker_bopts.ctl maker_exe.ctl 