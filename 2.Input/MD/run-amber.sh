#!/bin/bash
#SBATCH --job-name=md
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --mem=30GB
#SBATCH -M pudong
#SBATCH --gres=gpu:1
#SBATCH -p sci
#SBATCH -w gpu11

pmemd.cuda -O -i min1.in -o min1.out -p c_sol.prmtop -c c_sol.inpcrd -r min1.rst -ref c_sol.inpcrd
pmemd.cuda -O -i min2.in -o min2.out -p c_sol.prmtop -c min1.rst -r min2.rst 
pmemd.cuda -O -i heat.in -o md0.out -p c_sol.prmtop -c min2.rst -r md0.rst -ref min2.rst
pmemd.cuda -O -i md1.in -o md1.out -p c_sol.prmtop -c md0.rst -r md1.rst 
pmemd.cuda -O -i md2.in -o md2.out -p c_sol.prmtop -c md1.rst -r md2.rst -x md2.mdcrd
 


