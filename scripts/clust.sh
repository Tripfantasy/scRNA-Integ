#!/usr/bin/bash

#SBATCH --partition gpu # Tell SLURM to look for nodes in the 'GPU' partition
#SBATCH --nodes=1 # Only use one node
#SBATCH --ntasks=1 # Use one CPU
#SBATCH --gres=gpu:1 # Use one GPU
#SBATCH --job-name "GraphST"

cd /home/dm2763/cell-type
micromamba activate GraphST
/usr/bin/time -v python clust.py \
-r /home/dm2763/micromamba/envs/GraphST/lib/R -d /home/dm2763/cell-type/data1.h5ad -m cuda -t Stereo -n 22