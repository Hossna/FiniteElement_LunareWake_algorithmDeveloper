#!/bin/bash
#SBATCH --account=def-fdam
#SBATCH --time=30:00:00
#SBATCH --job-name=test
#SBATCH --output=%x-%j.out
#SBATCH --mem-per-cpu=2000M

echo "Current working directory is `pwd`"
echo "---------------------"
/home/hossna/project/hossna/FiniteElementCode/general_ObliqueB/topo



