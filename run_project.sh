#!/bin/bash
#SBATCH --job-name=julia_serial
#SBATCH --account=def-amykim13-ab
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Load Julia
module load julia/1.10.10

# Run Julia (single-threaded, serial)
julia run_project.jl \
    --temporal \
    --original \
    --time_limit 28800 \
    --method 1 \
    --routes [4,6,11,14,15,26,27,28,30,31]
