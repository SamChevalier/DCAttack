#!/bin/bash
#SBATCH --job-name=solve_118          # Job name
#SBATCH --output=out_%j.log          # Stdout and stderr log (%j = job ID)
#SBATCH --partition=standard         # Partition (queue) to use
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=4            # Number of CPUs per task
#SBATCH --mem=64G                    # Memory per node (or use --mem-per-cpu)
#SBATCH --time=12:00:00              # Walltime (hh:mm:ss)
#SBATCH --mail-type=END,FAIL         # Notifications
#SBATCH --mail-user=your_netid@uvm.edu  # Your UVM email

# Load required modules
module load julia
module load gurobi


# Your commands here
cd DCAttack/
git checkout schev
julia run_vacc_tests.jl
