for window in 2005_2024 2005_2025
do
  for level in 0_300 0_700 0_1000 0_2000 700_2000
  do
    sbatch emit.slurm $window $level
  done
done
