for level in 15_20 15_300 300_700 700_1000 700_1850 1800_1850
do
    sbatch combine.slurm ${level}
done
