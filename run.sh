for level in 15_20 15_300 300_700 700_1000 700_1850 1800_1850
do
    dir=/scratch/alpine/wimi7695/ohc_prod/potential_temperature/OP20260507
    sbatch combine.slurm ${dir}
done
