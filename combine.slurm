#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=ohc_gcos_emitter
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=katie.mills@colorado.edu
#SBATCH --mem=128GB
#SBATCH --exclude=blanca-g6-u21

module purge
module load anaconda/2023.09
conda activate ohc_crosscheck

cd $SLURM_SUBMIT_DIR
python combine.py /scratch/alpine/wimi7695/ohc_prod/results/OHC_*.nc --out /scratch/alpine/wimi7695/ohc_prod/results/ \
  --tag $1 --provenance-link https://github.com/argovis/ohc_map_emitter/blob/main/provenance/README.md \
  --levels 0_300,0_700,0_1000,700_2000,0_2000
