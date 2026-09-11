#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --job-name=ohc_map_emitter
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=katie.mills@colorado.edu
#SBATCH --mem=32GB
#SBATCH --exclude=blanca-g6-u21

module purge
module load anaconda/2023.09
conda activate ohc_crosscheck

cd $SLURM_SUBMIT_DIR
tag=260908-OP20260507
results=/scratch/alpine/wimi7695/ohc_prod/results_260908
window=$1
level=$2
python emit.py \
  ${results}/derive_${tag}_${window}_${level}.nc \
  --tag ${tag} \
  --provenance-link https://github.com/ocean-grid-processing/provenance/blob/main/README.md \
  --code-version https://github.com/ocean-grid-processing/ohc_map_emitter/releases/tag/1.0.0 \
  --out ${results}
