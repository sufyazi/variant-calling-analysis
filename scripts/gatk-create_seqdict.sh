#!/usr/bin/env bash
#PBS -l select=1:ncpus=12:mem=120GB
#PBS -l walltime=3:00:00
#PBS -j oe
#PBS -P personal
#PBS -q normal

module load singularity
singExec="singularity exec --bind /home/users/ntu/suffiazi/scratch:/mnt /home/users/ntu/suffiazi/apps/sifs/gatk.sif"

$singExec gatk --help

$singExec gatk CreateSequenceDictionary --REFERENCE /mnt/inputs/references/gatk4/hg38_GENCODE_primary_assembly_genome.fa