#!/usr/bin/env bash

module load singularity
singExec="singularity exec --bind /home/users/ntu/suffiazi/scratch:/mnt /home/users/ntu/suffiazi/apps/sifs/gatk.sif"

$singExec gatk --help

$singExec gatk CreateSequenceDictionary --REFERENCE /mnt/inputs/references/gatk4/hg38_GENCODE_primary_assembly_genome.fa