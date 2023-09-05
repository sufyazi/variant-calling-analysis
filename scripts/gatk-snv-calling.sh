#!/usr/bin/env bash

module load singularity
singExec="singularity exec --bind /home/users/ntu/suffiazi/scratch:/mnt /home/users/ntu/suffiazi/apps/sifs/gatk.sif"

$singExec gatk --help

$singExec gatk ValidateSamFile  -I /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/test_files/NA12878.chr17_69k_70k.dictFix.bam -MODE SUMMARY