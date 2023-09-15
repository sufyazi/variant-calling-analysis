#!/usr/bin/env bash
# shellcheck disable=SC1091

#PBS -l select=1:ncpus=12:mem=120GB
#PBS -l walltime=3:00:00
#PBS -j oe
#PBS -P 12003580
#PBS -q normal

# load conda environment
module load singularity
singExec="singularity exec --bind /home/users/ntu/suffiazi/scratch:/mnt /home/users/ntu/suffiazi/apps/sifs/gatk.sif"


# set up variables
bam=$BAMFILE
file_prefix=$PREFIX

# run gatk MarkIlluminaAdapters
$singExec gatk --java-options "-Xmx80G" MarkIlluminaAdapters \
        -I "${bam}" \
        -O "${file_prefix}"_markadapt.bam \
        -M "${file_prefix}"_markadapt-metrcs.txt \
        -TMP_DIR /mnt/outputs/tmp