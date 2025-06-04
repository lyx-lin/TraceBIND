#!/bin/bash 
#$ -N hint
#$ -l m_mem_free=5G
#$ -o job_output
#$ -j y

source miniconda3/bin/activate

rgt-hint footprinting --atac-seq --organism=hg38 --output-location=/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/degron/CTCF/treatment/HINT /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/degron/CTCF/treatment/outs/possorted_bam.bam /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/degron/CTCF/peaks_region.bed
 