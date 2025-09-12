#!/bin/bash 
#$ -N multifates_ot
#$ -l m_mem_free=30G
#$ -o job_output
#$ -j y

source miniconda3/bin/activate

samtools view -b /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson/snATAC/version_2.1/Control_6/outs/possorted_bam.bam chrM > /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/Control_6/mtDNA/chrM.bam
cd /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/Control_6/mtDNA
samtools sort -o sorted_chrM.bam chrM.bam
samtools index sorted_chrM.bam
sinto fragments -b sorted_chrM.bam -p 2 -f fragments_collapse_within.bed --collapse_within



