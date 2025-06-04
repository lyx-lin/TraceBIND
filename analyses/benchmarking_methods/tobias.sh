#!/bin/bash 
#$ -N Tobias
#$ -l m_mem_free=5G
#$ -o job_output
#$ -j y

source miniconda3/bin/activate

# TOBIAS ATACorrect --bam /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson/snATAC/version_2.1/Control_3/outs/possorted_bam.bam --genome /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/rgtdata/hg38/genome_hg38.fa --peaks //home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson/snATAC/version_2.1/Control_3/outs/peaks.bed --blacklist /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/shared/BlacklistFiles/hg38-blacklist.v2.bed --outdir /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/global_ref/Control_3/tobias --cores 2

TOBIAS FootprintScores --signal /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/global_ref/Control_3/tobias/possorted_bam_corrected.bw --regions /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/global_ref/Control_3/peaks_region.bed --output /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/global_ref/Control_3/tobias/TOBIAS_footprints.bw --cores 2

TOBIAS BINDetect --motifs /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt --signals /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/global_ref/Control_3/tobias/TOBIAS_footprints.bw --genome /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/rgtdata/hg38/genome_hg38.fa --peaks /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/HIFIA_unibind/Control_3/peaks_region.bed --outdir /home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/global_ref/Control_3/tobias --cores 2
