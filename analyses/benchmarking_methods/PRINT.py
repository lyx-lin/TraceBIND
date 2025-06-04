import pandas as pd
import numpy as np
import scprinter as scp
import os 
import gzip 

os.chdir('/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA/A2/mtDNA')

metadata = pd.read_csv("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson/metadata.csv")
metadata = metadata[metadata['library_id'] == 'A2']
barcodeGroups = metadata[['barcode', 'celltype']]

# Rename the columns
barcodeGroups.columns = ['barcodeID', 'group']

# Set the 'group' column to 1
barcodeGroups['group'] = 1
barcodes = barcodeGroups.iloc[:, 0]

# get regions
regions = pd.read_table("home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA/A2/mtDNA/peaks_region.bed", sep='\t')
regions.columns = ['Chromosome', 'Start', 'End']


printer = scp.pp.import_fragments(path_to_frags="/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA/A2/mtDNA/fragments_collapse_within.bed",
                    barcodes=barcodes,
                    savename='PRINT.h5ad',
                    genome=scp.genome.hg38,
                    min_num_fragments=0, min_tsse=0,
                    sorted_by_barcode=False)

printer.load_disp_model()
printer.load_bindingscore_model("TF",
                             scp.datasets.pretrained_TFBS_model)
printer.load_bindingscore_model("Nuc",
                             scp.datasets.pretrained_NucBS_model)

grouping, uniq_groups = scp.utils.df2cell_grouping(printer, barcodeGroups)

regions.iloc[:, 2] = np.max(regions.iloc[:, 2] - regions.iloc[:, 1]) + regions.iloc[:, 1]

scp.tl.get_binding_score(
    printer,
    grouping,
    uniq_groups,
    regions,
    model_key='TF',
    n_jobs=2, # nCores to use
    contextRadius=100,
    save_key="TF",
    backed=False,
    overwrite=True)

scp.tl.get_binding_score(
    printer,
    grouping,
    uniq_groups,
    regions,
    model_key='Nuc',
    n_jobs=2, # nCores to use
    contextRadius=100,
    save_key="Nuc",
    backed=False,
    overwrite=True)

tf = printer.bindingscoreadata['TF']
tf.write_h5ad("PRINT_TF.h5ad")

nuc = printer.bindingscoreadata['Nuc']
nuc.write_h5ad("PRINT_nuc.h5ad")
