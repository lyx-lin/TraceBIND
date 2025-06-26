# TraceBIND
## overview
The identification of footprints provides a powerful means to detect base-pair resolution signals of regulatory elements in chromatin accessibility data. We observe substantial variability in Tn5 transposase cleavage bias across individual samples, highlighting the need for sample-specific correction. TraceBIND is a R package that that corrects sample-specific Tn5 bias through mitochondria-based fine-tuning of PRINT’s deep learning model, and use it to identify TF and nucleosome footprints through a dynamic flanking window statistical scan. TraceBIND also enables sample-specific FDR-controlled p-value thresholds stratified by varying coverage, which is necessary because coverage can vary across orders of magnitude in ATAC-seq data.
## Tutorials 
For training a mitochondrial-finetuned sample-specific Tn5 bias model [Finetuning](https://github.com/lyx-lin/TraceBIND/blob/main/tutorial/tutorial_finetuning.ipynb). 

For identification of footprints [Identification](https://github.com/lyx-lin/TraceBIND/blob/main/tutorial/tutorial_footprint_identification.ipynb). 

For footprint-informed chromVar [footprint_chromvar](https://github.com/lyx-lin/TraceBIND/blob/main/tutorial/tutorial_chromvar.ipynb). 

