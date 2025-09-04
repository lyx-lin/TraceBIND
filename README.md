# TraceBIND

- [Overview](#overview)
- [Installation](#installation)
- [Tutorials](#tutorials)
- [License](#license)
  
<img src="https://github.com/lyx-lin/TraceBIND/blob/main/figures/tracebind_overview.png">

## Overview
The identification of footprints provides a powerful means to detect base-pair resolution signals of regulatory elements in chromatin accessibility data. We observe substantial variability in Tn5 transposase cleavage bias across individual samples, highlighting the need for sample-specific correction. TraceBIND is a R package that that corrects sample-specific Tn5 bias through mitochondria-based fine-tuning of PRINT’s deep learning model, and use it to identify TF and nucleosome footprints through a dynamic flanking window statistical scan. TraceBIND also enables sample-specific FDR-controlled p-value thresholds stratified by varying coverage, which is necessary because coverage can vary across orders of magnitude in ATAC-seq data.

## Installation
Before installing this package, make sure you have the following dependencies installed:
```bash
# install.packages("devtools")
devtools::install_github("yaowuliu/ACAT")
```

Then install this package:
```bash
# install.packages("devtools")
devtools::install_github("lyx-lin/TraceBIND", dependencies=TRUE)
```

Our Python package has been tested on python=3.11, 3.12. The requirements of python packages for traceBIND finetuning are listed in the [requirements](https://github.com/lyx-lin/TraceBIND/blob/main/requirements.txt), which can be done by:
```bash
pip install -r requirements.txt
```

## Tutorials 
For training a mitochondrial-finetuned sample-specific Tn5 bias model [Finetuning](https://github.com/lyx-lin/TraceBIND/blob/main/tutorials/tutorial_finetuning.ipynb). 
Before finetuning, please download [PRINT Tn5 bias model](https://github.com/HYsxe/PRINT/blob/main/data/shared/Tn5_NN_model.h5).

For identification of footprints [Identification](https://github.com/lyx-lin/TraceBIND/blob/main/tutorials/tutorial_footprint_identification.ipynb). 

For footprint-informed chromVar [footprint_chromvar](https://github.com/lyx-lin/TraceBIND/blob/main/tutorials/tutorial_chromvar.ipynb). 

The data used could be found [here](https://www.dropbox.com/scl/fo/zhmxfp0gxnmlgeo8jsmbv/AO3I75Lz6eP3Illn-eb0Zgc?rlkey=zkfi6c7c29eb11tbmcz80n8sf&st=2cstifvu&dl=0).
