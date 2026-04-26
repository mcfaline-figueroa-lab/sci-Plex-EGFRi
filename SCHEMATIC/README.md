**SCHEMATIC** (**S**ingle-cell **CHEM**ic**A**l **T**ranscriptom**I**cs-based chemotype **C**lassification)

Reusable template for running an MRVI-based single-cell transcriptomic classification of inhibitors. 
Configure paths and metadata columns once, then run the workflow end-to-end or resume from saved intermediate outputs.
Generalizable but some factors may be experiment specific e.g. how many cell models, inhibitors, concentrations, etc.

**Expected inputs**
- AnnData .h5ad file, or a Matrix Market counts file plus gene metadata and cell metadata.
- Cell metadata columns identifying biological group/sample, batch/replicate, and optionally labels/covariates.

**Dependencies**
- scvi-tools
- scanpy
- pandas
- numpy
- scipy
- seaborn
- matplotlib
