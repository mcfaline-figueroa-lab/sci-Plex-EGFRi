# sci-Plex-EGFRi

[sci-Plex-EGFRi](https://www.biorxiv.org/content/10.1101/2024.04.08.587960v1) is a study to elucidate the heterogenous pharmaco-transcriptomic landscape induced by targeting EGFR in models of glioblastoma. This github repository containing processing pipelines and analysis scripts, including files used for data generation, analysis and interpretation.

### Contents

##### Data and data processing scripts
The data presented in our manuscript are derived from 3 single-cell screens. All raw and processed files can be found on the National Center for Biotechnology Information Gene Expression Ominibus (NCBI GEO) repository under series GSE261618.

| Name        | Experiment           |GEO Accession  | Analysis Scripts |
| :-------------: |:-----------:| :----:| :---:|
| sci-Plex-EGFRi-ATCCpilot      | Pilot Screen | [GSM8147508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8147508)| sci-Plex-EGFRi_pilotscreen|
| sci-Plex-EGFRi-PDCLpilot| Pilot Screen      |  [GSM8147509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8147509) | sci-Plex-EGFRi_pilotscreen|
| sci-Plex-EGFRi-largescreen | Large Screen      | [GSM8147510](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8147510) | sci-Plex-EGFRi_largescreen |

##### Data processing and analysis
The scripts and code to process the data can be found in the [process\_from_raw](https://github.com/mcfaline-figueroa-lab/sci-Plex-EGFRi/tree/main/process_from_raw) folder. The scripts and code to analyze the data can be found in the folders for each experiment.
