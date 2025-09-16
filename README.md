# Real_etal_2025_Predator_Stress
Reproducibility code for VF Real, Madalena, et al. "The mouse gut microbiota responds to predator odor and predicts host behavior." _bioRxiv_ (2025): 2025-07. https://doi.org/10.1101/2025.07.01.662568 

## Raw data
Gut metagenomic and visceral adipose tissue (VAT) transcriptomic sequence data can be found in NCBI’s Sequence Read Archive under accession no. [PRJNA1271077](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1271077).

## RNAseq_data_processing.sh
**Bash script** used to quality control and map **RNAseq reads** to the _Mus musculus domesticus_ reference genome, in order to get a ReadsPerGene STAR output table that can be input into **Real_2025_analyses.Rmd** for further analyses.

## Metagenomic_data_processing.sh
**Bash script** used to quality control, assemble, and bin **gut metagenomic reads** into metagenome-assembled genomes (MAGs) using [MAGmaker](https://github.com/CUMoellerLab/sn-mg-pipeline), our lab's custom made snakemake pipeline. Outputs can be used in **Real_2025_analyses.Rmd** for further analyses.

## Stress_analyses.Rmd
**R markdown** file with code used to produce the figures and tables included in Real et al. 2025.

## RF_predict_phenotypes.py & RF_wrapper.sh
Use the **bash script wrapper** to run a **Python script** that uses **random forest models** to quantify the predictive power of the gut microbiota or the VAT transcriptome in regards to host behavior. Outputs can be used in **Stress_analyses.Rmd** for further analyses.
