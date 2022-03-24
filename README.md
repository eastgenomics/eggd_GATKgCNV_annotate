# Filter and annotate single sample CNV vcfs (DNAnexus Platform App)

Dx wrapper to run a custom script that filters CNV calls from the GATK _segments.vcf output file to panel regions and annotates the filtered CNV calls with gene symbol, transcript and exon number information.

This is the source code for an app that runs on the DNAnexus Platform.

## What does this app do?
Filters the GATK _segments.vcf file (output from the CNV calling app) with bcftools.
Annotates the filtered CNV vcf with gene symbol, transcript and exon number information.

## What are typical use cases for this app?
This app may be executed as the last step of the germline CNV calling workflow, after CNV calling of samples on a run.

## What data are required for this app to run?
This app requires:
* panel.bed file (from generate_bed)
* sample_segments.vcf file (from CNV_calling)
* exons BED file: sorted in chromosome order and containing the following columns: chromosome number (no 'chr'), start and end coordinates, gene symbol, (RefSeq) transcript ID and exon number (numerical)


Optional inputs:
* tiVisualise can be set TRUE/FALSE depending on whether sample copy ratios should be plotted against the run mean +/- standard deviation in a gCNV.bed file. If required, the sample's denoised copy ratios and the run mean and std should be provided as input files.
* toAnnotateVEP can be set TRUE/FALSE depending on whether calls should be annotated with information from VEP

## What does this app output?
* tsv of the filtered and annotated CNVs called in the sample
    * an empty file is provided if no CNVs were detected in the panel regions
* a sample gCNV.bed, if requested

## Dependencies
The app uses the htslib asset (stored in 001_Reference)
The app might include VEP or something in the future.

### This app was made by East GLH
