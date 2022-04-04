#!/bin/bash
# gCNVcaller_annotate

# Exit at any point if there is any error and output each line as it is executed (for debugging)
set -e -x -o pipefail

main() {

    echo "Installing packages"
    cd packages
    pip install -q pytz-* python_dateutil-* pysam-* numpy-* pandas-* pybedtools-* PyVCF-*
    cd ..

    # Download input files
    echo "Downloading input files"
    dx download "$sample_vcf" -o sample.vcf
    dx download "$panel_bed" -o panel.bed

    # A. Filter sample VCF with panel bed file
    echo "Filtering the sample VCF"
    # take the _intervals.vcf output from GATK and zip then index
    bgzip -c sample.vcf > sample.vcf.gz
    tabix -fp vcf sample.vcf.gz
    # and filter with the panel BED file
    bcftools filter -R panel.bed sample.vcf.gz > filtered.vcf

    # B. Annotate CNV calls with gene, transcript and exon number information
    echo "Annotating the filtered CNV calls"
    python3 annotate_calls.py filtered.vcf panel.bed

    # Create folders for the output files:
    fv=out/filtered_vcf/ && mkdir -p ${fv}
    cp filtered.vcf ${fv}/"${sample_vcf_prefix}_filtered.vcf"

    at=out/annotated_tsv/ && mkdir -p ${at}
    cp *_annotated_CNVs.tsv ${at}

    # C. Run VEP annotation:
    # bcftools   -e'GT ="0"' -f'%CHROM %POS [ %GT]\n' sample.vcf
    # bcftools annotate -a exon_annot.bed -c CHROM,POS,INFO/END,gene,transcript,exon_num -o annotated_intervals.vcf  sample_intervals.vcf

    echo "All scripts finished successfully, uploading output files to DNAnexus"

    # Upload results
    dx-upload-all-outputs --parallel
}
