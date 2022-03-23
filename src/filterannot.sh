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
    dx download "$exon_list" -o exons.bed

    # A. Filter sample VCF with panel bed file
    echo "Filtering the sample VCF"
    # take the _intervals.vcf output from GATK and zip then index
    bgzip -c sample.vcf > sample.vcf.gz
    tabix -fp vcf sample.vcf.gz
    # and filter with the panel BED file
    bcftools filter -R panel.bed sample.vcf.gz > filtered.vcf

    # B. Annotate CNV calls with gene, transcript and exon number information
    echo "Annotating the filtered CNV calls"
    python3 annotate_calls.py filtered.vcf exons.bed "$panel_bed_prefix"

    # Create folders for the output files:
    fv=out/filtered_vcf/ && mkdir -p ${fv}
    cp filtered.vcf ${fv}/"${sample_vcf_prefix}_filtered.vcf"

    at=out/annotated_tsv/ && mkdir -p ${at}
    cp *_annotated_CNVs.tsv ${at}

    if $toVisualise; then
        if [[ ! -z $sample_copy_ratios ]]; then
            if [[ ! -z $mean_copy_ratios ]]; then
                echo "Sample copy ratios and mean & std are provided"
                dx download "$sample_copy_ratios" -o sample_copy_ratios.tsv
                dx download "$mean_copy_ratios" -o mean_copy_ratios.tsv
                echo "Generating gcnv bed file for sample"
                python3 sample_gcnv_bed.py --panel "$panel_bed_prefix" \
                --sample sample_copy_ratios.tsv --mean_std mean_copy_ratios.tsv
            else
                echo "Run mean copy ratio file was not provided"
            fi
        else
            echo "Sample copy ratio file was not provided"
        fi
    fi

    vf=out/visualisation_files/ && mkdir -p ${vf}
    cp *_annotated_CNVs.tsv ${vf}

    # C. Run VEP annotation:
    # for vcf_file in inputs/*_segments.vcf; do
    #     bcftools query -e'GT ="0"' -f'%CHROM %POS [ %GT]\n' $vcf_file
    # done

    echo "All scripts finished successfully, uploading output files to DNAnexus"

    # Upload results
    dx-upload-all-outputs --parallel
}
