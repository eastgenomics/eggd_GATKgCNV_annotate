"""
    Script to parse CNV vcf file and annotate each call with exon information.
    
    Expected inputs:
        - a sample.vcf file
        - an exons.bed file for annotation with expected columns
            ["chrom", "start", "end", "gene", "transcript", "exon_num"]

    Sophie Ratkai 220317
"""
import sys

import pandas as pd
import vcf
import pybedtools as bedtools


def vcf2df(vcf_file):
    CNVcalls = pd.DataFrame()

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        CNV_call = {}
        call_sample = record.samples[0]  # vcf.model._Call
        CNV_call['sample'] = call_sample.sample
        CNV_call['chrom'] = record.CHROM  # str
        CNV_call['start'] = record.POS  # int
        CNV_call['end'] = record.INFO["END"]
        CNV_call['ID'] = record.ID  # str
        format_fields = record.FORMAT.split(':')  # GT:CN:NP:QA:QS:QSE:QSS
        for i, field in enumerate(format_fields):
            CNV_call[field] = call_sample.data[i]
        if int(CNV_call['CN']) == 2:
            continue  # skip normal regions
        elif int(CNV_call['CN']) < 2:
            CNV_call['CNV'] = 'DEL'
        else:
            CNV_call['CNV'] = 'DUP'
        CNVcalls = CNVcalls.append(CNV_call, ignore_index=True).astype(
                    {'start': 'int64', 'end': 'int64', 'CN': 'int64'})
    return CNV_call['sample'], CNVcalls


def annotate(calls_df, exons_df):
    """
        Annotate each CNV call with overlapping exon information
    """
    calls_bed = bedtools.BedTool.from_dataframe(calls_df)
    exon_bed = bedtools.BedTool.from_dataframe(exons_df)

    # Intersecting calls bed file with transcript_exon information
    # requires 100% overlap on exon within call coordinates
    calls_w_exons = calls_bed.intersect(exon_bed, loj=True)
    # convert pybedtools object to dataframe
    annotated_calls_df = calls_w_exons.to_dataframe(index_col=False,
        names=["call_chrom", "call_start", "call_end", "call_ID",
        "exon_chrom", "exon_start", "exon_end", "gene", "transcript", "exon"])
    annotated_calls_df['exon_length'] = \
        annotated_calls_df['exon_end'] - annotated_calls_df['exon_start']
    # print("annotated_calls_df")
    # print(annotated_calls_df)

    calls = annotated_calls_df["call_ID"].unique().tolist()
    genes = []  # list that will become the genes column in the dfexon
    transcripts = []
    exons = []  # list that will become the exons column in the df
    lengths = []

    for call in calls:  # Depending on the number of exon annotation and
                        # whether there's annotation available at all:
        annotation_df = annotated_calls_df[annotated_calls_df["call_ID"] == call]
        annotation_df = annotation_df.reset_index(drop=True)
        # print("annotation_df for sample {}, call {}".format(
        #     calls_df["sample"][0], call))
        if len(annotation_df) == 1:  # call is a single exon
            call_gene = annotation_df["gene"][0]
            call_transcript = annotation_df["transcript"][0]
            call_exon = annotation_df["exon"][0]
            call_length = annotation_df["exon_length"][0]
        else:  # this call covers multiple exons
            call_genes = annotation_df["gene"].unique().tolist()
            if len(call_genes) == 1:  # call covers multiple exons in 1 gene
                call_gene = call_genes[0]
                call_transcript = annotation_df["transcript"].tolist()[0]
                exon_nums = sorted(annotation_df["exon"].unique().tolist())
                start_exon = str(exon_nums[0])
                end_exon = str(exon_nums[-1])
                call_exon = "-".join([start_exon, end_exon])
                call_length = sum(annotation_df['exon_length'].tolist())
            else:  # this call covers multiple exons in multiple genes
                call_gene = ', '.join(call_genes)
                call_transcripts = annotation_df["transcript"].unique().tolist()
                call_transcript = ', '.join(call_transcripts)
                call_exons = []
                exon_lengths = []
                for transcript in call_transcripts:
                    exon_annot_df = annotation_df[
                        annotation_df['transcript'] == transcript][
                            ['exon', 'exon_length']]
                    exon_nums = sorted(exon_annot_df["exon"].unique().tolist())
                    if len(exon_nums) == 1:
                        call_exons.append(str(exon_nums[0]))
                    else:
                        start_exon = str(exon_nums[0])
                        end_exon = str(exon_nums[-1])
                        call_exons.append("-".join([start_exon, end_exon]))
                    exon_lengths.append(exon_annot_df['exon_length'].sum())
                call_exon = ', '.join(call_exons)
                call_length = sum(exon_lengths)

        genes.append(call_gene)
        transcripts.append(call_transcript)
        exons.append(call_exon)
        lengths.append(call_length)

    # create annotation dataframe
    CNV_annotation = pd.DataFrame.from_dict({"ID": calls, "gene": genes,
                "transcript": transcripts, "exon": exons, "length": lengths})
    return CNV_annotation


if __name__ == "__main__":

    # Parse command inputs: filtered vcf and exons.tsv
    vcf_file = sys.argv[1]
    exons_tsv = sys.argv[2]
    panel_name = sys.argv[3]

    # Parse VCF file
    sample_name, CNV_calls = vcf2df(vcf_file)
    # Ensure that filtered_vcf has variants, otherwise exit
    if len(CNV_calls) < 1:
        print(f"Sample {sample_name} has no CNVs")
        # No CNVs detected in the panel regions
        line = (
            f"There are no CNVs in the regions listed in {panel_name}"
            f" in sample {sample_name}"
        )
        with open(sample_name + "_annotated_CNVs.tsv", "w") as fh:
            fh.write(line)
        sys.exit(0)

    # else: CNV calls found in panel regions

    # Load the exon annotation bed file into a dataframe
    exons_df = pd.read_csv(exons_tsv, sep='\t', names=["chrom", "start", "end",
        "gene", "transcript", "exon_num"])

    # Add annotation: "gene","transcript","exon","length" columns
    CNV_annotation = annotate(
        CNV_calls[["chrom", "start", "end", "ID"]], exons_df
    )
    annotated_CNVs = CNV_calls.merge(CNV_annotation, on='ID')

    # Write annotated calls and counts to files
    annotated_CNVs.to_csv(sample_name + "_annotated_CNVs.tsv", sep='\t',
        encoding='utf-8', header=True, index=False
    )
