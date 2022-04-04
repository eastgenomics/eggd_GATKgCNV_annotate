"""
    Script to parse CNV vcf file and annotate each call with exon information.
    
    Expected inputs:
        - a sample.vcf file
        - an exons.bed file for annotation with expected columns
            ["chrom", "start", "end", "gene", "transcript", "exon_num"]

    created: 220317 Sophie Ratkai
    last modified: 220404 Sophie Ratkai
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
    call_exon_df = calls_w_exons.to_dataframe(index_col=False,
        names=["call_chrom", "call_start", "call_end", "ID",
        "exon_chrom", "exon_start", "exon_end", "gene", "transcript", "exon"])
    call_exon_df['exon_length'] = \
        call_exon_df['exon_end'] - call_exon_df['exon_start']
    # print("call_exon_df")
    # print(call_exon_df)

    unique_calls = call_exon_df["ID"].unique().tolist()  # segments
    genes = []  # list that will become the genes column in the df
    transcripts = []
    exons = []  # list that will become the exons column in the df
    exonic_starts = []
    exonic_ends = []
    lengths = []

    for call in unique_calls:  # Depending on the number of exon annotation and
                        # whether there's annotation available at all:
        annotation_df = call_exon_df[call_exon_df["ID"] == call]
        annotation_df = annotation_df.reset_index(drop=True)
        # print("annotation_df for sample {}, call {}".format(
        #     calls_df["sample"][0], call))
        if len(annotation_df) == 1:  # call is a single exon
            call_gene = annotation_df["gene"][0]
            call_transcript = annotation_df["transcript"][0]
            call_exon = annotation_df["exon"][0]
            exonic_start = annotation_df["exon_start"][0]
            exonic_end = annotation_df["exon_end"][0]
            call_length = annotation_df["exon_length"][0]
        else:  # this call covers multiple exons
            call_genes = annotation_df["gene"].unique().tolist()
            if len(call_genes) == 1:  # call covers multiple exons in 1 gene
                call_gene = call_genes[0]
                call_transcript = annotation_df["transcript"].tolist()[0]
                exon_nums = annotation_df["exon"].unique().tolist()
                start_exon = str(exon_nums[0])
                end_exon = str(exon_nums[-1])
                call_exon = "-".join([start_exon, end_exon])
                exonic_start = annotation_df["exon_start"][0]
                exonic_end = annotation_df["exon_end"].to_list()[-1]
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
                    exon_nums = exon_annot_df["exon"].unique().tolist()
                    if len(exon_nums) == 1:
                        call_exons.append(str(exon_nums[0]))
                    else:
                        start_exon = str(exon_nums[0])
                        end_exon = str(exon_nums[-1])
                        call_exons.append("-".join([start_exon, end_exon]))
                    exon_lengths.append(exon_annot_df['exon_length'].sum())
                call_exon = ', '.join(call_exons)
                exonic_start = annotation_df["exon_start"][0]
                exonic_end = annotation_df["exon_end"][0]
                call_length = sum(exon_lengths)

        genes.append(call_gene)
        transcripts.append(call_transcript)
        exons.append(call_exon)
        exonic_starts.append(exonic_start)
        exonic_ends.append(exonic_end)
        lengths.append(call_length)

    # create annotation dataframe
    CNV_annotation = pd.DataFrame.from_dict({"ID": unique_calls, "gene": genes,
                "transcript": transcripts, "exon": exons, "length": lengths,
                "exonic_start": exonic_starts, "exonic_end": exonic_ends})
    return CNV_annotation


if __name__ == "__main__":

    # Parse command inputs: filtered vcf and exons.tsv
    vcf_file = sys.argv[1]
    panel_bed = sys.argv[2]
    panel_name = panel_bed.strip('.bed')
    
    # Load the panel bed file into a dataframe
    exons_df = pd.read_csv(panel_bed, sep='\t')
    # Ensure that panel bed also has the annotation information
    assert len(exons_df.columns) == 6, (
        f"Panel bed does not have annotation information"
    )
    exons_df.columns = ["chrom", "start", "end", "gene", "transcript", "exon"]

    # Parse VCF file into a dataframe
    sample_name, CNV_calls = vcf2df(vcf_file)
    # Ensure that filtered_vcf has variants, otherwise create empty tsv
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
    outfile = "{}_{}_annotated_CNVs.tsv".format(sample_name, panel_name)
    with open(outfile, 'w') as fh:
        fh.write(f"Annotated CNV calls for sample {sample_name} in {panel_name} panel \n")

    # Add annotation: "gene","transcript","exon","length" columns
    call_annotation = annotate(
        CNV_calls[["chrom", "start", "end", "ID"]], exons_df
    )
    annotated_CNVs = CNV_calls.merge(call_annotation, on='ID')

    # display exonic range in a format that can be copied into IGV or databases
    annotated_CNVs['exonic range'] = annotated_CNVs.apply(
        lambda row: "-".join(
            [str(row.chrom), str(row.exonic_start), str(row.exonic_end)]),
        axis=1, result_type='expand'
    )

    # Modify the datatpye, order and name of columns for the output
    annotated_CNVs[['NP', 'QA', 'QS', 'QSE', 'QSS']] = \
        annotated_CNVs[['NP', 'QA', 'QS', 'QSE', 'QSS']].astype("int64")
    annotated_CNVs = annotated_CNVs.reindex(
        columns=['chrom', 'start', 'end', 'NP', 'CN', 'CNV',
                'QA', 'QS', 'QSS', 'QSE', 'gene', 'transcript', 'exon',
                'exonic range', 'length'])
    annotated_CNVs.columns = ['chromosome', 'segment start', 'segment end',
            '# intervals', 'copy number', 'CNV',
            'QA', 'QS', 'QS start', 'QS end',
            'gene', 'transcript', 'exon(s)', 'exonic range', 'affected bases']

    # Write annotated calls to file
    annotated_CNVs.to_csv(outfile, sep='\t', header=True, index=False,
        encoding='utf-8', mode='a'
    )
