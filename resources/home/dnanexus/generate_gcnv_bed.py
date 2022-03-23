"""
    Script to generate bed file from a run of cnv calls and intervals list.

    Expected inputs:
        - denoised copy ratio tsv files of samples to visualise
        - run name prefix

    Requires bgzip and tabix to be installed and on path.

    created: 211119 by Jethro Rainford
    last modified: 220323 by Sophie Ratkai
"""
import argparse
from pathlib import Path
import subprocess
from numpy import inner

import pandas as pd


def parse_args():
    """
        Parse command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--sample_name', help='name of sample to prefix output bed file with'
    )
    parser.add_argument(
        '--sample',
        help='denoised copy ratio file to generate bed file from'
    )
    parser.add_argument(
        '--mean_std',
        help='mean+std copy ratio values to generate bed file from'
    )

    args = parser.parse_args()

    return args


def write_outfile(copy_ratio_df, sample):
    """
        Write output bed files and compress with bgzip.
        Bed file with highlight_samples has random colouring of each sample
        trace to improve visibility.

        Args:
            - copy_ratio_df (df): df of all copy ratios to write
            - prefix (str): prefix for naming output file
    """
    outfile = "{}_copy_ratios.gcnv.bed".format(sample)

    with open(outfile, 'w') as fh:
        # colour mapping for tracks
        colours = {
            sample: 'red',
            'mean': 'blue',
            'plus_std': '#183B59',
            'minus_std': '#183B59',
            'plus_std2': '#317A86',
            'minus_std2': '#317A86',
        }

        highlight = ' '.join([f"highlight={x};{y}" for x, y in colours.items()])
        # this line is needed at the beginning to tell igv.js that it is a gcnv
        # bed as it can't automatically set this track type
        fh.write(f"track  type=gcnv height=500 {highlight} \n")

    copy_ratio_df.to_csv(outfile, sep='\t', index=False, mode='a')

    # compress & index
    subprocess.call("bgzip {}".format(Path(outfile)), shell=True)
    subprocess.call("tabix -f -S 2 -b 2 -e 3 {}.gz".format(Path(outfile)), shell=True)


def main():

    args = parse_args()
    # print(args.sample_name)

    # Read in the sample's copy ratio file
    sample_copy_ratio_df = pd.read_csv(
        args.sample, sep='\t', comment='@', header=0,
        names=['chr', 'start', 'end', args.sample_name]
    )
    intervals_df = sample_copy_ratio_df.iloc[:, :3].copy()

    # check if mean copy ratio file is provided as tsv or csv
    mean_std_copy_ratio_df = pd.read_csv(args.mean_std, sep='\t')
    if len(mean_std_copy_ratio_df.columns) < 2:
        mean_std_copy_ratio_df = pd.read_csv(args.mean_std)

    assert sample_copy_ratio_df.iloc[:, :3].equals(
        mean_std_copy_ratio_df.iloc[:, :3]), (
        f"Intervals do not match for sample and run mean files")

    copy_ratio_df = pd.merge(sample_copy_ratio_df, mean_std_copy_ratio_df,
                            on=['chr', 'start', 'end'])

    # write output bed file
    write_outfile(copy_ratio_df, args.sample_name)


if __name__ == "__main__":
    main()
