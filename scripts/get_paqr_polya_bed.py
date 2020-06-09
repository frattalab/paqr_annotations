#!/usr/bin/env python3


'''
Script to get 'clusters' file for PAQR
BED file from polyAsite website has following format:
Chromosome | Start | End | <chr>:<coord>:<strand> | average tags PM | Strand | % samples supporting |
n_supporting_protocols | average tags PM | cluster_annotation | upstream_polyA_signal

Corresponding column names when importing BED with pyranges:
['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ThickStart',
       'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes'}

PAQR's clusters file requires a 10 column BED file in following format:
Chromosome | Start | End | <chr>:<strand>:<coord>:<cluster_annotation> | n_supporting_protocols/samples | Strand | site_n_on_exon |
total_n_on_exon | <transcript_id>:<exon_number>:<n_exons_on_tr>:<exon_start>:<exon_end> | gene_id

Column mapping:
<chr>:<coord>:<strand> = Name
average tags PM = Score
% samples supporting = ThickStart
n_supporting_protocols = ThickEnd
average tags PM (2) = ItemRGB
cluster_annotation = BlockCount
upstream_polyA_signal = BlockSizes

'''

import pyranges as pyr
import pandas as pd
import sys
from get_compliant_genes import get_last_exons


# tr_gtf_path = "/home/sam/paqr_annotations/tests/new_join_non_overlapping.multi_polya.gencode.chr1.vM25.annotation.gtf"
# polya_bed_path = "/home/sam/paqr_annotations/atlas.clusters.2.0.GRCm38.96.bed"
# atlas_version = 2
# out = "/home/sam/paqr_annotations/tests/new_join_new_script_clusters.non_overlapping.multi_polya.gencode.chr1.vM25.annotation.gtf"


def tidy_chromosome_column(polya_bed_path=None):
    '''
    Read in polyA bed file
    Add 'chr' prefix to chromosome column
    return pyranges object
    '''
    bed_df = pyr.read_bed(f=polya_bed_path, as_df=True)
    bed_df['Chromosome'] = 'chr' + bed_df['Chromosome'].astype(str)

    return pyr.PyRanges(bed_df)


def join_by_intersect(pyranges1=None, pyranges2=None):
    '''
    Joins two given pyranges object, reporting only entries which have same strand overlap
    columns of two objects are appended to one another (pyranges1 first)
    '''
    df = pyranges1.join(pyranges2, strandedness="same", how=None)
    return df


def add_paqr_name(pyranges=None, col_name='paqr_name'):
    '''
    Add a column to pyranges consisting of
    <chr>:<strand>:<coord>:<cluster_annotation>
    '''
    pyranges = pyranges.as_df()

    def f(x): return ':'.join([x['Chromosome'], x['Strand'],
                               x['Name'].split(':')[1], x['BlockCount']])
    pyranges[col_name] = pyranges.apply(f, axis=1)

    return pyr.PyRanges(pyranges)


def add_paqr_long_name(pyranges=None, col_name='paqr_long_name'):
    '''
    Add a column to pyranges consisting of
    <transcript_id>:<exon_number>:<n_exons_on_tr>:<exon_start>:<exon_end>
    '''
    pyranges = pyranges.as_df()

    def f(x): return ':'.join([x['transcript_id'], str(x['exon_number']),
                               str(x['exon_number']), str(x['Start_b']), str(x['End_b'])])
    pyranges[col_name] = pyranges.apply(f, axis=1)

    return pyr.PyRanges(pyranges)


def add_n_along_exon(pyranges=None, col_name='n_along_exon'):
    '''
    Add a column containing the consecutive poly(A) site number along the exon
    (1 = most proximal)
    + strand - End coordinate = most 3' (smallest value = most proximal site in gene)
    - strand - Start coordinate = most 3' (smallest value = most distal site in gene)
    '''
    df = pyranges.as_df()

    def sort_sites_along_exon(x):
        # first reset_index call removes the original index of the group (e.g. row 4005 in df)
        # second reset_index call adds the sorted index as a column to the dataframe (the order along exon in each transcript)
        if (x.Strand == '+').all():
            return x.sort_values(by=['End']).reset_index(drop=True).reset_index()
        elif (x.Strand == '-').all():
            return x.sort_values(by=['Start'], ascending=False).reset_index(drop=True).reset_index()

    # adds column called #index with consective order of sites along exon
    # reset_index call removes the 'transcript_id' index from the dataframe (i.e. the grouping)
    df = df.groupby('transcript_id').apply(sort_sites_along_exon).reset_index(drop=True)

    # indexes are 0 based - add one so 1st site in exon has value 1 etc.
    df['index'] = df['index'].add(1)

    df = df.rename(columns={'index': col_name})

    # .sort_values(by=['End']).cumcount()
    # print(df.loc['ENSMUST00000000266.8', ['Start', 'Strand']
    #             ])  # .reset_index(drop=True).reset_index())

    return pyr.PyRanges(df)


def get_total_n_on_exon(pyranges=None, col_name='total_n_on_exon'):
    '''
    Add a column containing total number of polyA sites on given exon
    '''
    df = pyranges.as_df()
    # print("nrows: %s" % (len(df.index)))

    # Get maximum n_along_exon for each transcript
    # Returns df of transcript_id | n_along_exon
    a = df.groupby('transcript_id')['n_along_exon'].max(
    ).reset_index().rename(columns={'n_along_exon': col_name})

    # inner join df with a (adds column to df for given transcript_id of total_n_sites)
    df = pd.merge(df, a, on='transcript_id')

    # a = df.groupby('transcript_id').apply(get_n_sites)
    # = df.groupby('transcript_id')[
    #    'n_along_exon'].max().reset_index().n_along_exon

    return pyr.PyRanges(df)


def write_to_paqr_bed(pyranges=None, outfile=None, col_order=['Chromosome', 'Start', 'End', 'paqr_name', 'ThickEnd', 'Strand', 'n_along_exon', 'total_n_on_exon', 'paqr_long_name', 'gene_id']):
    '''
    Need columns in following order:
    Chromosome | Start | End | <chr>:<strand>:<coord>:<cluster_annotation> | n_supporting_protocols/samples | Strand | site_n_on_exon |
    total_n_on_exon | <transcript_id>:<exon_number>:<n_exons_on_tr>:<exon_start>:<exon_end> | gene_id
    To do this I need to subset pyranges for following columns:

    Chromosome | Start | End | paqr_name | ThickEnd | Strand | n_along_exon | total_n_on_exon | paqr_long_name | gene_id

    (ignore) Note: col order should not contain 'Chromosome', 'Start', 'End' at start (can't work out why but output funky without it?!
    '''
    # pyranges only contains columns specified in col_order (Note: DOES NOT ORDER THEM)
    pyranges = pyranges[col_order]

    # pyranges columns are always sorted in follow order ('mandatory fields')
    # Chromosome | Start | End | Strand | <other column names>

    # with default col_order - columns in pyranges are returned in following order
    # paqr_name | Chromosome | Start | End | Strand | ThickEnd | n_along_exon' | 'total_n_on_exon' | 'paqr_long_name' | 'gene_id'

    # BED12 column names would be:
    # Chromosome Start End Name Score Strand ThickStart ThickEnd ItemRGB BlockCount BlockSizes BlockStarts

    # to change names of columns to match this naming convention, need a list in following order
    # Name | Chromosome | Start | End | Strand | Score | ThickStart | ThickEnd | ItemRGB | BlockCount
    # print(pyranges.columns)
    # Then can use inbuilt to_bed function
    # pyranges_bed_colnames = ['Name', 'Chromosome', 'Start', 'End', 'Strand',
    #                         'Score', 'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount']
    # pyranges.columns = pyranges_bed_colnames
    # print(pyranges)
    # print(pyranges.columns)
    # pyranges.to_bed(path=outfile, keep=True)

    df = pyranges.as_df()
    df = df[col_order]
    df.to_csv(outfile, index=False, header=False, sep='\t')


def process_polya_bed(polyA_bed_path=None, atlas_version=2, outfile=None):
    '''

    '''
    if atlas_version == 1:
        # No need to add 'chr' prefix
        polya_bed = pyr.readers.read_bed(f=polya_bed_path)

        # read in GTF with multi-overlap transcripts
        tr_gtf = pyr.read_gtf(f=tr_gtf_path)
        tr_gtf = get_last_exons(ranges_obj=tr_gtf)

        polya_bed = join_by_intersect(pyranges1=polya_bed, pyranges2=tr_gtf)

        # short name column already in correct format - skip to adding long name
        polya_bed = add_paqr_long_name(pyranges=polya_bed)

        # add columns with order along exon & number of exons on transcript
        polya_bed = add_n_along_exon(pyranges=polya_bed)
        polya_bed = get_total_n_on_exon(pyranges=polya_bed)

        # first 6 already provided in format of version 1 BED file
        column_order = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand',
                        'n_along_exon', 'total_n_on_exon', 'paqr_long_name', 'gene_id']

        write_to_paqr_bed(pyranges=polya_bed, outfile=outfile, col_order=column_order)

    elif atlas_version == 2:
        # need to add chr prefix before joining
        polya_bed = tidy_chromosome_column(polya_bed_path=polya_bed_path)

        # read in GTF with multi-overlap transcripts
        tr_gtf = pyr.read_gtf(f=tr_gtf_path)
        tr_gtf = get_last_exons(ranges_obj=tr_gtf)

        polya_bed = join_by_intersect(pyranges1=polya_bed, pyranges2=tr_gtf)

        # version 2 doesn't have name column in formatting required for PAQR - need to add both short & long name
        polya_bed = add_paqr_name(pyranges=polya_bed)
        polya_bed = add_paqr_long_name(pyranges=polya_bed)

        # add columns with order along exon & number of exons on transcript
        polya_bed = add_n_along_exon(pyranges=polya_bed)
        polya_bed = get_total_n_on_exon(pyranges=polya_bed)

        # v2.0 has custom format
        column_order = ['Chromosome', 'Start', 'End', 'paqr_name', 'ThickEnd',
                        'Strand', 'n_along_exon', 'total_n_on_exon', 'paqr_long_name', 'gene_id']

        write_to_paqr_bed(pyranges=polya_bed, outfile=outfile, col_order=column_order)


# process_polya_bed(polyA_bed_path=polya_bed_path, atlas_version=atlas_version, outfile=out)


if __name__ == '__main__':

    tr_gtf_path = sys.argv[1]
    polya_bed_path = sys.argv[2]
    atlas_version = int(sys.argv[3])
    out = sys.argv[4]

    process_polya_bed(polyA_bed_path=polya_bed_path, atlas_version=atlas_version, outfile=out)

'''
    polya_bed = tidy_chromosome_column(polya_bed_path=polya_bed_path)
    # print(polya_bed)

    tr_gtf = pyr.read_gtf(f=tr_gtf_path)
    tr_gtf = get_last_exons(ranges_obj=tr_gtf)

    polya_tr_df = join_by_intersect(pyranges1=polya_bed, pyranges2=tr_gtf)
    # print(polya_tr_df)
    # print(polya_tr_df.columns)

    polya_tr_df = add_paqr_name(pyranges=polya_tr_df)
    polya_tr_df = add_paqr_long_name(pyranges=polya_tr_df)
    # print(polya_tr_df[['paqr_long_name']])
    # print(polya_tr_df.columns)

    polya_tr_df = add_n_along_exon(pyranges=polya_tr_df)

    # print(polya_tr_df[['transcript_id', 'n_along_exon']].head(n=8))
    # print(polya_tr_df.columns)

    polya_tr_df = get_total_n_on_exon(pyranges=polya_tr_df)
    # print(polya_tr_df[['total_n_on_exon']])
    # print(polya_tr_df.columns)

    write_to_paqr_bed(pyranges=polya_tr_df, outfile=out)
    '''
