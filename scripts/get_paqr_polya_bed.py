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
import sys
from get_compliant_genes import get_last_exons


tr_gtf_path = "/home/sam/paqr_annotations/tests/non_overlapping.multi_polya.gencode.chr1.vM25.annotation.gtf"
polya_bed_path = "/home/sam/paqr_annotations/atlas.clusters.2.0.GRCm38.96.bed"
out = "/home/sam/paqr_annotations/tests/clusters.non_overlapping.multi_polya.gencode.chr1.vM25.annotation.gtf"


def tidy_chromosome_column(polya_bed_path=None):
    '''
    Read in polyA bed file
    Add 'chr' prefix to chromosome column
    return pyranges object
    '''
    bed_df = pyr.read_bed(f=polya_bed_path, as_df=True)
    bed_df['Chromosome'] = 'chr' + bed_df['Chromosome'].astype(str)

    return pyr.PyRanges(bed_df)


polya_bed = tidy_chromosome_column(polya_bed_path=polya_bed_path)
# print(polya_bed)

tr_gtf = pyr.read_gtf(f=tr_gtf_path)
tr_gtf = get_last_exons(ranges_obj=tr_gtf)


def join_by_intersect(pyranges1=None, pyranges2=None):
    '''
    Joins two given pyranges object, reporting only entries which have same strand overlap
    columns of two objects are appended to one another (pyranges1 first)
    '''
    df = pyranges1.join(pyranges2, strandedness=None, how=None)
    return df


polya_tr_df = join_by_intersect(pyranges1=polya_bed, pyranges2=tr_gtf)
# print(polya_tr_df)
# print(polya_tr_df.columns)


def add_paqr_name_formatting(pyranges=None, col_name='paqr_name'):
    '''
    Add a column to pyranges consisting of
    <chr>:<strand>:<coord>:<cluster_annotation>
    '''
    pyranges = pyranges.as_df()

    def f(x): return x['Chromosome']+":" + x['Strand'] + \
        ":" + x['Name'].split(':')[1] + ":" + x['BlockCount']
    pyranges[col_name] = pyranges.apply(f, axis=1)

    # pyranges[col_name] = pyranges['Chromosome'].astype(str)+":" + pyranges['Strand'].astype(
    #    str) + ":" + split_by_colon(str=pyranges['Name'].astype(str), idx=1) + ":" + pyranges['BlockCount'].astype(str)

    return pyr.PyRanges(pyranges)


polya_tr_df = add_paqr_name_formatting(pyranges=polya_tr_df)
print(polya_tr_df[['paqr_name']])
print(polya_tr_df.columns)
