#!/usr/bin/env python3


'''
Script to get 'clusters' file for PAQR
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
print(polya_tr_df)
print(polya_tr_df.columns)
