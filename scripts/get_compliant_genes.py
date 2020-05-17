#!/usr/bin/env python3


'''
Script to get non-overlapping genes with at least two overlapping polyA sites
Takes a GENCODE GTF/GFF3 file as input
Outputs a GTF file only containing
'''

import sys
import pyranges as pyr
import pandas as pd


gtf = "/home/sam/paqr_annotations/tests/gencode.chr1.vM25.annotation.gtf"
polya_clusters = "/home/sam/paqr_annotations/atlas.clusters.2.0.GRCm38.96.bed"


def get_non_overlapping_genes(gtf_path=None):
    '''
    returns list of gene_ids of non-overlapping, protein-coding or lncRNA genes
    '''

    gtf_df = pyr.readers.read_gtf(f=gtf_path)
    # print(gtf_df.columns)
    gtf_df = gtf_df[gtf_df.Feature == 'gene']
    gtf_df = gtf_df[gtf_df.gene_type.isin(['protein_coding', 'lncRNA'])]

    # cluster function gives a common id to overlapping intervals
    # genes with common id/ id Count > 1 can be filtered out of pyranges
    gtf_df = gtf_df.cluster(strand=False, count=True)
    # print(gtf_df.columns)
    gtf_df = gtf_df[gtf_df.Count == 1]

    gene_id_list = gtf_df.gene_id.to_list()
    return gene_id_list


# print(get_non_overlapping_genes(gtf_path=gtf))

non_overlapping_genes = get_non_overlapping_genes(gtf_path=gtf)


def get_last_exons(ranges_obj):
    '''
    1. Select exons from Feature column
    2. Group by transcript_id
    3. Select last exon for each group (largest 'exon_number' value)
    '''
    df = ranges_obj.as_df()
    df = df[df['Feature'] == 'exon']

    # Select largest 'exon_number' row for each transcript_id - only solution I could get to work...
    # https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-count-in-groups-using-groupby
    df = df.sort_values('exon_number', ascending=False).drop_duplicates(['transcript_id'])

    return pyr.PyRanges(df)


def get_multi_polya_genes(gtf_path=None, subset_list=None, polya_bed_path=None):
    gtf_df = pyr.readers.read_gtf(f=gtf_path)
    # print(gtf_df.columns)
    gtf_df = gtf_df[gtf_df.gene_id.isin(subset_list)]

    polya_bed = pyr.readers.read_bed(f=polya_bed_path, as_df=True)
    # add chr prefix to Chromosome column
    polya_bed['Chromosome'] = 'chr' + polya_bed['Chromosome'].astype(str)

    polya_bed = pyr.PyRanges(polya_bed)
    print(polya_bed)
    gtf_last_exons = get_last_exons(gtf_df)
    print(gtf_last_exons)
    # 4. count_overlaps with BED file od polyA_sites
    gtf_last_exons = gtf_last_exons.count_overlaps(
        polya_bed, strandedness="same", keep_nonoverlapping=False)

    #print(gtf_last_exons[["transcript_id", "NumberOverlaps"]])

    # 5. Get list of transcript ids with at least two overlapping polyA_sites in terminal exon
    trs = gtf_last_exons[gtf_last_exons.NumberOverlaps >= 2]
    # print(trs)
    tr_list = set(trs.transcript_id.to_list())

    return tr_list


print(get_multi_polya_genes(gtf_path=gtf, subset_list=non_overlapping_genes, polya_bed_path=polya_clusters))
