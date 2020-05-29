#!/usr/bin/env python3


'''
Script to get non-overlapping genes with at least two overlapping polyA sites
Takes a GENCODE GTF/GFF3 file as input
Outputs a GTF file containing only transcript_ids with at least two overlapping polyA sites in last exon
'''

import sys
import pyranges as pyr
import pandas as pd


gtf = "/home/sam/paqr_annotations/tests/gencode.chr1.vM25.annotation.gtf"
polya_clusters = "/home/sam/paqr_annotations/atlas.clusters.2.0.GRCm38.96.bed"
atlas_version = 2
out = "/home/sam/paqr_annotations/tests/new_join_non_overlapping.multi_polya.gencode.chr1.vM25.annotation.gtf"
gtf_pyranges = pyr.readers.read_gtf(f=gtf)


def get_last_exons(ranges_obj):
    '''
    Returns last exon for each transcript
    1. Select exons from Feature column
    2. Group by transcript_id
    3. Select last exon for each group (largest 'exon_number' value)
    '''
    df = ranges_obj.as_df()
    df = df[df['Feature'] == 'exon']

    # Select largest 'exon_number' row for each transcript_id - only solution I could get to work...
    # https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-count-in-groups-using-groupby
    # sort by exon number value - last exon = top of sort
    # drop_duplicates keeps top row for each transcript_id - removes all rows but last exon
    df = df.sort_values('exon_number', ascending=False).drop_duplicates(['transcript_id'])

    return pyr.PyRanges(df)


def get_non_overlapping_genes(gtf_df=None):
    '''
    returns list of transcript_ids of non-overlapping, protein-coding or lncRNA genes
    non-overlapping = last exon of transcript doesn't overlap with coordinates of
    '''

    # print(gtf_df.columns)

    gtf_df = gtf_df[gtf_df.gene_type.isin(['protein_coding', 'lncRNA'])]

    # get last exon for each transcript
    exons = get_last_exons(ranges_obj=gtf_df)

    exons = exons.as_df()

    # print(gtf_df.dtypes)
    # remove entries with a transcript support level tag of 'NA' or 'NaN'
    # gtf_df = gtf_df.dropna(axis=0, subset=['transcript_support_level'])
    exons = exons[pd.to_numeric(exons['transcript_support_level'], errors='coerce').notnull()]

    # print(gtf_df.drop(gtf_df.columns[[22, 23, 24]], axis=1))
    # print(gtf_df.drop(gtf_df.columns[[22, 23, 24]], axis=1).groupby(
    #    'gene_id').get_group('ENSMUSG00000026131.20'))
    # print(gtf_df.groupby('gene_id').get_group('ENSMUSG00000000544.14'))

    def select_best_tsl_isoforms(x):
        '''
        gene_id has multiple isoforms with different last exons
        Only keeps isoforms with the highest transcript_support_level
        e.g. if gene has 4 isoforms, 3 have TSL of 1 and 1 of TSL 3
        isoform with TSL 3 is dropped (3 rows/isoforms returned)
        '''
        # only need to filter gene_id groups with > 1 transcript for last exon
        idx = x.groupby('gene_id')['transcript_support_level'].transform(
            min) == x['transcript_support_level']
        return x[idx]

        # sorts ascending order - TSL: 1 is best supported
        # return x.sort_values('transcript_support_level').drop_duplicates(['transcript_id'])

    # Get best supported transcripts for each gene (transcript isoforms can have same last exon)
    exons = select_best_tsl_isoforms(x=exons)
    # print(gtf_df.drop(gtf_df.columns[[22, 23, 24]], axis=1))
    # print(testing.drop(testing.columns[[22, 23, 24]], axis=1).groupby(
    #    'gene_id').get_group('ENSMUSG00000026131.20'))

    # cluster function gives a common id to overlapping intervals
    # genes with common id/ id Count > 1 can be filtered out of pyranges
    # gtf_df = pyr.PyRanges(gtf_df)
    # gtf_df = gtf_df.cluster(strand=False, count=True)

    # print(gtf_df[['gene_id', 'Cluster', 'Count']])

    # Interest is finding exons that can be unambigously assigned to a SINGLE gene
    # cluster function may give common ID to isoforms of the same gene that overlap in last exon
    # if cluster group contains more than one gene it should be filtered out (evaluate False in filter)

    # gtf_df = gtf_df.as_df()
    # gtf_df = gtf_df.groupby('Cluster').filter(lambda x: len(list(set(x['gene_id']))) == 1)
    # print(gtf_df)
    # gtf_df = gtf_df[gtf_df.Count == 1]

    exons = pyr.PyRanges(exons)
    # pyranges of all protein-coding and lncRNA genes in GTF
    gtf_df = gtf_df[gtf_df.Feature == 'gene']
    # join two pyranges - if last exon overlaps with gene annotation then coordinates are joined
    # only overlaps on the same strand are kept
    joined = exons.join(gtf_df, strandedness="same", how=None)

    # print(joined.columns)
    # print(joined)
    joined = joined.as_df()
    # print(joined[['gene_id', 'gene_id_b']])

    # looking to filter out cases where exon gene_id differs to gene_id_b from gtf_df
    # i.e. exon overlaps with coordinates of a different gene
    joined = joined.groupby('transcript_id').filter(
        lambda x: (x['gene_id'] == x['gene_id_b']).all())

    # print(joined)
    transcript_id_list = joined.transcript_id.to_list()
    print(len(transcript_id_list))
    return transcript_id_list


non_overlapping_genes = get_non_overlapping_genes(gtf_df=gtf_pyranges)

print("number of last exons that do not overlap with different gene on the same strand is %s" %
      (len(non_overlapping_genes)))


def get_multi_polya_transcripts(gtf_df=None, subset_list=None, polya_bed_path=None, polya_version=2):

    # print(gtf_df.columns)
    gtf_df = gtf_df[gtf_df.transcript_id.isin(subset_list)]

    if polya_version == 1:
        polya_bed = pyr.readers.read_bed(f=polya_bed_path)

    elif polya_version == 2:
        polya_bed = pyr.readers.read_bed(f=polya_bed_path, as_df=True)
        # add chr prefix to Chromosome column (overlap won't work without same chromosome names)
        polya_bed['Chromosome'] = 'chr' + polya_bed['Chromosome'].astype(str)
        polya_bed = pyr.PyRanges(polya_bed)

    # print(polya_bed)
    gtf_last_exons = get_last_exons(gtf_df)
    # print(gtf_last_exons[['transcript_id', 'gene_id']])

    # 4. count_overlaps with BED file of polyA_sites
    gtf_last_exons = gtf_last_exons.count_overlaps(
        polya_bed, strandedness="same", keep_nonoverlapping=True)

    # print(gtf_last_exons[['transcript_id', 'gene_id', 'NumberOverlaps']])
    # print(gtf_last_exons[["transcript_id", "NumberOverlaps"]])
    # print(gtf_last_exons.as_df().groupby('gene_id').get_group(
    #    'ENSMUSG00000090394.8'))

    '''
    def get_best_supported_transcript(x):
        x['n_na'] = x.groupby('gene_id').apply(
            lambda x: x.isnull().sum(axis=1)).reset_index(drop=True)
        idx = x.groupby('gene_id')['n_na'].transform(min) == x['n_na']
        return x[idx]
    '''
    # Trying to minimise overlapping transcripts for each gene
    # Filter for 'best annotated transcript' in line with 'transcript_support_level' filter previously
    # 'best annotated' = fewest 'NaNs' for each transcript
    # if same number both are kept

    # gtf_last_exons = get_best_supported_transcript(gtf_last_exons.as_df())

    # n_exons_per_gene = gtf_last_exons.groupby('gene_id').size()
    # print(n_exons_per_gene[n_exons_per_gene > 1])
    # grouped = gtf_last_exons.as_df().groupby('gene_id').apply(lambda x: x.isnull().sum(axis=1))
    # grouped = grouped.idxmin()
    # print(grouped)

    # 5. Get list of transcript ids with at least two overlapping polyA_sites in terminal exon
    trs = gtf_last_exons[gtf_last_exons.NumberOverlaps >= 2]
    # print(trs)
    trs = trs.as_df()
    tr_list = list(set(trs.transcript_id.to_list()))
    # print(len(tr_list))
    return tr_list


multi_overlap_transcripts = get_multi_polya_transcripts(gtf_df=gtf_pyranges,
                                                        subset_list=non_overlapping_genes, polya_bed_path=polya_clusters, polya_version=atlas_version)

print("number of transcripts with multiple overlapping polyA_sites is %s" %
      (len(multi_overlap_transcripts)))


def write_overlapping_gtf(gtf_df=None, subset_list=None, outfile=None):
    '''
    Subsets gtf for transcripts containing at least two overlapping polyA_sites
    write to gtf
    '''
    gtf_df = gtf_df[gtf_df.transcript_id.isin(subset_list)]
    gtf_df.to_gtf(path=outfile,)


write_overlapping_gtf(gtf_df=gtf_pyranges, subset_list=multi_overlap_transcripts, outfile=out)
print("gtf containing transcipts with multiple overlapping polyA_sites in last exon written to %s" % (out))


'''
if __name__ == '__main__':
    gtf = sys.argv[1]
    polya_clusters = sys.argv[2]
    atlas_version = int(sys.argv[3])
    out = sys.argv[4]

    gtf_pyranges = pyr.readers.read_gtf(f=gtf)

    # list of transcript ids that do not overlap with other genes on the same strand
    non_overlapping_genes = get_non_overlapping_genes(gtf_df=gtf_pyranges)

    print("number of last exons that do not overlap with different gene on the same strand is %s" %
          (len(non_overlapping_genes)))

    # list of transcript ids with at least two overlapping polyA_sites in last exon
    multi_overlap_transcripts = get_multi_polya_transcripts(
        gtf_df=gtf_pyranges, subset_list=non_overlapping_genes, polya_bed_path=polya_clusters, polya_version=atlas_version)

    print("number of transcripts with multiple overlapping polyA_sites is %s" %
          (len(multi_overlap_transcripts)))

    write_overlapping_gtf(gtf_df=gtf_pyranges, subset_list=multi_overlap_transcripts, outfile=out)
    print("gtf containing transcipts with multiple overlapping polyA_sites in last exon written to %s" % (out))
'''
