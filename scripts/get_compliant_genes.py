#!/usr/bin/env python3


'''
Script to get non-overlapping genes with at least two overlapping polyA sites
Input: GENCODE GTF file & BED file of polyA_sites from PolyASite database
Output: GTF file containing transcripts in which
-their terminal exon does not overlap with the coordinates of a DIFFERENT gene
-their terminal exon contains at least two overlapping poly(A) sites (coordinates from BED file)
1. Filter GTF for 'protein_coding' & 'lncRNA' transcripts only
2. Optional filtering for transcript support level
3. overlap terminal exons with coordinates of all 'gene' Features in GTF
4. Remove terminal exons that overlap with coordinates of a different gene
5. For each valid terminal exon - count number of overlaps with poly(A) site clusters in provided BED file
6. Output GTF of valid transcripts with at least two overlapping poly(A) sites in the terminal exon
'''

import sys
import pyranges as pyr
import pandas as pd


#gtf = "/home/sam/paqr_annotations/tests/gencode.chr1.vM25.annotation.gtf"
#polya_clusters = "/home/sam/paqr_annotations/atlas.clusters.2.0.GRCm38.96.bed"
#atlas_version = 2
#remove_na = "yes"
#gene_best_isoforms = "yes"
#minimum_tsl = 5
#strip_ver_no = "yes"
#out = "/home/sam/paqr_annotations/tests/new_join_non_overlapping.multi_polya.gencode.chr1.vM25.annotation.gtf"

#gtf_pyranges = pyr.readers.read_gtf(f=gtf)
#remove_na = yes_no2Boolean(remove_na)
#gene_best_isoforms = yes_no2Boolean(gene_best_isoforms)
#strip_ver_no = yes_no2Boolean(strip_ver_no)


def yes_no2Boolean(x):
    '''
    x = "yes" returns True
    x = "no" returns False
    '''
    if x.lower() == "yes":
        return True
    elif x.lower() == "no":
        return False


def remove_na_tsl(pyranges):
    '''
    remove entries with a transcript support level tag of 'NA' or 'NaN'
    '''
    pyranges_df = pyranges.as_df()
    # if NA cannot be converted to integer - given null value
    pyranges_df = pyranges_df[pd.to_numeric(
        pyranges_df['transcript_support_level'], errors='coerce').notnull()]

    return pyr.PyRanges(pyranges_df)


def filter_min_tsl(pyranges, min=5):
    '''
    remove transcripts with transcript_support_level value less than or equal to min
    min should be integer between 1 & 5
    '''
    if min not in [n for n in range(1, 6)]:
        raise ValueError("min must be integer between 1 & 5. Please edit config.yaml")
    else:
        if min == 5:
            return pyranges
        else:
            pyranges = pyranges[pd.to_numeric(pyranges.transcript_support_level) <= min]
            return pyranges


def select_best_tsl_isoforms(pyranges):
    '''
    gene_id has multiple isoforms with different last exons
    Only keeps isoforms with the highest transcript_support_level
    e.g. if gene has 4 isoforms, 3 have TSL of 1 and 1 of TSL 3
    isoform with TSL 3 is dropped (3 rows/isoforms returned)
    '''
    # only need to filter gene_id groups with > 1 transcript for last exon
    x = pyranges.as_df()
    idx = x.groupby('gene_id')['transcript_support_level'].transform(
        min) == x['transcript_support_level']
    x = x[idx]

    return pyr.PyRanges(x)


def get_last_exons(ranges_obj):
    '''
    Returns last exon for each transcript
    1. Select exons from Feature column
    2. Group by transcript_id
    3. Select last exon for each group (largest 'exon_number' value)
    '''
    df = ranges_obj.as_df()
    df = df[df['Feature'] == 'exon']

    # set exon_number dtype to int32 (integer) (selecting max value doesn't work as intended - if >= 10 exons selects 9th)
    df = df.astype({'exon_number': 'int32'})

    # for each transcript_id, select the largest exon number as entry to retain
    idx = df.groupby('transcript_id')['exon_number'].transform(
        max) == df['exon_number']
    df = df[idx]
    return pyr.PyRanges(df)

    # Select largest 'exon_number' row for each transcript_id - only solution I could get to work...
    # https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-count-in-groups-using-groupby
    # sort by exon number value - last exon = top of sort
    # drop_duplicates keeps top row for each transcript_id - removes all rows but last exon
    #df = df.sort_values('exon_number', ascending=False).drop_duplicates(['transcript_id'])

    return pyr.PyRanges(df)


def get_non_overlapping_genes(gtf_df=None, rm_na_tsl=True, gene_best_tsl=True, min_tsl=5):
    '''
    returns list of transcript_ids of non-overlapping, protein-coding or lncRNA genes
    non-overlapping = last exon of transcript doesn't overlap with coordinates of a DIFFERENT gene

    1. Get pyranges object containing terminal exons of all protein_coding & lncRNA transcripts
    2. Filter these transcripts for transcript_support_level (if set)
    3. Subset gtf_df to contain coordinates of 'gene' features ONLY
    4. Perform same-stranded join of terminal exons pyranges & subsetted gtf_df pyranfes
    (columns & row values of overlapping interval added to row )
    5a. If TE's gene_id is different to subsetted gtf_df gene_id ('gene_id_b')
    --> This terminal exon overlaps with coordinates of a DIFFERENT gene
    5b. Filter out these terminal exon rows from the joined dataframe
    6. Return list of transcript_ids from filtered joined dataframe
    '''

    # 1a. protein-coding & lncRNA entries only
    exons = gtf_df[gtf_df.gene_type.isin(['protein_coding', 'lncRNA'])]
    print("number of distinct genes after filtering for protein coding and lncRNA transcripts is %s" % (
        len(set(exons.gene_id.to_list()))))

    # 1b. PyRanges object storing last exons for each transcript
    exons = get_last_exons(ranges_obj=exons)
    print("terminal exons extracted for %s distinct genes" % (len(set(exons.gene_id.to_list()))))
    print("terminal exons extracted for %s distinct transcipts" %
          (len(set(exons.transcript_id.to_list()))))

    # 2. TSL filtering
    if rm_na_tsl is True:
        # can do other filters
        print("number of distinct transcripts before removing TSL:NA is %s" %
              (len(set(exons.transcript_id.to_list()))))
        exons = remove_na_tsl(exons)
        print("number of distinct transcripts after removing TSL:NA is %s" %
              (len(set(exons.transcript_id.to_list()))))
        print("(sanity check) number of distinct genes after removing TSL:NA isoforms is %s" % (
            len(set(exons.gene_id.to_list()))))

        # filter for minimum transcript support level (all transcripts)
        exons = filter_min_tsl(pyranges=exons, min=min_tsl)
        print("number of distinct transcripts after filtering for minimum TSL:%s is %s" %
              (min_tsl, len(set(exons.transcript_id.to_list()))))
        print("(sanity check) number of distinct genes after filtering for minimum TSL:%s is %s" % (min_tsl,
                                                                                                    len(set(exons.gene_id.to_list()))))

        if gene_best_tsl is True:
            print("number of transcripts before selecting gene best supported isoforms is %s" %
                  (len(set(exons.transcript_id.to_list()))))

            exons = select_best_tsl_isoforms(exons)

            print("number of transcripts after selecting gene best supported isoforms is %s" %
                  (len(set(exons.transcript_id.to_list()))))
            print("(sanity check) number of distinct genes after selecting gene best supported isoforms is %s" % (
                len(set(exons.gene_id.to_list()))))

        elif gene_best_tsl is False:
            print("no gene-level best supported isoform filtering performed. If this was unexpected, double-check config.yaml")
            pass

    elif rm_na_tsl == False:
        # cannot do other TSL filters (having problems with sorting with NA & int values)
        print("No transcript support level pre-filtering performed. If this was unexpected, double-check config.yaml")
        pass

    # Interest is finding exons that can be unambigously assigned to a SINGLE gene
    # join two pyranges - if last exon overlaps with gene annotation then coordinates are joined
    # only overlaps on the same strand are kept
    # Then remove entries where gene_id of last exon is different to overlapping gene_id (gene_id_b)

    # 3. GENE entries only of GTF
    gtf_df = gtf_df[gtf_df.Feature == 'gene']
    print("total number of 'gene' entries (to check for overlapping) in GTF is %s" %
          (len(set(gtf_df.gene_id.to_list()))))

    # 4. Same stranded join of two objects - how = None means only keep overlapping intervals
    # (expect exon's coordinates to overlap with coordinates of its corresponding gene)
    exons = exons.join(gtf_df, strandedness="same", how=None)
    exons = exons.as_df()

    # 5 looking to filter out cases where exon gene_id differs to gene_id_b from gtf_df
    exons = exons.groupby('transcript_id').filter(
        lambda x: (x['gene_id'] == x['gene_id_b']).all())

    # 6.
    transcript_id_list = exons.transcript_id.to_list()
    # print(len(transcript_id_list))
    print("number of terminal exons that do not overlap with different gene on the same strand is %s" %
          (len(set(transcript_id_list))))

    return transcript_id_list


#non_overlapping_genes = get_non_overlapping_genes(gtf_df=gtf_pyranges)

# print("number of last exons that do not overlap with different gene on the same strand is %s" %
#      (len(non_overlapping_genes)))


def get_multi_polya_transcripts(gtf_df=None, subset_list=None, polya_bed_path=None, polya_version=2):
    '''
    Returns list of transcript ids with >=2 polyA_sites from polyA_bed_path overlapping last exon
    '''

    # print(gtf_df.columns)
    gtf_df = gtf_df[gtf_df.transcript_id.isin(subset_list)]
    print("number of transcripts being checked for overlapping PolyASite poly(A) sites is %s" %
          (len(set(gtf_df.transcript_id.to_list()))))

    if polya_version == 1:
        print("Processing poly(A) site BED file as if it has PolyASite 1.0 formatting")
        polya_bed = pyr.readers.read_bed(f=polya_bed_path)

    elif polya_version == 2:
        print("Processing poly(A) site BED file as if it has PolyASite 2.0 formatting")
        polya_bed = pyr.readers.read_bed(f=polya_bed_path, as_df=True)
        # add chr prefix to Chromosome column (overlap won't work without same chromosome names)
        polya_bed['Chromosome'] = 'chr' + polya_bed['Chromosome'].astype(str)
        polya_bed = pyr.PyRanges(polya_bed)

    # print(polya_bed)
    gtf_last_exons = get_last_exons(gtf_df)
    print("terminal exons extracted for %s distinct transcipts" %
          (len(set(gtf_last_exons.transcript_id.to_list()))))

    # print(gtf_last_exons[['transcript_id', 'gene_id']])

    # 4. count_overlaps with BED file of polyA_sites (adds column called NumberOverlaps)
    gtf_last_exons = gtf_last_exons.count_overlaps(
        polya_bed, strandedness="same", keep_nonoverlapping=True)

    # def get_best_supported_transcript(x):
    #    '''
    #    # Trying to minimise overlapping transcripts for each gene
    #    # Filter for 'best annotated transcript' in line with 'transcript_support_level' filter previously
    #    # 'best annotated' = fewest 'NaNs' for each transcript
    #    # if same number both are kept
    #    '''
    #    x['n_na'] = x.groupby('gene_id').apply(
    #        lambda x: x.isnull().sum(axis=1)).reset_index(drop=True)
    #    idx = x.groupby('gene_id')['n_na'].transform(min) == x['n_na']
    #    return x[idx]

    # gtf_last_exons = get_best_supported_transcript(gtf_last_exons.as_df())

    # 5. Get list of transcript ids with at least two overlapping polyA_sites in terminal exon
    trs = gtf_last_exons[gtf_last_exons.NumberOverlaps >= 2]
    print("number of transcripts with at least two overlapping poly(A) sites in their terminal exons is %s" % (
        len(set(trs.transcript_id.to_list()))))
    print("number of genes with at least 1 transcript with multiple,terminal-exon-overlapping polyA sites is %s" %
          (len(set(trs.gene_id.to_list()))))
    # print(trs)
    trs = trs.as_df()
    tr_list = list(set(trs.transcript_id.to_list()))
    # print(len(tr_list))
    return tr_list


# multi_overlap_transcripts = get_multi_polya_transcripts(gtf_df=gtf_pyranges,
#                                                        subset_list=non_overlapping_genes, polya_bed_path=polya_clusters, polya_version=atlas_version)

# print("number of transcripts with multiple overlapping polyA_sites is %s" %
#      (len(multi_overlap_transcripts)))


def write_overlapping_gtf(gtf_df=None, subset_list=None, strip_ver_no=True, outfile=None):
    '''
    Subsets gtf for transcripts containing at least two overlapping polyA_sites
    write to gtf
    '''
    gtf_df = gtf_df[gtf_df.transcript_id.isin(subset_list)]

    if strip_ver_no == True:
        # ENSMUST00000082908.1 --> ENSMUST00000082908 (expand to two columns, select first column)
        gtf_df.transcript_id = gtf_df.transcript_id.str.split('.', expand=True)[0]
        gtf_df.gene_id = gtf_df.gene_id.str.split('.', expand=True)[0]

    elif strip_ver_no == False:
        pass

    gtf_df.to_gtf(path=outfile,)


#write_overlapping_gtf(gtf_df=gtf_pyranges, subset_list=multi_overlap_transcripts, outfile=out)
#print("gtf containing transcipts with multiple overlapping polyA_sites in last exon written to %s" % (out))


if __name__ == '__main__':
    gtf = sys.argv[1]
    polya_clusters = sys.argv[2]
    atlas_version = int(sys.argv[3])
    remove_na = sys.argv[4]
    gene_best_isoforms = sys.argv[5]
    minimum_tsl = int(sys.argv[6])
    strip_version_no = sys.argv[7]
    out = sys.argv[8]

    # convert CLI strings of 'yes' or 'no' to True or False
    remove_na = yes_no2Boolean(remove_na)
    gene_best_isoforms = yes_no2Boolean(gene_best_isoforms)
    strip_version_no = yes_no2Boolean(strip_version_no)

    gtf_pyranges = pyr.readers.read_gtf(f=gtf)
    print("Number of distinct genes in %s is %s" % (gtf, len(set(gtf_pyranges.gene_id.to_list()))))

    # list of transcript ids that do not overlap with other genes on the same strand
    non_overlapping_genes = get_non_overlapping_genes(
        gtf_df=gtf_pyranges, rm_na_tsl=remove_na, gene_best_tsl=gene_best_isoforms, min_tsl=minimum_tsl)

    # list of transcript ids with at least two overlapping polyA_sites in last exon
    multi_overlap_transcripts = get_multi_polya_transcripts(
        gtf_df=gtf_pyranges, subset_list=non_overlapping_genes, polya_bed_path=polya_clusters, polya_version=atlas_version)

    # print("number of transcripts with multiple overlapping polyA_sites is %s" %
    #      (len(multi_overlap_transcripts)))

    write_overlapping_gtf(gtf_df=gtf_pyranges, subset_list=multi_overlap_transcripts,
                          strip_ver_no=strip_version_no, outfile=out)
    print("gtf containing transcipts with multiple overlapping polyA_sites in last exon written to %s" % (out))
