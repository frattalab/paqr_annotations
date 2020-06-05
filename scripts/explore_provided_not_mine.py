#!/usr/bin/env python3


'''
Noticed differences between my workflow (vM14 & PolyASite 1.0) & provided annotations (vM14 & PolyASite 1.0)
I find ~3.5-4k genes in provided annotations do not pass through my workflow
A very small proportion can be explained by genes that do not have the 'prot_coding' or 'lncRNA' tags (~270)
    - all of which do not overlap with another gene & do have at least two poly(A) sites (i.e. pass my workflow)
What about the other 3.5k?
'''

import pyranges as pr
import pandas as pd
import os
# going to modify get_non_overlapping_genes, get_last_exons & get_multi_polya_transcripts functions slightly
from get_compliant_genes import filter_min_tsl, select_best_tsl_isoforms, remove_na_tsl

gtf = "/home/sam/paqr_annotations/gencode.vM14.annotation.gtf"
polyA_bed = "/home/sam/paqr_annotations/mm10.PolyASite.v1.0.clusters.bed"
provided_not_mine_list = "/home/sam/paqr_annotations/analysis/gene_list_provided_not_in_mine.txt"
mine_not_provided_list = "/home/sam/paqr_annotations/analysis/gene_list_mine_not_in_provided.txt"
output_dir = "/home/sam/paqr_annotations/comparisons/"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

'''
# ignore Changes to both - set number of cores to 4

Changes to get_non_overlapping_genes
1. Option to output joined df of last exons & gene coordinates
2. Option to output filtered df where check if gene_id == gene_id_b
3. Accepts subset_list

Changes to get_multi_polya_transcripts
1. Option to perform count_overlaps in non-strand-specific manner (i.e. ignore strand)
2. Can set NumberOverlaps filter value (default = 2)
'''


def get_last_exons(ranges_obj, as_numeric=False):
    '''
    Returns last exon for each transcript
    1. Select exons from Feature column
    2. Group by transcript_id
    3. Select last exon for each group (largest 'exon_number' value)
    '''
    df = ranges_obj.as_df()
    df = df[df['Feature'] == 'exon']

    if as_numeric is False:
        # Select largest 'exon_number' row for each transcript_id - only solution I could get to work...
        # https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-count-in-groups-using-groupby
        # sort by exon number value - last exon = top of sort
        # drop_duplicates keeps top row for each transcript_id - removes all rows but last exon
        df = df.sort_values('exon_number', ascending=False).drop_duplicates(['transcript_id'])
        return pr.PyRanges(df)
    elif as_numeric is True:
        df = df.sort_values(pd.to_numeric(
            df['exon_number'], errors='coerce'), ascending=False).drop_duplicates(['transcript_id'])
        return pr.PyRanges(df)


def get_non_overlapping_genes(gtf_df=None, subset_list=['protein_coding', 'lncRNA'], rm_na_tsl=True, gene_best_tsl=True, min_tsl=5, output_joined_df=False, output_filtered_df=False, alt_last_function=False):
    '''
    returns list of transcript_ids of non-overlapping, protein-coding or lncRNA genes
    non-overlapping = last exon of transcript doesn't overlap with coordinates of a DIFFERENT gene

    Extra options
    output_joined_df - if True returns pyranges object of last exons joined with overlapping gene_coordinates
    output_filtered_df - if True returns pyranges object of joined df with last exons overlapping with different genes removed


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
    if subset_list[0] == 'protein_coding':
        # filtering for gene_type
        exons = gtf_df[gtf_df.gene_type.isin(subset_list)]

    else:
        # assume filtering for gene_id
        exons = gtf_df[gtf_df.gene_id.isin(subset_list)]

    print("number of distinct genes after filtering for protein coding and lncRNA transcripts is %s" % (
        len(set(exons.gene_id.to_list()))))

    # 1b. PyRanges object storing last exons for each transcript
    if alt_last_function is True:
        exons = alt_get_last_exons(pyranges=exons)

    elif alt_last_function is False:
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
    exons = exons.join(gtf_df, strandedness="same", how=None, nb_cpu=1)

    if output_joined_df is True:
        return exons

    exons = exons.as_df()

    # 5 looking to filter out cases where exon gene_id differs to gene_id_b from gtf_df
    exons = exons.groupby('transcript_id').filter(
        lambda x: (x['gene_id'] == x['gene_id_b']).all())

    if output_filtered_df is True:
        return pr.PyRanges(exons)
    # 6.
    transcript_id_list = exons.transcript_id.to_list()
    # print(len(transcript_id_list))
    print("number of terminal exons that do not overlap with different gene on the same strand is %s" %
          (len(set(transcript_id_list))))

    return transcript_id_list


def get_multi_polya_transcripts(gtf_df=None, subset_list=None, polya_bed_path=None, polya_version=2, strand_specific=True, min_NOverlaps=2, alt_last_function=False):
    '''
    Returns list of transcript ids with >=2 polyA_sites from polyA_bed_path overlapping last exon
    '''

    # print(gtf_df.columns)
    gtf_df = gtf_df[gtf_df.transcript_id.isin(subset_list)]
    print("number of transcripts being checked for overlapping PolyASite poly(A) sites is %s" %
          (len(set(gtf_df.transcript_id.to_list()))))

    if polya_version == 1:
        print("Processing poly(A) site BED file as if it has PolyASite 1.0 formatting")
        polya_bed = pr.readers.read_bed(f=polya_bed_path)

    elif polya_version == 2:
        print("Processing poly(A) site BED file as if it has PolyASite 2.0 formatting")
        polya_bed = pr.readers.read_bed(f=polya_bed_path, as_df=True)
        # add chr prefix to Chromosome column (overlap won't work without same chromosome names)
        polya_bed['Chromosome'] = 'chr' + polya_bed['Chromosome'].astype(str)
        polya_bed = pr.PyRanges(polya_bed)

    # print(polya_bed)
    if alt_last_function is True:
        gtf_last_exons = alt_get_last_exons(gtf_df)
    elif alt_last_function is False:
        gtf_last_exons = get_last_exons(gtf_df)

    print("terminal exons extracted for %s distinct transcipts" %
          (len(set(gtf_last_exons.transcript_id.to_list()))))

    # print(gtf_last_exons[['transcript_id', 'gene_id']])

    # 4. count_overlaps with BED file of polyA_sites (adds column called NumberOverlaps)
    if strand_specific is True:
        gtf_last_exons = gtf_last_exons.count_overlaps(
            polya_bed, strandedness="same", keep_nonoverlapping=True)

    elif strand_specific is False:
        # ignore strand information for overlap
        gtf_last_exons = gtf_last_exons.count_overlaps(
            polya_bed, strandedness=False, keep_nonoverlapping=True)

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
    trs = gtf_last_exons[gtf_last_exons.NumberOverlaps >= min_NOverlaps]
    print("number of transcripts with at least two overlapping poly(A) sites in their terminal exons is %s" % (
        len(set(trs.transcript_id.to_list()))))
    print("number of genes with at least 1 transcript with multiple,terminal-exon-overlapping polyA sites is %s" %
          (len(set(trs.gene_id.to_list()))))
    # print(trs)
    trs = trs.as_df()
    tr_list = list(set(trs.transcript_id.to_list()))
    # print(len(tr_list))
    return tr_list


def get_trs_for_gene_list(pyranges=None, gene_list=None):
    '''
    Subset pyranges for gene_ids in gene_list, return list of transcript_ids (all belonging to genes in gene_list)
    '''
    pyranges = pyranges[pyranges.Feature == 'transcript']
    # print(pyranges)
    pyranges = pyranges[pyranges.gene_id.isin(gene_list)]
    # print(pyranges)
    # print(pyranges.columns)

    trs = list(set(pyranges.transcript_id.to_list()))
    return trs


# ---------------------------
# Read in files
# ---------------------------

with open(provided_not_mine_list) as infile:
    provided_not_mine_genes = [line.rstrip() for line in infile]

with open(mine_not_provided_list) as infile:
    mine_not_provided_genes = [line.rstrip() for line in infile]

# print(provided_not_mine_genes)

print("number of genes present in my workflow's clusters file but not in provided clusters file is %s" %
      (len(mine_not_provided_genes)))

print("number of genes present in provided clusters file but not in my workflow's clusters file is %s" % (
    len(provided_not_mine_genes)))

gtf_pyranges = pr.readers.read_gtf(f=gtf,)
# print(gtf_pyranges)
# print(gtf_pyranges.columns)

# strip version numbers
gtf_pyranges.transcript_id = gtf_pyranges.transcript_id.str.split('.', expand=True)[0]
gtf_pyranges.gene_id = gtf_pyranges.gene_id.str.split('.', expand=True)[0]

gtf_pyranges_provided = gtf_pyranges[gtf_pyranges.gene_id.isin(provided_not_mine_genes)]

print(gtf_pyranges_provided[gtf_pyranges_provided.transcript_id == 'ENSMUST00000081551'].print(n=30))
print(gtf_pyranges_provided.dtypes)
# check get_last_exons_function - think it may be failing for transcripts with >= 10 exons


def alt_get_last_exons(pyranges):
    pyranges = pyranges.as_df()
    pyranges = pyranges[pyranges.Feature == 'exon']

    # set exon_number dtype to int32 (integer) (selecting max value doesn't work as intended - if >= 10 exons selects 9th)
    pyranges = pyranges.astype({'exon_number': 'int32'})

    # for each transcript_id, select the largest exon number as entry to retain
    idx = pyranges.groupby('transcript_id')['exon_number'].transform(
        max) == pyranges['exon_number']
    pyranges = pyranges[idx]
    return pr.PyRanges(pyranges)


existing_last_exons = get_last_exons(gtf_pyranges_provided)
new_f_last_exons = alt_get_last_exons(gtf_pyranges_provided)

print(existing_last_exons[['transcript_id', 'exon_number']])
print(new_f_last_exons[['transcript_id', 'exon_number']])

print(existing_last_exons[existing_last_exons.transcript_id == 'ENSMUST00000081551'].print())
print(new_f_last_exons[new_f_last_exons.transcript_id == 'ENSMUST00000081551'].print())

# print(existing_last_exons[existing_last_exons.transcript_id == 'ENSMUST00000081551'].print(n=30))

#trs_to_check = ['ENSMUST00000193812','ENSMUST00000082908','ENSMUST00000192857','ENSMUST00000161581','ENSMUST00000183034','ENSMUST00000097786','ENSMUST00000187582','ENSMUST00000191048']
#existing_to_print = existing_last_exons[existing_last_exons.transcript_id.isin(trs_to_check)]
#new_f_to_print = new_f_last_exons[new_f_last_exons.transcript_id.isin(trs_to_check)]

#print("existing function\n")
# print(existing_to_print[['transcript_id','exon_number']])
#print("new alternative function\n")
# print(new_f_to_print[['transcript_id','exon_number']])

# pandas to join the two dataframes by transcript_id


def join_2_pyranges_pandas_like(x, y):
    x = x.as_df()
    y = y.as_df()

    z = x.join(y.set_index('transcript_id'), on='transcript_id', rsuffix='_b')

    return pr.PyRanges(z)


existing_new_f_joined = join_2_pyranges_pandas_like(existing_last_exons, new_f_last_exons)

# get cases where exon numbers do not match up
existing_new_f_joined = existing_new_f_joined.subset(lambda df: pd.to_numeric(
    df.exon_number, errors='coerce') != pd.to_numeric(df.exon_number_b, errors='coerce'))

print(existing_new_f_joined[['transcript_id', 'exon_number', 'Start_b', 'End_b', 'Strand_b', 'exon_number_b']].print(n=20))


# PROBLEM IS BECAUSE GET LAST EXONS WAS SELECTING MATHCING VALUE BASED ON OBJECT DATATYPE NOT AN INTEGER DATATYPE
# if transcript had >= 10 exons - exon 9 was selected :(

# function check - does drop_duplicates approach work if I convert to numeric first?
#convert_numeric_existing_last_exons = get_last_exons(gtf_pyranges_provided,as_numeric=True)
#print("checking existing function if convert to numeric inside it\n")
# print(convert_numeric_existing_last_exons[['transcript_id','exon_number']])


# Run polyA overlap on provided_not_mine_transcripts using existing and alternative get_last_exons function
print("--------\nchecking for provided_not_mine last exon poly(A) overlap if use existing get_last_exons function\n")

provided_not_mine_transcripts = get_trs_for_gene_list(
    pyranges=gtf_pyranges, gene_list=provided_not_mine_genes)

get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, subset_list=provided_not_mine_transcripts, polya_bed_path=polyA_bed, polya_version=1)


print("--------\nchecking for provided_not_mine last exon poly(A) overlap if use alternative get_last_exons function\n")
get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, subset_list=provided_not_mine_transcripts, polya_bed_path=polyA_bed, polya_version=1, alt_last_function=True)

# check how many would not overlap with other genes on the same strand
print("--------\ncheck how many with updated last exons function would be valid with my worklow (non-overlapping with a different gene)\n")

provided_not_mine_valid_alt_lastExons_f_transcripts = get_non_overlapping_genes(
    gtf_df=gtf_pyranges, subset_list=provided_not_mine_genes, rm_na_tsl=False, alt_last_function=True)

provided_not_mine_valid_alt_lastExons_f_multi_polya_transcripts = get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, polya_bed_path=polyA_bed, polya_version=1, subset_list=provided_not_mine_valid_alt_lastExons_f_transcripts, alt_last_function=True)

'''

# 1. DO THE GENES IN PROVIDED ANNOTATIONS BUT NOT MY WORKFLOW HAVE 'protein_coding' or 'lncRNA' gene_type tags?
print("\nAssessing whether provided_not_mine genes lack a 'protein_coding' or 'lncRNA' gene_type tag\n")

gtf_pyranges_provided = gtf_pyranges[gtf_pyranges.gene_id.isin(provided_not_mine_genes)]
not_prot_lncRNA = gtf_pyranges_provided[~gtf_pyranges_provided.gene_type.isin(
    ['protein_coding', 'lncRNA'])]

# print(not_prot_lncRNA)
# print(not_prot_lncRNA.columns)

print("number of genes without 'protein_coding' or 'lncRNA' gene_type tag is %s" %
      (len(set(not_prot_lncRNA.gene_id.to_list()))))

# Would genes without protein_coding or lncRNA gene_type tags pass my workflow if considered?
print("\nCHECKING GENES WITHOUT 'protein_coding' OR 'lncRNA' GENE_TYPE TAGS ONLY (do they pass my workflow?)\n(no TSL filtering considered)\n")

not_prot_lncRNA_valid_transcripts = get_non_overlapping_genes(
    gtf_df=gtf_pyranges, subset_list=not_prot_lncRNA.gene_id.to_list(), rm_na_tsl=False)

not_prot_lnc_multi_polyA_transcripts = get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, polya_bed_path=polyA_bed, polya_version=1, subset_list=not_prot_lncRNA_valid_transcripts)

print(len(not_prot_lnc_multi_polyA_transcripts))

with open(os.path.join(output_dir, "no_prot_coding_lncRNA_tag_non_overlapping_multi_polyA_transcripts.txt"), "w") as outfile:
    for transcript in not_prot_lnc_multi_polyA_transcripts:
        outfile.write("%s\n" % (transcript))


# 2. CHECK IF PROVIDED BUT NOT MINE GENES HAVE >2 OVERLAPPING POLY(A) SITES IN LAST EXON
# Skip non_overlapping_genes function, get list of transcript_ids for provided_not_mine_list
print("\n CHECKING WHETHER PROVIDED_NOT_MINE GENES HAVE MULTIPLE OVERLAPPING POLYA SITES IN LAST EXON\n")

provided_not_mine_transcripts = get_trs_for_gene_list(
    pyranges=gtf_pyranges, gene_list=provided_not_mine_genes)

print("number of transcripts for genes (n genes =%s) in %s is %s" %
      (len(provided_not_mine_genes), provided_not_mine_list, len(provided_not_mine_transcripts)))

# print(provided_not_mine_transcripts)
# run get_multi_polya_transcripts using provided_not_mine_transcripts as subset_list
multi_polya_provided_not_mine_transcripts = get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, subset_list=provided_not_mine_transcripts, polya_bed_path=polyA_bed, polya_version=1)

# print("number of transcripts in which last exon has multiple poly(A) sites is %s" %
#      (len(multi_polya_provided_not_mine_transcripts)))

# How many of these multi-poly(A) genes are in the non -tagged protein_coding or lncRNA?
print("\nchecking overlap between multi_polyA provided_not_mine transcripts list and non-tagged (no protein_coding or lncRNA tag) multi_polyA transcripts list")
print("number of multi_polya_provided_not_mine_transcripts that belong to transcripts without gene_type tags is %s" % (len(
    [transcript for transcript in multi_polya_provided_not_mine_transcripts if transcript in not_prot_lnc_multi_polyA_transcripts])))
print("number of multi_polya_provided_not_mine_transcripts that belong to transcripts with valid gene_type tags is %s" % (len(
    [transcript for transcript in multi_polya_provided_not_mine_transcripts if transcript not in not_prot_lnc_multi_polyA_transcripts])))

with open(os.path.join(output_dir, "provided_not_mine_prot_lncRNA_tagged_multi_polyA_transcripts.txt"), "w") as outfile:
    for transcript in multi_polya_provided_not_mine_transcripts:
        if transcript not in not_prot_lnc_multi_polyA_transcripts:
            outfile.write("%s\n" % (transcript))

# 2b. CHECK IF PROVIDED BUT NOT MINE GENES HAVE >2 OVERLAPPING POLY(A) SITES IN LAST EXON IF OVERLAP IGNORES STRAND
print("\n------------------\nCHECKING IF PROVIDED_NOT_MINE GENES HAVE >2 OVERLAPPING POLY(A) SITES IN LAST EXON IF OVERLAP IGNORES STRAND\n")

ignore_strand_multi_polya_provided_not_mine_transcripts = get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, subset_list=provided_not_mine_transcripts, strand_specific=False, polya_bed_path=polyA_bed, polya_version=1)

# provided_not_mine transcripts that have multiple poly(A) sites if you perform overlap whilst ignoring strand (only - if in strand-specific set they are ignored)
with open(os.path.join(output_dir, "provided_not_mine_prot_lncRNA_tagged_ignore_strand_multi_polyA_transcripts.txt"), "w") as outfile:
    for transcript in ignore_strand_multi_polya_provided_not_mine_transcripts:
        if transcript not in multi_polya_provided_not_mine_transcripts:
            outfile.write("%s\n" % (transcript))


# print("number of transcripts in which last exon has multiple poly(A) sites (ignoring strand for overlap) is %s" %
#      (len(ignore_strand_multi_polya_provided_not_mine_transcripts)))


# 3. CHECKING WHETHER NUMBER OF OVERLAPS >= 2 IS AFFECTING WHETHER PROVIDED_NOT_MINE GENES ARE MISSED WITH MY WORKFLOW
# i.e. I require 2 overlapping polyAsite sites, but provided considers annotation 3'end as a poly(A) site
print("\nCHECKING WHETHER NUMBER OF OVERLAPS >= 2 IS AFFECTING WHETHER PROVIDED_NOT_MINE GENES ARE MISSED WITH MY WORKFLOW\ni.e. I require 2 overlapping polyAsite sites, but provided considers annotation 3'end as a poly(A) site")
print("(Considering overlap as strand-specific again)\n")

nOverlaps_1_multi_polya_provided_not_mine_transcripts = get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, subset_list=provided_not_mine_transcripts, polya_bed_path=polyA_bed, polya_version=1, min_NOverlaps=1)


with open(os.path.join(output_dir, "provided_not_mine_prot_lncRNA_tagged_nOverlaps_1_multi_polyA_transcripts.txt"), "w") as outfile:
    for transcript in nOverlaps_1_multi_polya_provided_not_mine_transcripts:
        # not in list where provided_not_mine have >=2 overlaps with polyAsite sites
        if transcript not in multi_polya_provided_not_mine_transcripts:
            outfile.write("%s\n" % (transcript))

# Repeat the nOverlaps = 1 but do overlap when ignoring strand
print("\n---------------\nREPEAT ONLY 1 OVERLAPPING POLYASITE REQUIRED BUT PERFORM OVERLAP IGNORING STRAND INFORMATION\n")
ignore_strand_nOverlaps_1_multi_polya_provided_not_mine_transcripts = get_multi_polya_transcripts(
    gtf_df=gtf_pyranges, subset_list=provided_not_mine_transcripts, strand_specific=False, polya_bed_path=polyA_bed, polya_version=1, min_NOverlaps=1)

# How many genes & transcripts NEVER HAVE overlapping poly(A) sites regardless of how relaxed parameters are?
print("\n---------------\nCHECKING HOW MANY GENES AND TRANSCRIPTS NEVER HAVE OVERLAPPING POLYA SITES REGARDLESS OF LOOSENESS OF OVERLAPPING PARAMETERS\n")

# file containining transcripts in which I could not find any overlapping poly(A) sites (REGARDLESS OF OPTIONS/PARAMETERS used)
# ignore_strand_nOverlaps_1_multi_polya_provided_not_mine_transcripts captures all possibilities in other transcript lists
# (overlap regardless of strand, includes non-tagged, only one PolyASite poly(A) site)
provided_not_mine_never_any_multi_polyA_transcripts = [
    transcript for transcript in provided_not_mine_transcripts if transcript not in ignore_strand_nOverlaps_1_multi_polya_provided_not_mine_transcripts]

#print(len(provided_not_mine_never_any_multi_polyA_transcripts))
print("number of transcripts belonging to genes in provided_not_mine_gene_list that do not have overlapping poly(A) sites regardless of above searches is %s" %
      (len(set(provided_not_mine_never_any_multi_polyA_transcripts))))

print("(sanity check) - number of transcripts in 'gtf_pyranges_provided' is %s" % (len(set(gtf_pyranges_provided.transcript_id.to_list()))))

#Get gene_ids corresponding to transcripts in ignore_strand_nOverlaps_1_multi_polya_provided_not_mine_transcripts
provided_not_mine_any_multi_polyA_genes = gtf_pyranges_provided[gtf_pyranges_provided.transcript_id.isin(ignore_strand_nOverlaps_1_multi_polya_provided_not_mine_transcripts)].gene_id.to_list()

print("number of genes corresponding to transcripts with multiple poly(A) sites by any parameters is %s" % (len(set(provided_not_mine_any_multi_polyA_genes))))

#get list of gene_ids in which never any overlap (i.e. not in provided_not_mine_never_any_multi_polyA_genes)
provided_not_mine_never_any_multi_polyA_genes = list(set(gtf_pyranges_provided[~gtf_pyranges_provided.gene_id.isin(provided_not_mine_any_multi_polyA_genes)].gene_id.to_list()))

print("(sanity check) number of genes in which none of its transcripts ever have multi_polyA transcripts (any parameters) is %s " % (len(set(provided_not_mine_never_any_multi_polyA_genes))))

#write never any multiple polyA_genes to file
with open(os.path.join(output_dir,"provided_not_mine_never_any_multi_polyA_tr_gene_ids.txt"),"w") as outfile:
    for gene in provided_not_mine_never_any_multi_polyA_genes:
        outfile.write("%s\n" % gene)

#write never any multiple polyA transcripts to file
with open(os.path.join(output_dir,"provided_not_mine_never_any_multi_polyA_transcript_ids.txt"),"w") as outfile:
    for transcript in provided_not_mine_never_any_multi_polyA_transcripts:
        outfile.write("%s\n" % transcript)


#Write GTF containing all never_any_multi_polyA_transcripts to file
gtf_pyranges_provided = gtf_pyranges_provided[gtf_pyranges_provided.transcript_id.isin(provided_not_mine_never_any_multi_polyA_transcripts)]
gtf_pyranges_provided.to_gtf(path=os.path.join(output_dir,"gencode.mm10.vM14.atlas1.provided.not.mine.never.any.multi.polyA.transcripts.gtf"))

#write GTF containing terminal exons of all_never_any_multi_polyA_transcripts to file
gtf_pyranges_provided = get_last_exons(ranges_obj=gtf_pyranges_provided)
print("(sanity check) - number of terminal exons of never_any_multi_polyA_transcripts is %s (should = number of transcripts)" % (len(set(gtf_pyranges_provided.transcript_id.to_list()))))

gtf_pyranges_provided.to_gtf(path=os.path.join(output_dir,"gencode.mm10.vM14.atlas1.provided.not.mine.never.any.multi.polyA.Texons.gtf"))

#test = gtf_pyranges_provided[gtf_pyranges_provided.Feature == 'transcript']
#test = test[~test.transcript_id.isin(
#    ignore_strand_nOverlaps_1_multi_polya_provided_not_mine_transcripts)]

#print("n genes in test is %s" % len(set(test.gene_id.to_list())))
#provided_not_mine_never_any_multi_polyA_genes = gtf_pyranges_provided[gtf_pyranges_provided.transcript_id.isin(
#    provided_not_mine_never_any_multi_polyA_transcripts)].gene_id.to_list()

# print(provided_not_mine_never_any_multi_polyA_transcripts)
'''
