# PAQR Annotations Workflow

Snakemake pipeline to generate PAQR-compliant transcript annotations using GENCODE gene annotations of choice & PolyASite database (1.0/2.0 release) poly(A) sites

## General Notes

The workflow is functional (and its output runs cleanly through PAQR (KAPAC not tested yet...)) but is still a work in progress (tidying up, multithreading support...). Any feedback on code, functionality, clarity of instructions & bug fixes would be much appreciated!

See [linked issue on PAQR Github for discussion on generating the new annotations](https://github.com/zavolanlab/PAQR_KAPAC/issues/15)

## Details
Per recommendations in [following issue](https://github.com/zavolanlab/PAQR_KAPAC/issues/15), and details specified in the [original PAQR/KAPAC paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1415-3), PAQR input annotation files must meet the following criteria:
 - transcripts belong to 'protein_coding' or 'lncRNA' genes
 - The terminal exon coordinates of a transcript of a given gene must not overlap with another gene (on the same strand)
 - Non-overlapping terminal exons must contain at least two overlapping poly(A) sites defined by the PolyASite database

This workflow satisfies the non-overlapping criteria by identifying the terminal exon for each transcript and checking whether its coordinates overlap with the *gene* (selected tag in 'Feature' column of GTF file) coordinates of a **different gene**.

I have also included additional, customisable filters for the **'transcript support level' (TSL)** of annotations. [See Ensembl website for descriptions of different TSL flags](https://m.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#tsl). These filters include:
 - Filter out transcripts with an 'NA' TSL flag
 - Minimum TSL threshold
 - Select the 'best supported isoforms' for each gene (see config.yaml file for a more verbose explanation)

Under the examples directory, I have provided 'transcript' and 'cluster' annotation files generated using the GENCODE mouse vM25 GTF file ('reference chromosomes only' GTF file, see website for more details) & mouse PolyASite 2.0 release (Human potentially on the way (no promises!)). I have successfully ran both steps of PAQR (GitHub as of May 2020) with these annotation files. All transcripts have a transcript_support_level of 1, but you are able to generate files with more relaxed thresholds using this workflow if you are not happy with this.  

As a side note, this workflow is 'backwards compatible' with the PolyASite v1.0 release (wanted to sanity check my approach). Details on how my workflow lines up with provided annotations will follow in due course...

## Installation & Dependencies
This pipeline makes use of the following packages & version numbers:
- Python 3.6
- pyranges 0.0.77
- ucsc-genepredtobed 377 (Bioconda recipe version - originates from UCSC Tools)
- ucsc-gtftogenepred 377 (Bioconda recipe version - originates from UCSC Tools)
- Snakemake 5.17.0

Simplest way to solve dependencies for these packages is to use Conda. If you do not have a conda installation on your machine, see the [Anaconda website for installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Note: if you have Snakemake available in your path, it is not essential to install the Conda environment first (pass `--use-conda` flag to Snakemake when calling the pipeline). If you are happy to not install first, you can skip to the *Configuring & Running* step.

Once you have a Conda installation, you can set up a functional conda environment using the provided *paqr_annotations_env.yaml* environment file. Type the following command into your Terminal:
`conda env create -f paqr_annotations_env.yaml`

You'll likely receive some prompts from conda - work through these until the installation is complete. You should then double check the installation has worked correctly.
1. An environment named 'paqr_annotations' appears in list of environments
`conda env list`

2. You can activate the environment with the following command:
`conda activate paqr_annotations`

## Configuring and Running

As input to the workflow, you will need to download a **GENCODE GTF** file for your species of choice, as well as the **BED** file of poly(A) sites from the [PolyASite database](https://www.polyasite.unibas.ch/). You then need to edit the **config.yaml** file provided with this repository to point to your downloaded files, and to modify other parameters as you wish (more details in config.yaml file - if unclear feel free to open up an issue!).

Once you are happy with your configurations, I recommend performing a dry run to double-check filepaths, valid options etc. To do this type the following command into your Terminal:
`snakemake -n -p -s get_paqr_annotations.Snakefile`

If everything looks in order, you can then run the workflow with the following command:
`snakemake --cores 2 -s get_paqr_annotations.Snakefile`

Snakemake requires you to set a max number of cores (*--cores <max_number>*). You can put what you like here, but note I haven't made the pyranges commands multithreaded yet (on the to-do list) so it won't really make a difference.

If you already had Snakemake available in your path (and did not install the conda environment yourself previously), Snakemake can do it for you if you pass the *--use-conda* flag:
`snakemake --use-conda --cores 2 -s get_paqr_annotations.Snakefile`


## Acknowledgements

Many thanks to [Ralf Schmidt](https://github.com/koljaLanger) for his helpful clarifications and advice.

## To dos
1. How to handle transcripts with exactly the same terminal exon (carry all in annotations file (as doing now) or pick one - how to decide?)
2. Writing to BED via tab-separated instead of write_bed - losing correct 'basing' of coordinates?
..* could change column names to expected BED names & write to file via write_bed...
..* ignore this - to_bed seems very funky - clusters coordinates BED already so should be fine...
3. Multithreading (work out how to do some of my pandas functions without converting to DFs and back)
..* pip install -U ray (add to env.yaml - package for Multithreading support)
4. Option to exclude overlaps with other genes regardless of strand
5. replace sys.argv with argparse (lazy...)
6. Add specific details on how new release compares to old polyAsite release, how well my annotations overlap with provided annotations, drop off with TSL filters etc.
