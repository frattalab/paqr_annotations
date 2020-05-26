```PAQR Annotations Pipeline

Snakemake pipeline to generate PAQR-compliant transcript annotations with updated reference transcript set & polyAsite 2.0 release poly(A) sites

Starts with a Gencode transcript annotation GTF
Per recommendations in issue X:
-Genes that overlap on the same strand are excluded (default - can also exclude genes overlapping regardless of strand)
-Terminal exon of gene must contain at least two overlapping poly(A) sites (from polyAsite BED file


``To dos
1. Option to strip version numbers from transcript & gene_ids
2. Options to include/exclude TSL filtering (& impose minimum TSL filter for overlapping transcripts)
3. Multithreading (work out how to do some of my pandas functions without converting to DFs and back)
