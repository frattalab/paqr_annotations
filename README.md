#PAQR Annotations Pipeline

Snakemake pipeline to generate PAQR-compliant transcript annotations with updated reference transcript set & polyAsite 2.0 release poly(A) sites

Starts with a Gencode transcript annotation GTF
Per recommendations in issue X:
-Genes that overlap on the same strand are excluded (default - can also exclude genes overlapping regardless of strand)
-Terminal exon of gene must contain at least two overlapping poly(A) sites (from polyAsite BED file


##To dos
1. How to handle transcripts with exactly the same terminal exon (carry all in annotations file or pick one - how to decide?)
2. Writing to BED via tab-separated instead of write_bed - losing correct 'basing' of coordinates?
  - could change column names to expected BED names & write to file via write_bed...
3. Multithreading (work out how to do some of my pandas functions without converting to DFs and back)
