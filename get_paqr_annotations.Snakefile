import os
import yaml
configfile: 'config.yaml'

# Ensure output dir created prior to initiating pipeline
os.system("mkdir -p {0}".format(config['output_dir']))


rule all:
    input:
        os.path.join(config['output_dir'],"config.yaml"),
        #os.path.join(config['output_dir'], config['transcripts_name']),
        #os.path.join(config['output_dir'], config['clusters_name'])


rule get_non_overlap_multi_polyA_genes:
    input:
        gtf = config['gencode_gtf'],
        polyA_bed = config['polyA_bed']

    output:
        tr_gtf = os.path.join(config['output_dir'],
                              config['output_subdir'], config['transcripts_gtf_name'])
    params:
        script = os.path.join(config['scripts_dir'], "get_compliant_genes.py"),
        atlas_version = config['polyA_atlas_version'],
        rm_na_tsl = config['remove_na_transcript_support_level'],
        filter_best_isoforms = config['gene_best_supported_transcripts'],
        minimum_tsl = config['minimum_transcript_support_level'],
        strip_version_number = config['strip_version_number']

    conda: "paqr_annotations_env.yaml"

    log:
        os.path.join(config['output_dir'], "get_non_overlap_multi_polyA_genes.log")

    shell:
        '''
        python {params.script} \
        {input.gtf} \
        {input.polyA_bed} \
        {params.atlas_version} \
        {params.rm_na_tsl} \
        {params.filter_best_isoforms} \
        {params.minimum_tsl} \
        {params.strip_version_number} \
        {output.tr_gtf} \
        > {log}
        '''

rule get_transcripts_BED12:
    input:
        tr_gtf = os.path.join(config['output_dir'],
                              config['output_subdir'], config['transcripts_gtf_name'])
    output:
        tr_bed = os.path.join(config['output_dir'], config['transcripts_name'])

    params:
        tmp_genePred = os.path.join(config['output_dir'], config['output_subdir'], "temp.genePred")

    conda: "paqr_annotations_env.yaml"

    shell:
        '''
        gtfToGenePred {input.tr_gtf} {params.tmp_genePred}
        genePredToBed {params.tmp_genePred} {output.tr_bed}
        '''

rule get_clusters_BED:
    input:
        tr_gtf = os.path.join(config['output_dir'],
                              config['output_subdir'], config['transcripts_gtf_name'])
    output:
        clusters_bed = os.path.join(config['output_dir'], config['clusters_name']),

    params:
        script = os.path.join(config['scripts_dir'], "get_paqr_polya_bed.py"),
        polyA_bed = config['polyA_bed'],
        atlas_version = config['polyA_atlas_version']

    conda: "paqr_annotations_env.yaml"

    shell:
        '''
        python {params.script} \
        {input.tr_gtf} \
        {params.polyA_bed} \
        {params.atlas_version} \
        {output.clusters_bed}
        '''

rule copy_config:
    input:
        trs_bed = os.path.join(config['output_dir'], config['transcripts_name']),
        clusters_bed = os.path.join(config['output_dir'], config['clusters_name'])

    output:
        os.path.join(config['output_dir'],"config.yaml")

    shell:
        '''
        cp config.yaml {output}
        '''	
