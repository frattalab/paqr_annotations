import os
import yaml
configfile: 'config.yaml'

# Ensure output dir created prior to initiating pipeline
os.system("mkdir -p {0}".format(config['output_dir']))


rule all:
    input:
        os.path.join(config['output_dir'], config['transcripts_name']),
        os.path.join(config['output_dir'], config['clusters_name'])


rule get_non_overlap_multi_polyA_genes:
    input:
        gtf = config['gencode_gtf'],
        polyA_bed = config['polyA_bed']

    output:
        tr_gtf = os.path.join(config['output_dir'],
                              config['output_subdir'], config['transcripts_gtf_name'])
    params:
        script = os.path.join(config['scripts_dir'], "get_compliant_genes.py"),

    conda: "paqr_annotations_env.yaml"

    shell:
        '''
        python {params.script} \
        {input.gtf} \
        {input.polyA_bed} \
        {output.tr_gtf}
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
        gtftoGenePred {input.tr_gtf} {params.tmp_genePred}
        genePredToBed {params.tmp_genePred} {output.tr_bed}
        '''

rule get_clusters_BED:
    input:
        tr_gtf = os.path.join(config['output_dir'],
                              config['output_subdir'], config['transcripts_gtf_name'])
    output:
        clusters_bed = os.path.join(config['output_dir'], config['clusters_name'])

    params:
        script = os.path.join(config['scripts_dir'], "get_paqr_polya_bed.py"),
        polyA_bed = config['polyA_bed']

    conda: "paqr_annotations_env.yaml"

    shell:
        '''
        python {params.script} \
        {input.tr_gtf} \
        {params.polyA_bed} \
        {output.clusters_bed}
        '''
