# These are flags you must include - Two memory and one runtime.
# Runtime is either seconds or hours:min:sec

#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=18:00:00 

#These are optional flags but you probably want them in all jobs

#$ -S /bin/bash
#$ -j y
#$ -N paqr_annotations
#$ -cwd
#$ -pe smp 2
#$ -R y

#The code you want to run now goes here.

hostname
date

if [ "$1" != "" ]; then
    RUN_NAME=$1
else
    RUN_NAME=$""
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p $FOLDER
cp config.yaml $FOLDER/${RUN_NAME}_config.yaml

snakemake -s get_paqr_annotations.Snakefile \
--use-conda \
--jobscript cluster_qsub.sh \
--cluster-config cluster.yaml \
--cluster-sync "qsub -R y -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -o $FOLDER" \
-j 50 \
--rerun-incomplete \
--latency-wait 45 

date
