#!/usr/bin/env bash
export HDF5_DISABLE_VERSION_CHECK=2
FAST5_LIST="test/input_files/fast5/loman100.txt" # path to fast5 files
FASTA_OUTPUT="test/output_files/fasta"
STATS_OUTPUT="test/output_files/stats"
REF="/storageNGS/ngs3/projects/other/Nanopore/PublicData/LomanLab_MAP-006/ecoli_mg1655.fa"
GRAPHMAP=/home/ibis/gregor.sturm/nanopore/tools/graphmap/graphmap
mkdir -p $FASTA_OUTPUT
mkdir -p $STATS_OUTPUT


## Basecalling
./basecall.py --template test/output_files/models/model6-1.template_median68pA.pickle \
    --complement test/output_files/models/model6-1.complement_median68pA_pop2.pickle \
    --filelist $FAST5_LIST --output $FASTA_OUTPUT/called.fa

## Convert metrichor basecalling to fasta
cat $FAST5_LIST | xargs poretools fasta > $FASTA_OUTPUT/metrichor.fa 

## Make random reference
./validate.py randomize -f $FASTA_OUTPUT/metrichor.fa > $FASTA_OUTPUT/random.fa

## split fasta files in template and complement
for strand in template complement; do
    for dataset in metrichor called random; do
        ./validate.py filter -k $strand -f $FASTA_OUTPUT/$dataset.fa > $FASTA_OUTPUT/$dataset.$strand.fa
    done
done

## Make stats for all three fasta files
for strand in template complement; do
    for dataset in metrichor called random; do
         ./validate.py align --fasta $FASTA_OUTPUT/$dataset.$strand.fa --reference $REF \
            --sam_basename $FASTA_OUTPUT/$dataset.$strand --graphmap $GRAPHMAP
         ./validate.py stats --bam $FASTA_OUTPUT/$dataset.$strand.bam \
            --reference $REF > $STATS_OUTPUT/$dataset.$strand.stats
    done
done
