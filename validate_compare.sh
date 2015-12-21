FAST5 = "/test/input_files/fast5/loman10.txt" # path to fast5 files

## Basecalling

## Make random reference

## Convert metrichor basecalling to fasta
cat $FAST5 | xargs poretools fasta > metrichor.fa

## Make stats for all three fasta files
./validate.py called.fa > called.stats
./validate.py metrichor.fa > metrichor.stats
./validate.py random.fa > random.stats