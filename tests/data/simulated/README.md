# commands for sequence data simulation

## download reads
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752113/real_wes_reads_probes.fastq.gz

## download the reqseq profile
wget https://github.com/schmeing/ReSeq-profiles/raw/master/profiles/Hs-Nova-TruSeq.reseq

## Randomly replace all Ns in the genome being sequenced using the following command:
reseq replaceN -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -R GCA_000001405.15_GRCh38_no_alt_analysis_set_noNs.fasta

## Generate systematic errors for the genome
reseq illuminaPE -r GCA_000001405.15_GRCh38_no_alt_analysis_set_noNs.fasta -s Hs-Nova-TruSeq.reseq \
--stopAfterEstimation --writeSysError GCA_000001405.15_GRCh38_no_alt_analysis_set_noNs_syserrors.fq
