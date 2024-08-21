## generate anonymous reads from real data

This README describes the steps to generate anonymous reads from real data.

## anoyomize the reads
gatk ReadAnonymizer -I data/APA81-T.merged.dedup.bqsr.apa-genes.bam -O tests/data/real/not_contaminated.bam -R reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
gatk ReadAnonymizer -I data/APA18-T.merged.dedup.bqsr.apa-genes.bam -O tests/data/real/contaminated.bam -R reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

### convert to fastq

picard SamToFastq -I tests/data/real/not_contaminated.bam -F tests/data/real/not_contaminated.R1.fastq -F2 tests/data/real/not_contaminated.R2.fastq --VALIDATION_STRINGENCY LENIENT
picard SamToFastq -I tests/data/real/contaminated.bam -F tests/data/real/contaminated.R1.fastq -F2 tests/data/real/contaminated.R2.fastq --VALIDATION_STRINGENCY LENIENT

### rename the reads in the fastq files
seqtk rename tests/data/real/not_contaminated.R1.fastq | gzip - > tests/data/real/not_contaminated.renamed.R1.fastq.gz
seqtk rename tests/data/real/not_contaminated.R2.fastq | gzip - > tests/data/real/not_contaminated.renamed.R2.fastq.gz
seqtk rename tests/data/real/contaminated.R1.fastq | gzip - > tests/data/real/contaminated.renamed.R1.fastq.gz
seqtk rename tests/data/real/contaminated.R2.fastq | gzip - > tests/data/real/contaminated.renamed.R2.fastq.gz

### make interleaved fastq files
seqtk mergepe tests/data/real/not_contaminated.renamed.R1.fastq.gz tests/data/real/not_contaminated.renamed.R2.fastq.gz | gzip - > tests/data/real/not_contaminated.renamed.interleaved.fastq.gz
seqtk mergepe tests/data/real/contaminated.renamed.R1.fastq.gz tests/data/real/contaminated.renamed.R2.fastq.gz | gzip - > tests/data/real/contaminated.renamed.interleaved.fastq.gz

### align the reads to the reference genome using bwa
bwa mem -t 8 reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna tests/data/real/not_contaminated.renamed.R1.fastq.gz tests/data/real/not_contaminated.renamed.R2.fastq.gz | samtools sort -@8 -o tests/data/real/not_contaminated.renamed.bam -
bwa mem -t 8 reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna tests/data/real/contaminated.renamed.R1.fastq.gz tests/data/real/contaminated.renamed.R2.fastq.gz | samtools sort -@8 -o tests/data/real/contaminated.renamed.bam -

### index the bam files
samtools index tests/data/real/not_contaminated.renamed.bam
samtools index tests/data/real/contaminated.renamed.bam

### remove the intermediate files
rm tests/data/real/not_contaminated.R1.fastq tests/data/real/not_contaminated.R2.fastq tests/data/real/contaminated.R1.fastq tests/data/real/contaminated.R2.fastq tests/data/real/contaminated.bam tests/data/real/not_contaminated.bam tests/data/real/contaminated.bai tests/data/real/not_contaminated.bai
