# COVIDAnalysisPipeline

COVID analysis pipeline for Single-end and Paired-end reads.
The pipeline uses two different approaches: 1. Aligning directly on microbial genome. 2. Removing human-mapped reads and then aligning specifically on SARS-CoV2 genome.
The pipeline also includes SARS-CoV2 variants anad preparing data for phylogeny analysis.

Pipeline includes:

1. FastQC : for quality filter of reads
2. Trimmomatic: Trimming reads
3. KRAKEN2: Aligning reads on microbial genomes and identifying reads mapped on SARS-CoV2 genome.
4. Filtering out Human Reads using HISAT2
5. Extracting unmapped reads from Human Alignment.
6. Mapping the unmapped reads on SARS CoV2 genome
7. Converting Sam to Bam
8. Variant calling using VARSCAN.
9. Asssembly of COV2 reads using MEGAHIT and SPAdes.
