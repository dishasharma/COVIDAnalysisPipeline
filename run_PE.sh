#!/bin/bash 
start=`date +%s`
usage()
{
  echo "Usage: ./run_PE.sh -d <Directory of fastq files> -t <Directory where tools are installed> -r <Human_Reference_Genome.fa> -c <Covid_refernce_genome.fa> -z <suffix of samples (_R1.fastq.gz/_R2.fastq.gz or _R1.fastq/_R2.fastq) -n <output folder name> -g <GISAID fasta path> -m <mao file> -h <help>"
  exit 2
}

while getopts d:t:r:c:z:n:g:m:h: option 
do 
 case "${option}" 
 in 
 d) DIRECTORY=${OPTARG};;
 t) TOOLS_DIRECTORY=${OPTARG};;
 r) HUMAN_REFERENCE=${OPTARG};;
 c) COVID_REFERENCE=${OPTARG};;
 z) FASTQ_ZIP=${OPTARG};;
 n) OUT=${OPTARG};;
 g) GISAID_PATH=${OPTARG};;
 m) MAO_PATH=${OPTARG};;
 h|?) usage ;; esac
done

mkdir Fastqc_output_$OUT
mkdir Trimmed_output_$OUT
mkdir Kraken_output_$OUT
mkdir Krona_output_$OUT
mkdir Hisat2_human_output_$OUT
mkdir Hisat2_covid_output_$OUT
mkdir Variant_calling_output_$OUT
mkdir DenovoAssembly_output_$OUT
mkdir Phylogeny_output_$OUT

FileList="$(ls $DIRECTORY/*$FASTQ_ZIP | awk '{ print $1}' | awk -F'/' '{ print $2}' | grep "_R1" | awk -F'_R1' '{ print $1}')"

##picard CollectMultipleMetrics I=NCDC2105.bam O=alignmentsummary R=../genome/GRCh38.p13.genome.fa

for i in $FileList
do
	echo "Starting analysis for "$i
	echo -e $DIRECTORY"\t"$i
	FASTQC_COMMAND_R1=`echo -e $DIRECTORY"\t"$i"\t"$TOOLS_DIRECTORY"\t"$FASTQ_ZIP"\t"$OUT | awk '{ print $3"/FastQC/fastqc "$1"/"$2"_R1_"$4" -o Fastqc_output_"$5}'`
        echo $FASTQC_COMMAND_R1
	#eval "$FASTQC_COMMAND_R1"
        FASTQC_COMMAND_R2=`echo -e $DIRECTORY"\t"$i"\t"$TOOLS_DIRECTORY"\t"$FASTQ_ZIP"\t"$OUT | awk '{ print $3"/FastQC/fastqc "$1"/"$2"_R2_"$4" -o Fastqc_output_"$5}'`
        echo $FASTQC_COMMAND_R2
        #eval "$FASTQC_COMMAND_R2"
	#TRIM_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$TOOLS_DIRECTORY | awk '{ print "TrimmomaticSE "$1"/"$2"*.fastq.gz Trimmed_output/"$2"_trimmed.fastq ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"}'`
        TRIM_COMMAND_PE=`echo -e $DIRECTORY"\t"$i"\t"$TOOLS_DIRECTORY"\t"$FASTQ_ZIP"\t"$OUT | awk '{ print "java -jar "$3"/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 "$1"/"$2"_R1_"$4" "$1"/"$2"_R2_"$4" Trimmed_output_"$5"/"$2"_fwd_paired.fastq Trimmed_output_"$5"/"$2"_fwd_unpaired.fastq Trimmed_output_"$5"/"$2"_rev_paired.fastq Trimmed_output_"$5"/"$2"_rev_unpaired.fastq ILLUMINACLIP:adapter_run2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:100"}'`
	echo $TRIM_COMMAND_PE
	#eval "$TRIM_COMMAND_PE" 
        export KRAKEN2_DB_PATH=$TOOLS_DIRECTORY/kraken2-master
	KRAKEN2_COMMAND=`echo -e $i"\t"$OUT"\t"$TOOLS_DIRECTORY | awk '{ print "kraken2 --db minikraken_8GB_20200312 --threads 16 --report Kraken_output_"$2"/"$1".kreport --paired Trimmed_output_"$2"/"$1"_fwd_paired.fastq Trimmed_output_"$2"/"$1"_rev_paired.fastq > Kraken_output_"$2"/"$1"_trimmed.kraken"}'`
        echo $KRAKEN2_COMMAND
        #eval $KRAKEN2_COMMAND
        KRONA_COMMAND=`echo -e $i"\t"$OUT | awk '{ print "ktImportTaxonomy -s 3 -t 4 -o Krona_output_"$2"/"$1".html Kraken_output_"$2"/"$1"_trimmed.kraken"}'`
        echo $KRONA_COMMAND
        #eval $KRONA_COMMAND
	HISAT_HUMAN_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$HUMAN_REFERENCE"\t"$OUT"\t"$TOOLS_DIRECTORY | awk '{ print $5"/hisat2-2.1.0/hisat2 -x "$3" -1 Trimmed_output_"$4"/"$2"_fwd_paired.fastq -2 Trimmed_output_"$4"/"$2"_rev_paired.fastq -S Hisat2_human_output_"$4"/"$2"_human.sam -p 16 --dta-cufflinks --summary-file Hisat2_human_output_"$4"/"$2"_human.log"}'`
        echo $HISAT_HUMAN_COMMAND
        #eval $HISAT_HUMAN_COMMAND
	HISAT_COVIDALL_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$OUT"\t"$TOOLS_DIRECTORY | awk '{ print $5"/hisat2-2.1.0/hisat2 -x "$3" -1 Trimmed_output_"$4"/"$2"_fwd_paired.fastq -2 Trimmed_output_"$4"/"$2"_rev_paired.fastq -S Hisat2_covid_output_"$4"/"$2"_covidall.sam -p 16 --dta-cufflinks --summary-file Hisat2_covid_output_"$4"/"$2"_covidall.log"}'`
        echo $HISAT_COVIDALL_COMMAND
        #eval $HISAT_COVIDALL_COMMAND
        SAM2BAM_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -S -b -@ 8 Hisat2_human_output_"$3"/"$1"_human.sam | "$2"/samtools-1.10/samtools sort -@ 8 - -o Hisat2_human_output_"$3"/"$1"_human_sort.bam -O BAM"}'`
        echo $SAM2BAM_COMMAND
	#eval $SAM2BAM_COMMAND
        BAM_MAP1_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -u -f 4 -F 264 Hisat2_human_output_"$3"/"$1"_human_sort.bam > temp1.bam"}'`
        BAM_MAP2_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -u -f 8 -F 260 Hisat2_human_output_"$3"/"$1"_human_sort.bam > temp2.bam"}'`
        BAM_MAP3_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -u -f 12 -F 256 Hisat2_human_output_"$3"/"$1"_human_sort.bam > temp3.bam"}'`
        echo $BAM_MAP1_COMMAND
	echo $BAM_MAP2_COMMAND
        echo $BAM_MAP3_COMMAND
        #eval $BAM_MAP1_COMMAND
        #eval $BAM_MAP2_COMMAND
        #eval $BAM_MAP3_COMMAND
        BAM_MERGE_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools merge -u - temp[123].bam | "$2"/samtools-1.10/samtools sort -n - -o Hisat2_human_output_"$3"/"$1"_unmapped.bam"}'`
        echo $BAM_MERGE_COMMAND
        #eval $BAM_MERGE_COMMAND
        #rm temp*.bam
	BAM2FASTQ_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/bam2fastq/bam2fastq Hisat2_human_output_"$3"/"$1"_unmapped.bam --output Hisat2_human_output_"$3"/"$1"#.fq --force"}'`
        echo $BAM2FASTQ_COMMAND
        #eval $BAM2FASTQ_COMMAND
        HISAT_COVID_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$OUT"\t"$TOOLS_DIRECTORY | awk '{ print $5"/hisat2-2.1.0/hisat2 -x "$3" -1 Hisat2_human_output_"$4"/"$2"_1.fq -2 Hisat2_human_output_"$4"/"$2"_2.fq  -S Hisat2_covid_output_"$4"/"$2"_covidunmapped.sam -p 16 --dta-cufflinks --summary-file Hisat2_covid_output_"$4"/"$2"_covidunmapped.log"}'`
        echo $HISAT_COVID_COMMAND
        #eval $HISAT_COVID_COMMAND
        SAM2BAM_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -S -b -@ 8 Hisat2_covid_output_"$3"/"$1"_covidunmapped.sam | "$2"/samtools-1.10/samtools sort -@ 8 - -o Hisat2_covid_output_"$3"/"$1"_covid_sort.bam"}'`
        echo $SAM2BAM_COMMAND
        #eval $SAM2BAM_COMMAND
        COVID_FLAGSTAT_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools flagstat Hisat2_covid_output_"$3"/"$1"_covid_sort.bam > Hisat2_covid_output_"$3"/"$1"_flagstat.txt"}'`
        echo $COVID_FLAGSTAT_COMMAND
        #eval $COVID_FLAGSTAT_COMMAND
	COVID_MPILEUP=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT"\t"$COVID_REFERENCE | awk '{ print $2"/bcftools/bcftools mpileup -f "$4" Hisat2_covid_output_"$3"/"$1"_covid_sort.bam | "$2"/bcftools/bcftools call -c | "$2"/samtools-master/bcftools/vcfutils.pl vcf2fq > Variant_calling_output_"$3"/"$1"_consensus.fq"}'`
	echo $COVID_MPILEUP
	#eval $COVID_MPILEUP
	FASTQTOFASTA_COMMAND=`echo -e $i"\t"$OUT"\t"$TOOLS_DIRECTORY | awk '{ print $3"/seqtk-master/seqtk seq -aQ64 -q20 -n N Variant_calling_output_"$2"/"$1"_consensus.fq > Variant_calling_output_"$2"/"$1"_consensus.fasta"}'`
	echo $FASTQTOFASTA_COMMAND
	#eval $FASTQTOFASTA_COMMAND
	VARIANT_CALLING_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT"\t"$COVID_REFERENCE | awk '{ print $2"/bcftools/bcftools mpileup -f "$4" Hisat2_covid_output_"$3"/"$1"_covid_sort.bam | "$2"/bcftools/bcftools call -cv -Ob > Variant_calling_output_"$3"/"$1"_variant_bcftools.bcf"}'`
	echo $VARIANT_CALLING_COMMAND
	#eval $VARIANT_CALLING_COMMAND
	BCF2VCF_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/bcftools/bcftools view Variant_calling_output_"$3"/"$1"_variant_bcftools.bcf > Variant_calling_output_"$3"/"$1"_variant_bcftools.vcf"}'`
	echo $BCF2VCF_COMMAND
	#eval $BCF2VCF_COMMAND
	VARSCAN_PILEUP_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT"\t"$COVID_REFERENCE | awk '{ print $2"/samtools-1.10/samtools mpileup -f "$4" Hisat2_covid_output_"$3"/"$1"_covid_sort.bam > Hisat2_covid_output_"$3"/"$1"_covid_varscan.pileup"}'`
	echo $VARSCAN_PILEUP_COMMAND
	#eval $VARSCAN_PILEUP_COMMAND
	VARSCAN_VCF_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print "java -jar "$2"/varscan/VarScan.v2.4.4.jar mpileup2cns Hisat2_covid_output_"$3"/"$1"_covid_varscan.pileup --output-vcf 1 --variants Variant_calling_output_"$3"/"$1"_covid_varscan.vcf > Variant_calling_output_"$3"/"$1"_covid_varscan.vcf"}'`
	echo $VARSCAN_VCF_COMMAND
	#eval $VARSCAN_VCF_COMMAND
	MEGAHIT_ASSEMBLY_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 Hisat2_human_output_"$3"/"$1"_1.fq -2 Hisat2_human_output_"$3"/"$1"_2.fq -o DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_megahit"}'`
	echo $MEGAHIT_ASSEMBLY_COMMAND
	eval $MEGAHIT_ASSEMBLY_COMMAND
	SPADES_ASSEMBLY_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{print "python "$2"/SPAdes-3.14.0-Linux/bin/spades.py -t 8 -1 Hisat2_human_output_"$3"/"$1"_1.fq -2 Hisat2_human_output_"$3"/"$1"_2.fq -o DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_spades"}'`
	echo $SPADES_ASSEMBLY_COMMAND
	eval $SPADES_ASSEMBLY_COMMAND
	
	MEGAHITONCOVID_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $4"/hisat2-2.1.0/hisat2 -x "$3" -f DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_megahit/final.contigs.fasta -S DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_megahit/"$2"_covidall.sam -p 16 --dta-cufflinks --summary-file DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_megahit/"$2"_covidall.log"}'`
        echo $MEGAHITONCOVID_COMMAND
        echo "Aligning reads over COVID Reference Genome"
        eval $MEGAHITONCOVID_COMMAND
		
	SPADESONCOVID_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $4"/hisat2-2.1.0/hisat2 -x "$3" -f DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_spades/contigs.fasta -S DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_spades/"$2"_covidall.sam -p 16 --dta-cufflinks --summary-file DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_spades/"$2"_covidall.log"}'`
        echo $SPADESONCOVID_COMMAND
        echo "Aligning reads over COVID Reference Genome"
        eval $SPADESONCOVID_COMMAND
QUAST_MEGAHIT_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$COVID_REFERENCE"\t"$OUT | awk '{ print "python "$2"/quast/quast.py DenovoAssembly_output_"$4"/"$1"_unmappedDenovo_megahit/final.contigs.fa -r "$3" -1 Trimmed_output_"$4"/"$1"_fwd_paired.fastq -2 Trimmed_output_"$4"/"$1"_rev_paired.fastq -o DenovoQuality_output_"$4"/"$1"_megahit_quast_output"}'`
	echo $QUAST_MEGAHIT_COMMAND
	#eval $QUAST_MEGAHIT_COMMAND
	QUAST_SPADES_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$COVID_REFERENCE"\t"$OUT | awk '{ print "python "$2"/quast/quast.py DenovoAssembly_output_"$4"/"$1"_unmappedDenovo_spades/contigs.fasta -r "$3" -1 Trimmed_output_"$4"/"$1"_fwd_paired.fastq -2 Trimmed_output_"$4"/"$1"_rev_paired.fastq -o DenovoQuality_output_"$4"/"$1"_spades_quast_output"}'`
        echo $QUAST_SPADES_COMMAND
        #eval $QUAST_SPADES_COMMAND
	#CONSENSUS_GISAID_SEQUENCE=`echo $i"\t"$GISAID_PATH"\t"$OUT | awk '{ print "cat DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_megahit/final.contigs.fa DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_spades/contigs.fasta "$2" > Phylogeny_output_"$3"/"$1"_gisaid_cov2020_sequences.fasta"}'`
	#echo $CONSENSUS_GISAID_SEQUENCE
	#eval $CONSENSUS_GISAID_SEQUENCE 
	#MAFFT_COMMAND=`echo -e $i"\t"$OUT | awk '{ print "mafft --thread 4 Phylogeny_output_"$2"/"$1"_gisaid_cov2020_sequences.fasta > Phylogeny_output_"$2"/"$1"_gisaid_cov2020_sequences_MSA.fasta"}'`
	#echo $MAFFT_COMMAND
	#eval $MAFFT_COMMAND
	#MEGACC_COMMAND=`echo -e $i"\t"$TOOLS_DIRETORY"\t"$OUT | awk '{ print $2"/megacc -a <infer_NJ_nucleotide.mao> -d Phylogeny_output_"$3"/"$1"_gisaid_cov2020_sequences_MSA.fasta -o Phylogeny_output_"$3/"$1"_gisaid_cov2020_tree"}'`
	#echo $MEGACC_COMMAND
	#eval $MEGACC_COMMAND
done

end=`date +%s`

runtime=$((end-start))

echo ------------------------------------------------------------------
echo All done
echo Time Taken: $runtime seconds
echo ------------------------------------------------------------------



