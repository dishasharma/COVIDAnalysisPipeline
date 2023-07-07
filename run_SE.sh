#!/bin/bash 
start=`date +%s`
usage()
{
  echo "Usage: ./run_SE.sh -d <Directory of fastq files> -t <Directory where tools are installed> -r <Human_Reference_Genome.fa> -c <Covid_refernce_genome.fa> -g <GISAID Genome Path> -m <MAO file path> -n <output folder name> -h <help>"
  exit 2
}

while getopts d:t:r:c:g:n:m:h: option 
do 
 case "${option}" 
 in 
 d) DIRECTORY=${OPTARG};;
 t) TOOLS_DIRECTORY=${OPTARG};;
 r) HUMAN_REFERENCE=${OPTARG};;
 c) COVID_REFERENCE=${OPTARG};;
 g) GISAID_PATH=${OPTARG};;
 m) MAO_PATH=${OPTARG};;
 n) OUT=${OPTARG};;
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

FileList="$(ls $DIRECTORY/*.fastq.gz | awk '{ print $1}' | awk -F'.fastq' '{ print $1}' | awk -F'/' '{ print $2}')"


for i in $FileList
do
	echo $i
	echo -e $DIRECTORY"\t"$i
	#echo "fastqc "$i
        #echo $FASTQC_OUTPUT       
	FASTQC_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $3"/FastQC/fastqc "$1"/"$2".fastq* -o Fastqc_output_"$4}'`
        echo $FASTQC_COMMAND
	echo "Running Quality check for sequencing reads"
	#eval "$FASTQC_COMMAND"
        TRIM_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print "java -jar "$3"/Trimmomatic-0.39/trimmomatic-0.39.jar SE "$1"/"$2".fastq.gz Trimmed_output_"$4"/"$2"_trimmed.fastq ILLUMINACLIP:adapter_run2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50"}'`
        echo $TRIM_COMMAND
	echo "Trimming Bad Quality Reads"
	eval "$TRIM_COMMAND"
	export KRAKEN2_DB_PATH=$TOOLS_DIRECTORY/kraken2-master
        KRAKEN2_COMMAND=`echo -e $i"\t"$OUT | awk '{ print "kraken2 --db minikraken_8GB_20200312 --threads 16 --report Kraken_output_"$2"/"$1".kreport Trimmed_output_"$2"/"$1"_trimmed.fastq > Kraken_output_"$2"/"$1"_trimmed.kraken"}'`
        echo $KRAKEN2_COMMAND
        echo "Analyzing species diversity using Kraken2"
	eval $KRAKEN2_COMMAND
        KRONA_COMMAND=`echo -e $i"\t"$OUT | awk '{ print "ktImportTaxonomy -s 3 -t 4 -o Krona_output_"$2"/"$1".html Kraken_output_"$2"/"$1"_trimmed.kraken"}'`
        echo $KRONA_COMMAND
       	echo "Performing Taxonomy analysis using Krona"
	eval $KRONA_COMMAND
	HISAT_HUMAN_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$HUMAN_REFERENCE"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $4"/hisat2-2.1.0/hisat2 -x "$3" -U Trimmed_output_"$5"/"$2"_trimmed.fastq -S Hisat2_human_output_"$5"/"$2"_human.sam -p 16 --dta-cufflinks --summary-file Hisat2_human_output_"$5"/"$2"_humanUnpaired.log"}'`
        echo $HISAT_HUMAN_COMMAND
        echo "Aligning reads over Human Reference Genome"
	eval $HISAT_HUMAN_COMMAND
	HISAT_COVIDALL_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $4"/hisat2-2.1.0/hisat2 -x "$3" -U Trimmed_output_"$5"/"$2"_trimmed.fastq -S Hisat2_covid_output_"$5"/"$2"_covidall.sam -p 16 --dta-cufflinks --summary-file Hisat2_covid_output_"$5"/"$2"_covidall.log"}'`
        echo $HISAT_COVIDALL_COMMAND
        echo "Aligning reads over COVID Reference Genome"
        eval $HISAT_COVIDALL_COMMAND
        SAM2BAM_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -S -b -@ 8 Hisat2_human_output_"$3"/"$1"_human.sam | "$2"/samtools-1.10/samtools sort -@ 8 - -o Hisat2_human_output_"$3"/"$1"_human_sort.bam -O BAM"}'`
        echo $SAM2BAM_COMMAND
	echo "Converting Aligned SAM to BAM"
	eval $SAM2BAM_COMMAND
        BAM_MAP1_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -u -f 4 -F 264 Hisat2_human_output_"$3"/"$1"_human_sort.bam > temp1.bam"}'`
        BAM_MAP2_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -u -f 8 -F 260 Hisat2_human_output_"$3"/"$1"_human_sort.bam > temp2.bam"}'`
        BAM_MAP3_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -u -f 12 -F 256 Hisat2_human_output_"$3"/"$1"_human_sort.bam > temp3.bam"}'`
        echo $BAM_MAP1_COMMAND
	echo $BAM_MAP2_COMMAND
        echo $BAM_MAP3_COMMAND
        echo "Extractig Unmapped Reads"
	eval $BAM_MAP1_COMMAND
        eval $BAM_MAP2_COMMAND
        eval $BAM_MAP3_COMMAND
        BAM_MERGE_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools merge -u - temp[123].bam | "$2"/samtools-1.10/samtools sort -n - -o Hisat2_human_output_"$3"/"$1"_unmapped.bam"}'`
        echo $BAM_MERGE_COMMAND
        echo "Merging Unmapped BAMs"
	eval $BAM_MERGE_COMMAND
        rm temp*.bam
	BAM2FASTQ_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/bam2fastq/bam2fastq Hisat2_human_output_"$3"/"$1"_unmapped.bam --output Hisat2_human_output_"$3"/"$1"#.fq --force"}'`
        echo $BAM2FASTQ_COMMAND
        echo "Converting BAM to FASTQ"
	eval $BAM2FASTQ_COMMAND
	#rm Hisat2_human_output/*_1.fq Hisat2_human_output/*_2.fq
        HISAT_COVID_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $4"/hisat2-2.1.0/hisat2 -x "$3" -U Hisat2_human_output_"$5"/"$2"_M.fq -S Hisat2_covid_output_"$5"/"$2"_covid.sam -p 16 --dta-cufflinks --summary-file Hisat2_covid_output_"$5"/"$2"_covidUnpaired.log"}'`
        echo $HISAT_COVID_COMMAND
        echo "Aligning Unmapped reads to COVID reference genome"
	eval $HISAT_COVID_COMMAND
        SAM2BAM_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools view -S -b -@ 8 Hisat2_covid_output_"$3"/"$1"_covid.sam | "$2"/samtools-1.10/samtools sort -@ 8 - -o Hisat2_covid_output_"$3"/"$1"_covid_sort.bam"}'`
        echo $SAM2BAM_COMMAND
        echo "Converting Unmapped SAM to BAM"
	eval $SAM2BAM_COMMAND
        COVID_FLAGSTAT_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/samtools-1.10/samtools flagstat Hisat2_covid_output_"$3"/"$1"_covid_sort.bam > Hisat2_covid_output_"$3"/"$1"_flagstat.txt"}'`
        echo $COVID_FLAGSTAT_COMMAND
        echo "Generating Alignment statistics for COVID reference genome"
	eval $COVID_FLAGSTAT_COMMAND
	COVID_MPILEUP=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT"\t"$COVID_REFERENCE | awk '{ print $2"/bcftools/bcftools mpileup -f "$4" Hisat2_covid_output_"$3"/"$1"_covid_sort.bam | "$2"/bcftools/bcftools call -c | "$2"/samtools-master/bcftools/vcfutils.pl vcf2fq > Variant_calling_output_"$3"/"$1"_consensus.fq"}'`
	echo $COVID_MPILEUP
	echo "Generating Consensus sequence for COVID reference genome"
	eval $COVID_MPILEUP
	FASTQTOFASTA_COMMAND=`echo -e $i"\t"$OUT"\t"TOOLS_DIRECTORY | awk '{ print $3"/seqtk-master/seqtk seq -aQ64 -q20 -n N Variant_calling_output_"$2"/"$1"_consensus.fq > Variant_calling_output_"$2"/"$1"_consensus.fasta"}'`
	echo $FASTQTOFASTA_COMMAND
	echo "Consensus sequence Fastq to Fasta"
	eval $FASTQTOFASTA_COMMAND
	VARIANT_CALLING_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT"\t"$COVID_REFERENCE | awk '{ print $2"/bcftools/bcftools mpileup -f "$4" Hisat2_covid_output_"$3"/"$1"_covid_sort.bam | "$2"/bcftools/bcftools call -cv -Ob > Variant_calling_output_"$3"/"$1"_variant_bcftools.bcf"}'`
	echo $VARIANT_CALLING_COMMAND
	echo "Identifying variants for COVID Reference genome"
	eval $VARIANT_CALLING_COMMAND
	BCF2VCF_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/bcftools/bcftools view Variant_calling_output_"$3"/"$1"_variant_bcftools.bcf > Variant_calling_output_"$3"/"$1"_variant_bcftools.vcf"}'`
	echo $BCF2VCF_COMMAND
	echo "Generating Variant calling using BCFtools"
	eval $BCF2VCF_COMMAND
	VARSCAN_PILEUP_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT"\t"$COVID_REFERENCE | awk '{ print $2"/samtools-1.10/samtools mpileup -f "$4" Hisat2_covid_output_"$3"/"$1"_covid_sort.bam > Hisat2_covid_output_"$3"/"$1"_covid_varscan.pileup"}'`
	echo $VARSCAN_PILEUP_COMMAND
	echo "Generating Pileup for Variant calling using BCFtools"
	eval $VARSCAN_PILEUP_COMMAND
	VARSCAN_VCF_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print "java -jar "$2"/varscan/VarScan.v2.4.4.jar mpileup2cns Hisat2_covid_output_"$3"/"$1"_covid_varscan.pileup --output-vcf 1 --variants Variant_calling_output_"$3"/"$1"_covid_varscan.vcf > Variant_calling_output_"$3"/"$1"_covid_varscan.vcf"}'`
	echo $VARSCAN_VCF_COMMAND
	echo "Generating VCF using Varscan"
	eval $VARSCAN_VCF_COMMAND
	MEGAHIT_ASSEMBLY_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $2"/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -r Hisat2_human_output_"$3"/"$1"_M.fq -o DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_megahit"}'`
	echo $MEGAHIT_ASSEMBLY_COMMAND
	echo "Assembling reads using MEGAHIT"
	eval $MEGAHIT_ASSEMBLY_COMMAND
	SPADES_ASSEMBLY_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{print "python "$2"/SPAdes-3.14.0-Linux/bin/spades.py -t 8 -s Hisat2_human_output_"$3"/"$1"_M.fq -o DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_spades"}'`
	echo $SPADES_ASSEMBLY_COMMAND
	echo "Assembling reads using SPAdes"
	eval $SPADES_ASSEMBLY_COMMAND
	
	
	MEGAHITONCOVID_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $4"/hisat2-2.1.0/hisat2 -x "$3" -f DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_megahit/final.contigs.fasta -S DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_megahit/"$2"_covidall.sam -p 16 --dta-cufflinks --summary-file DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_megahit/"$2"_covidall.log"}'`
        echo $MEGAHITONCOVID_COMMAND
        echo "Aligning reads over COVID Reference Genome"
        eval $MEGAHITONCOVID_COMMAND
		
	SPADESONCOVID_COMMAND=`echo -e $DIRECTORY"\t"$i"\t"$COVID_REFERENCE"\t"$TOOLS_DIRECTORY"\t"$OUT | awk '{ print $4"/hisat2-2.1.0/hisat2 -x "$3" -f DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_spades/contigs.fasta -S DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_spades/"$2"_covidall.sam -p 16 --dta-cufflinks --summary-file DenovoAssembly_output_"$1"/"$2"_unmappedDenovo_spades/"$2"_covidall.log"}'`
        echo $SPADESONCOVID_COMMAND
        echo "Aligning reads over COVID Reference Genome"
        eval $SPADESONCOVID_COMMAND

	QUAST_MEGAHIT_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$COVID_REFERENCE"\t"$OUT | awk '{ print "python "$2"/quast/quast.py DenovoAssembly_output_"$4"/"$1"_unmappedDenovo_megahit/final.contigs.fa -r "$3" --single Trimmed_output_"$4"/"$1"_trimmed.fastq -o DenovoQuality_output_"$4"/"$1"_megahit_quast_output"}'`
	echo $QUAST_MEGAHIT_COMMAND
	echo "Quality check for assembly generated from MEGAHIT"
	eval $QUAST_MEGAHIT_COMMAND
	QUAST_SPADES_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$COVID_REFERENCE"\t"$OUT | awk '{ print "python "$2"/quast/quast.py DenovoAssembly_output_"$4"/"$1"_unmappedDenovo_spades/contigs.fasta -r "$3" --single Trimmed_output_"$4"/"$1"_trimmed.fastq -o DenovoQuality_output_"$4"/"$1"_spades_quast_output"}'`
        echo $QUAST_SPADES_COMMAND
        echo "Quality check for assembly generated from SPAdes"
	eval $QUAST_SPADES_COMMAND
	MEGAHIT_CAT=`echo -e $i"\t"$GISAID_PATH"\t"$OUT | awk '{ print "cat  DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_megahit/final.contigs.fa DenovoAssembly_output_"$3"/"$1"_unmappedDenovo_spades/contigs.fasta "$2" > Phylogeny_output_"$3"/"$1"_GISAID_MEGAHIT.fa"}'`
	echo $MEGAHIT_CAT
	echo "merging GISAID sequence with generated contigs from megahit"
	#eval $MEGAHIT_CAT
	MEGAHIT_MAFFT_COMMAND=`echo -e $i"\t"$OUT | awk '{ print "mafft --thread 4 Phylogeny_output_"$2"/"$1"_GISAID_MEGAHIT.fa > Phylogeny_output_"$2"/"$1"_megahit_gisaid_cov2020_sequences_MSA.fasta"}'`
	echo $MEGAHIT_MAFFT_COMMAND
	echo "Phylogeny analysis using MAFFT"
	#eval $MEGAHIT_MAFFT_COMMAND
	MEGAHIT_MEGACC_COMMAND=`echo -e $i"\t"$TOOLS_DIRECTORY"\t"$MAO_PATH"\t"$OUT | awk '{ print $2"/"megacc -a "$3" -d Phylogeny_output_"$4"/"$1"_megahit_gisaid_cov2020_sequences_MSA.fasta -o Phylogeny_output_"$4"/"$1"_megahit_gisaid_cov2020_tree"}'`
	echo $MEGAHIT_MEGACC_COMMAND
	echo "Generating tree for phylogeny analysis"
	#eval $MEGAHIT_MEGACC_COMMAND
done

end=`date +%s`

runtime=$((end-start))

echo ------------------------------------------------------------------
echo All done
echo Time Taken: $runtime seconds
echo ------------------------------------------------------------------



