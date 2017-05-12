#! /bin/bash

## Setup the working directory
INDIR="/media/ab/RData/PRJNA369803_RNA-seq/"
OUTDIR=$INDIR

## Setup genome information
BT2INDEX="/media/ab/data/Research_Data/Indexes/bowtie2/hg38_Bt2/hg38"
HT2INDEX="/media/ab/data/Research_Data/Indexes/hisat2/hg38_Ht2/hg38"
CHROMSIZES="/media/ab/data/Research_Data/Genomes/hg38.chrom.sizes"
BEDGENES="/media/ab/data/Research_Data/Genomes/hg38.RefSeqGenes.bed"
GTFGENES="/media/ab/data/Research_Data/Genomes/hg38.RefSeqGenes.gtf"
GENOME="Hg38"

# Trimming settings
SE_TRIMADAPTERS="/media/ab/data/anaconda2/pkgs/trimmomatic-0.36-3/share/trimmomatic/adapters/TruSeq3-SE.fa"
PE_TRIMADAPTERS="/media/ab/data/anaconda2/pkgs/trimmomatic-0.36-3/share/trimmomatic/adapters/TruSeq3-PE.fa"
NX_TRIMADAPTERS="/media/ab/data/anaconda2/pkgs/trimmomatic-0.36-3/share/trimmomatic/adapters/NexteraPE-PE.fa"
tr_coding="-phred33"
tr_threads="2"
tr_leading="LEADING:3"
tr_trailing="TRAILING:3"
tr_slidingwindow="SLIDINGWINDOW:4:15"
tr_minlength="MINLEN:36"
tr_trimlog="" # if you want a log of all trimming, use "-trimlog $OUTDIR$datasetID"/"$datasetID".trim.log""
NEXTERA="TRUE"

# FastQC settings
qc_threads="2"
qc_contam=""

# Bowtie2 settings
bt_sensitivity="--very-sensitive"
bt_threads="2"
bt_maxfraglength="-X500"   # Set this to "-X2000" for ATAC-seq, "" or "-X500" for anything else

# Hisat2 settings
ht_threads="2"
known_splice_sites="--known-splicesite-infile /media/ab/data/Research_Data/Indexes/hisat2/hg38.hisat2.ss"
QC_filter="--qc-filter"

# Samtools settings
st_threads="2"

# FeatureCount settings
fc_type="exon"	# The feature of the GTF file to do the initial counting on
fc_attribute="gene_id" # The column of the GTF file that should be used to summarize counts
fc_level="" #use -f to count at the exon level, instead of the gene level
fc_multimap="-M" # use -M to count multi-mapping reads, or "" to ignore them
fc_threads="-T 2"
fc_minfragmentlength="-d 20"
fc_maxfragmentlength="-D 600"

# Here are the functions

trimmomatic_ () {
	if [ $read_type == "SE" ]
	then
		echo "Running Trimmomatic on unpaired reads ..."
		tr_clipping="ILLUMINACLIP:"$SE_TRIMADAPTERS":2:30:10"
		trimmomatic $read_type \
		$tr_coding \
		-threads $tr_threads \
		$tr_trimlog \
		$INDIR$filename \
		$OUTDIR$datasetID"/"$datasetID".trimmed.fastq.gz" \
		$tr_clipping \
		$tr_leading \
		$tr_trailing \
		$tr_slidingwindow \
		$tr_minlength 2>&1 | \
		tee $OUTDIR$datasetID"/"$datasetID".trimstats.log"
	fi

	if [ $read_type == "PE" ]
	then
		echo "Running Trimmomatic on paired-end reads ..."
		if [ $NEXTERA == "TRUE" ]
		then
			PE_TRIMADAPTERS=$NX_TRIMADAPTERS
			echo "Using Nextera paired-end adapters ..."			
		fi
		tr_clipping="ILLUMINACLIP:"$PE_TRIMADAPTERS":2:30:10"
		trimmomatic $read_type \
		$tr_coding \
		-threads $tr_threads \
		$tr_trimlog \
		$INDIR$filename_1 \
		$INDIR$filename_2 \
		$OUTDIR$datasetID"/"$datasetID_1".trimmed.fastq.gz" \
		$OUTDIR$datasetID"/"$datasetID_1".unpaired.fast.gz" \
		$OUTDIR$datasetID"/"$datasetID_2".trimmed.fastq.gz" \
		$OUTDIR$datasetID"/"$datasetID_2".unpaired.fastq.gz" \
		$tr_clipping \
		$tr_leading \
		$tr_trailing \
		$tr_slidingwindow \
		$tr_minlength 2>&1 | \
		tee $OUTDIR$datasetID"/"$datasetID".trimstats.log"
	fi
	gzip $OUTDIR$datasetID"/"$datasetID".trim.log"
	echo "Trimmomatic complete"
}

fastqc_ () {	
	if [ $read_type == "SE" ]
	then
		echo "Running FastQC on unpaired reads ..."		
		fastqc \
		-t $qc_threads \
		$qc_contam \
		-o $OUTDIR$datasetID"/" \
		$OUTDIR$datasetID"/"$datasetID".trimmed.fastq.gz"
	fi	

	if [ $read_type == "PE" ]
	then
		echo "Running FastQC on paired-end reads ..."
		fastqc \
		-t $qc_threads \
		$qc_contam \
		-o $OUTDIR$datasetID"/" \
		$OUTDIR$datasetID"/"$datasetID_1".trimmed.fastq.gz" \
		$OUTDIR$datasetID"/"$datasetID_2".trimmed.fastq.gz"
	fi	
	echo "FastQC complete"
}	

bowtie2_ () {
	if [ $read_type == "SE" ]
	then
		echo "Running Bowtie2 on unpaired reads ..."		
		bowtie2 \
		--un-gz $OUTDIR$datasetID"/"$datasetID"_unaligned.fastq.gz" \
		-x $BT2INDEX \
		$bt_sensitivity \
		-p $bt_threads \
		-U $OUTDIR$datasetID"/"$datasetID".trimmed.fastq.gz" \
		-S $OUTDIR$datasetID"/"$datasetID".sam" 2>&1 | \
		tee $OUTDIR$datasetID"/"$datasetID".bowtie2.log"
	fi

	if [ $read_type == "PE" ]
	then
		echo "Running Bowtie2 on paired-end reads ..."		
		bowtie2 \
		--un-conc-gz $OUTDIR$datasetID"/"$datasetID"_unaligned.fastq.gz" \
		$bt_maxfraglength \
		-x $BT2INDEX \
		$bt_sensitivity \
		-p $bt_threads \
		-1 $OUTDIR$datasetID"/"$datasetID_1".trimmed.fastq.gz" \
		-2 $OUTDIR$datasetID"/"$datasetID_2".trimmed.fastq.gz" \
		-S $OUTDIR$datasetID"/"$datasetID".sam" 2>&1 | \
		tee $OUTDIR$datasetID"/"$datasetID".bowtie2.log"
	fi
	echo "Bowtie2 complete"
}

hisat2_ () {
	if [ $read_type == "SE" ]
	then
		echo "Running Hisat2 on unpaired reads ..."		
		hisat2 \
		--un-gz $OUTDIR$datasetID"/"$datasetID"_unaligned.fastq.gz" \
		-x $HT2INDEX \
		-p $ht_threads \
		$known_splice_sites \
		$QC_filter \
		-U $OUTDIR$datasetID"/"$datasetID".trimmed.fastq.gz" \
		-S $OUTDIR$datasetID"/"$datasetID".sam" 2>&1 | \
		tee $OUTDIR$datasetID"/"$datasetID".hisat2.log"
	fi

	if [ $read_type == "PE" ]
	then
		echo "Running Hisat2 on paired-end reads ..."		
		hisat2 \
		--un-conc-gz $OUTDIR$datasetID"/"$datasetID"_unaligned.fastq.gz" \
		-x $HT2INDEX \
		-p $ht_threads \
		$known_splice_sites \
		$QC_filter \
		-1 $OUTDIR$datasetID"/"$datasetID_1".trimmed.fastq.gz" \
		-2 $OUTDIR$datasetID"/"$datasetID_2".trimmed.fastq.gz" \
		-S $OUTDIR$datasetID"/"$datasetID".sam" 2>&1 | \
		tee $OUTDIR$datasetID"/"$datasetID".hisat2.log"
	fi
	echo "Hisat2 complete"
	
}


samtools_bam_ () {
	echo "Running Samtools sort ..."
	samtools view \
	-b -u \
	-@ $st_threads \
	$OUTDIR$datasetID"/"$datasetID".sam" | \
	samtools sort \
	-@ $st_threads \
	-O bam > \
	$OUTDIR$datasetID"/"$datasetID".bam"

	# Remove the SAM file
	rm $OUTDIR$datasetID"/"$datasetID".sam"
	echo "Samtools sort complete"
}

samtools_index_ () {
	echo "Running index and flagstat ..."
	samtools index \
	$OUTDIR$datasetID"/"$datasetID".bam"
	samtools flagstat \
	-@ $st_threads \
	$OUTDIR$datasetID"/"$datasetID".bam" > \
	$OUTDIR$datasetID"/"$datasetID".flagstats.log"
	echo "Samtools index and flagstat complete"
}

picard_ () {
	echo "Running Picard to remove duplicates ..."
	picard MarkDuplicates \
	I=$OUTDIR$datasetID"/"$datasetID".bam" \
	O=$OUTDIR$datasetID"/"$datasetID".rmdup.bam" \
	M=$OUTDIR$datasetID"/"$datasetID".marked_dup_metrics.log" \
	REMOVE_DUPLICATES=true
	rm $OUTDIR$datasetID"/"$datasetID".bam"
	mv $OUTDIR$datasetID"/"$datasetID".rmdup.bam" $OUTDIR$datasetID"/"$datasetID".bam"
	echo "Picard duplicates removal complete"
}

htseq_ () {
	echo "Running Htseq-count ..."
	htseq-count \
	-f bam \
	-r pos \
	$OUTDIR$datasetID"/"$datasetID".bam" \
	$GTFGENES > \
	$OUTDIR$datasetID"/"$datasetID".htseq-counts.txt"
	echo "Htseq-count complete"
}

featurecounts_ () {
	echo "Running subread featureCounts ..."
	featureCounts \
	$fc_multimap \
	-t $fc_type \
	-g $fc_attribute \
	$fc_level \
	$fc_threads \
	$fc_minfragmentlength \
	$fc_maxfragmentlength \
	-a $GTFGENES \
	-o $OUTDIR$datasetID"/"$datasetID".featureCounts.txt" \
	$OUTDIR$datasetID"/"$datasetID".bam"
	echo "subread featureCounts complete"
}

bedtools_ () {	
	echo "Running GenomeCov to create a bedgraph file ..."
	bedtools genomecov -bg -ibam \
	$OUTDIR$datasetID"/"$datasetID".bam" \
	-strand + \
	-scale "1.0" \
	-g $CHROMSIZES > \
	$OUTDIR$datasetID"/"$datasetID".bedgraph_tmp"

	bedtools genomecov -bg -ibam \
	$OUTDIR$datasetID"/"$datasetID".bam" \
	-strand - \
	-scale "-1.0" \
	-g $CHROMSIZES >> \
	$OUTDIR$datasetID"/"$datasetID".bedgraph_tmp"

	bedtools sort -i \
	$OUTDIR$datasetID"/"$datasetID".bedgraph_tmp" > \
	$OUTDIR$datasetID"/"$datasetID".bedgraph"

	# Remove the unsorted bedgraph
	rm $OUTDIR$datasetID"/"$datasetID".bedgraph_tmp"
	echo "GenomeCov complete"
}

igv_ () {
	echo "Running IGVTools toTDF ..."
	igvtools toTDF \
	$OUTDIR$datasetID"/"$datasetID".bedgraph" \
	$OUTDIR$datasetID"/"$datasetID".tdf" \
	$CHROMSIZES
	echo "IGVTools toTDF complete"
}

rseqc_ () {
	echo "Running RSeQC ..."
	python /media/ab/data/anaconda2/pkgs/rseqc-2.6.4-py27r3.3.1_1/bin/read_distribution.py \
	-i $OUTDIR$datasetID"/"$datasetID".bam" \
	-r $BEDGENES > \
	$OUTDIR$datasetID"/"$datasetID".IERatio.log"
	echo "RSeQC complete"
}

# Get the list of files from the input directory
file_list=$(ls -p $INDIR | grep -v "/")
echo "These are the files in the directory" $INDIR":" $file_list

# Run the pipeline
for filename in $file_list
do 
	if [ ${filename: -4} == ".txt" ] || [ ${filename: -3} == ".sh" ] || [ ${filename: -4} == ".log" ] || [ ${filename: -1} == "~" ]
	then
		echo "Skipping" $filename "because it is an invalid file"		
		continue	
	fi

	read_type="SE"
	# Get the dataset root name
	filename=$(basename $filename)
	datasetID=${filename%%.*} 

	# If this is the second set of read pairs
	if [ ${datasetID: -2} == "_2" ]
	then
		echo "Skipping" $filename "because it is the 2nd set of paired-end reads"		
		continue	
	fi
	
	# Check if these are the first of a set of read pairs	
	if [ ${datasetID: -2} == "_1" ]
	then
		datasetID_1=$datasetID
		datasetID_2=${datasetID:0:-2}"_2"
		extension=${filename#*.}		
		filename_1=$datasetID_1"."$extension
		filename_2=$datasetID_2"."$extension
		datasetID=${datasetID:0:-2}		
		read_type="PE"
		echo "This is a paired-end read dataset"
		echo "Filenames:" $filename_1 $filename_2	
		echo "Datasets:" $datasetID_1 $datasetID_2
	fi
	
	if [ $read_type == "SE" ]
	then
		echo "This is an unpaired read dataset"
		echo "Filename:" $filename 
		echo "Dataset:" $datasetID
	fi

	## Make the output directory
	mkdir $OUTDIR$datasetID

	## Call the functions
#	trimmomatic_
#	fastqc_
#	bowtie2_
#	hisat2_
#	samtools_bam_
#	picard_
#	samtools_index_
#	htseq_
	featurecounts_
#	bedtools_
#	igv_
#	rseqc_
	echo "Analysis of" $datasetID "is complete"
done



### Building an index with Hisat
# python /media/ab/data/anaconda2/pkgs/hisat2-2.0.5-py27_1/bin/hisat2_extract_exons.py ./Genomes/hg38.RefSeqGenes.gtf > hg38.hisat2.exons
# python /media/ab/data/anaconda2/pkgs/hisat2-2.0.5-py27_1/bin/hisat2_extract_splice_sites.py ../../Genomes/hg38.RefSeqGenes.gtf > hg38.hisat2.ss
# hisat2-build --exon ../hg38.hisat2.exons --ss ../hg38.hisat2.ss --seed 2 /media/ab/data/Research_Data/Genomes/hg38.fa /media/ab/data/Research_Data/Indexes/hisat2/hg38_Ht2/hg38


