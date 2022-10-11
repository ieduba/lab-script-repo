#! /bin/bash
# Get inputs and set up variables

# Set the script to exit upon error
set -o errexit


# Get parameters from command line
trimfastqdir=$1 # first input is the directory holding the raw fastq.gz files
libsfile=$2 # second input is file of libraries: ls | cut -d"_" -f1 | uniq | grep <library prexix>  > libs.txt to create library 
ref=$3 # third input is the genome index  reference (optional)..... Grab the genome from the genomes box # to build index: bowtie2-build <fastagenome> <name of index>. Input path to index/name of index (path/X.1.bt2) 
#blacklist=$4 # fourth input is the genome reference for chr sizes (optional)
#fastaref=$5 # fifth input is the FASTA genome reference (optional)

# if no libs file exit
if [ -z "$libsfile" ]
then
    echo "Library name file required"
    exit
fi

# if no ref then use the human hg19 genome reference
if [ -z "$ref" ]
then
	echo "No provided index genome. Providing reference hg19"	
   	ref="/rugpfs/fs0/risc_lab/scratch/iduba/genomes/hg19/Bowtie2Index/hg19"
fi

if [ -z "$blacklist"]
then 
	blacklist=/rugpfs/fs0/risc_lab/scratch/iduba/senescence-data/LSexp1-ATACseq/hg19-blacklist.v2.bed
fi

# if no fastaref then use the human hg19 fasta reference
if [ -z "$fastaref" ]
then
    fastaref="/rugpfs/fs0/risc_lab/scratch/iduba/genomes/hg19/Bowtie2Index/hg19.fa"
fi 


# Set hard-coded parameters for filtering. 
mapq=30 # minimum mapping quality cutoff  

######################################################ALIGN###################################################################################

#for filename in `cat "$libsfile"` ; do
	
	# change to library directory again
#	cd  "$trimfastqdir"/"$filename"
#	echo ref= $ref
	#Assign variables
#	trimRead1="$filename"*R1*.trim.fastq
#        trimRead2="$filename"*R2*.trim.fastq
	# create bam and log file names
#	out0=`basename $trimRead1 | sed 's/.fastq/.sam/g' | sed 's/_R1//g'`
#	out1=`basename $trimRead1 | sed 's/.fastq/.bam/g' | sed 's/_R1//g'`
#  	out2=`basename $trimRead1 | sed 's/.fastq/.alignlog/g' | sed 's/_R1//g'`
#	echo out0= $out0
	# call alignment with 10 cores and output your bam to out1 as sam file
#	(bowtie2  -p 10 -x $ref -1 $trimRead1 -2 $trimRead2 -S /rugpfs/fs0/risc_lab/scratch/iduba/$out0) 2> /rugpfs/fs0/risc_lab/scratch/iduba/$out2


	# go back to parent directory
#	cd ..

#done

### REORGANIZED FILES MANUALLY THEN CONVERTED SAM TO BAM WITH HOMERANALYSIS.SH
       
###############################################################DO ENRICHMENT############################################################################################
for filename in `cat "$libsfile"`; do


        cd "$filename"
        pwd
        # declare bam file
        bamFile=`ls *.trim.bam`
#        mkdir -p tmp

        # sort bam file for filtering using picard SortSam. Appoint it to variable 'sorted'

        sorted=`echo $bamFile | sed 's/.bam/.st.bam/'`

#        picard SortSam  I=$bamFile  O=$sorted  SORT_ORDER=coordinate

        # index samfile 
#        samtools index $sorted
        ## remove mitochondria, chrY, and unwanted scaffolds
#        echo "sh02a_filter_bam.sh: Removing reads from unwanted chromosomes and scaffolds"

        # change name for files filtered with all or only mitochondrial reads removal
        out1=`echo $sorted | sed 's/.bam/.all.bam/'`
        out2=`echo $sorted | sed 's/.bam/.chrM.bam/'`

        # This is the step that does the filtering. It is looking for the corresponding lines and choosing those
#        chrs=`samtools view -H *.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v _gl | grep -v Y | grep -v hap | grep -v random | grep -v v1 | grep -v v2`

        # output filtered results
#        samtools view -b $sorted `echo $chrs` > $out1

#        echo "Filtering file $bamFile by script sh02a_filter_bam.sh" >> filtering.log
        # filter by blacklist

        # if a blacklist is provided: filter by the blacklist
        if [ -f "$blacklist" ]
        then
                echo "sh02a_filter_bam.sh: Removing blacklisted reads from $blacklist"
                out3=`echo $out1 | sed 's/.bam/.blft.bam/'`
		echo $out1
		echo $out3
                # find elements in blacklist that are within your files by using bedtools and keep just those. Output to new file out3 declare above
                bedtools intersect -v -abam $out1 -b $blacklist -wa > temp.bam

                samtools view -bh -f 0x2 temp.bam -o $out3

                echo "Blacklist filtered using file $blacklist." >> filtering.log
        else
                echo "sh02a_filter_bam.sh: Blacklist file not found or specified. Skipping blacklist filter."
                out3=$out1;
                echo "Did not filter by blacklist." >> filtering.log
        fi


        # filter by mapping quality as set above
	if [ -n "$mapq" ]
        then
                echo "sh02a_filter_bam.sh: Removing low quality reads"
                out4=`echo $out3 | sed 's/.bam/.qft.bam/'`
                echo $out4
        # only include reads with set mapping quality
                samtools view -bh -f 0x2 -q $mapq $out3 -o $out4
                echo "Filtered with mapping quality filter $mapq." >> filtering.log
        else
                echo "sh02a_filter_bam.sh: Mapping quality threshold not specified, quality filter skipped"
                out4=$out3
                echo "Did not filter by mapping quality." >> filtering.log
        fi

        # Create histogram of the size f reads using picard 
        echo "Histogram with duplicates"        
        picard CollectInsertSizeMetrics I=$out4 O=hist_data_withdups.log H=hist_graphwithdups.pdf W=1000 STOP_AFTER=50000


        echo "Duplicates"
        # remove duplicates
        echo "sh02b_remove_dups_estimate_diversity.sh: Removing duplicates"
        out5=`echo "$out4" | sed 's/.bam/.rmdup.bam/'`
        echo $out5       

        picard MarkDuplicates I=$out4 O=$out5 METRICS_FILE=dups.log REMOVE_DUPLICATES=true
        #java -jar -Djava.io.tmpdir=`pwd`/tmp /app/picard/picard-tools-1.117/MarkDuplicates.jar INPUT=$p1 OUTPUT=$out1 METRICS_FILE=dups.log REMOVE_DUPLICATES=true


        # index
        echo "index"
        samtools index $out5

        # histogram file after duplicates
        picard CollectInsertSizeMetrics I=$out5 O=hist_data_withoutdups.log H=hist_graphwithoutdups.pdf W=1000 STOP_AFTER=50000


        #clean up directory
        set +e
        mkdir -p 00_source
        mv *.chrM.bam 00_source/
        mv *.all.bam 00_source/
        mv *.st.bam 00_source/
        mv *.rmdup.bam 00_source/
        mv *.st.bam.bai 00_source/

        echo "sh02a_filter_bam.sh: Done"
        cd ..

done

####################################################CALL PEAKS################################################################################

for filename in `cat "$libsfile"`; do


        cd "$filename"
        bamFile=00_source/"$filename"*.rmdup.bam

        set -e -o nounset
        # help  
        if [ -z "$bamFile" ]
        then
                echo "This will call peaks, assumes pre-filtered bam"
                echo "<BAM> Required bam input file"
                echo "<BED> Optional blacklist regions"
        exit
        fi

        # make peakCalls dir
        peakDir='peakCalls_singles'
        if [ -d $peakDir ]
        then
                echo "Found peakCalls directory"
        else
                mkdir $peakDir
        fi

        # call peaks
                echo "Calling Peaks..."
                out2=`basename $bamFile | sed 's/.rmdup.bam//'`

                # use macs two to call peaks in rawfile
                macs2 callpeak --nomodel -t $bamFile -n $peakDir/$out2 --nolambda --keep-dup all --call-summits --slocal 10000

        cd ..
done



