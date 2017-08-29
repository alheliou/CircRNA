


### the following varaibles do not need to be changed
DIR=$(echo $(pwd))"/"
BLASTn=../ncbi-blast-2.2.26+/bin/blastn
BLAST=../ncbi-blast-2.2.26+/bin
PICARD=../picard-tools-1.119/
i=0
CONCAT=""
SEQ_ID=$(echo $(python seqid.py $REF))
REF_SIZE=$(echo $(python seqlength.py $REF))
GBK_NEW=$OUTPUT".gbk"

if [ ! -b "$GBK_NEW" ]
then
	cp $GBK_REF $GBK_NEW
fi

echo "building a database from the reference genome"
cd $BLAST
./makeblastdb -in $REF -out $REF".db" -dbtype nucl

cd $DIR

for INPUT in $( find $INPUTDIR -iname "*.fastq")
do
preflength=${#INPUTDIR}
inputlength=${#INPUT}
OUTPUT=$OUTPUTDIR${INPUT: -$(( $inputlength - $preflength)): $(( $inputlength - $preflength))-6}
##$( find $INPUTDIR -iname "*.fastq")
CONCAT=$CONCAT$OUTPUT"_circ.txt "


### 1rst step: conversion from bam to fastq if needed
if [ ${INPUT: -4} = ".bam" ]
then 	
	echo "1rst step: from bam to fastq"
	cd $PICARD
	java -jar SamToFastq.jar INPUT=$INPUT FASTQ=$OUTPUT".fastq"

elif [ ${INPUT: -6} = ".fastq" ] || [ ${INPUT: -6} = ".fasta" ]
then 
	echo "1rst step: nothing to do the input is already in fastq or fasta format"
	 cp $INPUT $OUTPUT${INPUT: -6}	
else 
	echo "Input is not in a good format, it should be a bam or a fastq file"
fi

### 2nd step: fastq to fasta
if [ ${INPUT: -6} = ".fastq" ]
then 
	echo "conversion from fastq to fasta"
	cd $DIR
	../fastqtofasta $OUTPUT".fastq" $OUTPUT".fasta"
	sed -i -e "s/@//g" $OUTPUT".fasta"
##	rm $OUTPUT".fasta-e"
elif [ ${INPUT: -6} = ".fasta" ]
then
	echo "nothing to do the INPUT is already in fasta format"
fi

cd $DIR

### 3rd step: blastn
echo "3rd step blastn, it may takes a while"
cd $BLAST 

./blastn -task megablast -outfmt '5' -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -num_threads $NB_THREAD -out  $OUTPUT".xml" -query $OUTPUT".fasta"  -db $REF".db"


### 4th step: reads selection and writing in bam format
echo "4th step: reads selection"
cd $DIR
python blast_analysis.py $OUTPUT".xml" $OUTPUT"_outc.bam" $OUTPUT"_outl.bam" $REF $OUTPUT".fasta"

### 5th step: sorting and indexing by samtools
samtools sort $OUTPUT"_outc.bam" -o $OUTPUT"_outc.sorted.bam"
samtools index $OUTPUT"_outc.sorted.bam"
samtools sort $OUTPUT"_outl.bam" -o $OUTPUT"_outl.sorted.bam"
samtools index $OUTPUT"_outl.sorted.bam"

echo "Number of reads that align in a linear way from input file: "$INPUT
LIN=$(echo $(samtools view -c $OUTPUT"_outl.sorted.bam"))
echo $LIN

echo "Number of circular reads from input file: "$INPUT
CIRC=$(echo $(samtools view -c $OUTPUT"_outc.sorted.bam"))
echo $CIRC
CIRC=$((CIRC / 2 ))
echo "CIRC/2="$CIRC

echo "Total number of reads mapped: "$(( $CIRC + $LIN))

### 6th step : stat_circ.py, determining in which loci the circular reads match
echo "6th step: determining in which loci the circular reads match"
samtools view $OUTPUT$"_outc.sorted.bam" > $OUTPUT".tmp"


echo "python stat_circ.py "$OUTPUT".tmp" $OUTPUT"_circ.txt" $GBK_REF $GBK_NEW
echo $i
python stat_circ.py $OUTPUT".tmp" $OUTPUT"_circ.txt" $GBK_REF $GBK_NEW $REF_SIZE
rm $OUTPUT".tmp"

let i++
done
### we are now analysing all the files together

### 7th step : analyse_stat_circ_file.py 
python analyse_circ.py $CONCAT > $OUTPUT"_file.txt"
### 8th step : analyse_finale_fileV2.py
python analyse_file.py $OUTPUT"_file.txt" $OUTPUT $OUTPUT".out" $SEQ_ID $REF



