#!/bin/bash
#SBATCH --job-name=yeastmetatranscript


ACCESSIONLIST="~/list.txt"
YEASTREF="~/masked_merged.fna"
OTHERREF="~/out_merged.fna"
ADAPTERS="~/bbmap/37.28/bin/resources/adapters.fa" #for trimming adapters from bbduk

mkdir raw_data trimmed_data x-mapped fasta blast
RAW="~/raw_data/"
TRM="~/trimmed_data/"
XMAP="~/MATSS/x-mapped/"
FASTA="~/fasta/"
BLAST="~/blast/"

module load aspera/3.7.2 bbmap/37.28 trimgalore/0.4.4 blast/2.6.0+

##############################################################################
############################ DOWNLOAD AND EXTRACT ############################
##############################################################################

#download using michael gerth sra_download.pl script
#available at https://github.com/gerthmicha/perlscripts/blob/master/sra_download.pl
perl sra_download.pl $ACCESSIONLIST

#move all files in these directories to one directory
ls -d ERR* SRR* > subfolders.txt

while read f
do
      mv ${f}/* $RAW
done <subfolders.txt

rm -r SRR* DRR* ERR*

##############################################################################
###################### CHECK AND FORMAT ENCODING  ############################
##############################################################################

#33 (Sanger), 64 (Illumina)
ls $RAW > filenames.txt

while read f
do
      testformat.sh $RAW/$f | echo `cut -f1` $f >> encoding.txt
done <filenames.txt

#reformat illumina to sanger encoding
grep "illumina" encoding.txt | sed 's/[^ ]* //' > 64enc_filenames.txt

while read f
do
	 mv $RAW/$f enc_tmp/
 	 reformat.sh in=enc_tmp/$f out=$RAW/$f qin=64 qout=33
done<64enc_filenames.txt

rm -r enc_tmp

##############################################################################
######################## QUALITY/LENGTH TRIMMING  ############################
##############################################################################

ls $RAW/*.fastq.gz | xargs -n1 basename | sed -n 's/\.fastq.gz$//p' > raw_filenames.txt

while read f
do
        bbduk.sh in=$RAW/$f.fastq.gz out=$TRM/$f-trimmed.fq.gz ref=$ADAPTERS minlen=50 qtrim=f
done <raw_filenames.txt

##############################################################################
########################## MAP WITH BBSPLIT ##################################
##############################################################################

#build a reference for bbsplit
bbsplit.sh build=2 ref_x=$YEASTREF ref_y=$OTHERREF

ls $TRM/*.fq.gz | xargs -n1 basename | sed -n 's/\-trimmed.fq.gz$//p' > trim_filenames.txt

while read f
do
        bbsplit.sh build=2 ambiguous2=toss threads=auto fast=t in=$TRM/$f-trimmed.fq.gz basename=$f-%.fq
done <trim_filenames.txt

mv *x.fq $XMAP
rm *y.fq

##############################################################################
########################## BLAST X-MAPPED READS ##############################
##############################################################################

#convert x-mapped fastq to fasta
ls $XMAP/*.fq | xargs -n1 basename | sed -n 's/\-merged_trimmed-x.fq$//p' > xmapped.txt

while read f
do
	reformat.sh in=$XMAP/${f}_trimmed-x.fq out=$FASTA/$f.fasta
done <xmapped.txt


#create blast db
makeblastdb -in fungi.fna -dbtype nucl -title fungi

ls $FASTA/*.fasta | xargs -n1 basename | sed -n 's/\.fasta$//p' > listfasta.txt

while read f
do
	blastn -db fungi.fna -query $FASTA/$f.fasta -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send evalue" -perc_identity 95 >> $BLAST/$f.txt
done <listfasta.txt

cat $BLAST/*.txt >> blastoutput.txt

Rscript top_hit.R

grep -f schizosaccharomyces_NC.txt blastfiltered.txt >> blastyeast.txt

##############################################################################
################################# END ########################################
##############################################################################


