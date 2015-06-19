#!/bin/bash

#This script takes the output of HAMSTR and performs several "cleaning" steps to remove groups and sequences that are not suitable for phylogenomic analysis.
#The final product of this script is a set of trimmed amino acid alignments representing putatively orthologous groups.
#Note which programs must be included in $PATH
#Run this script from the working directory containing the input fasta files (extension must be .fa) where each file represents one putatively orthologous group (according to HaMStR).
#HAMSTR output fasta files must be parsed using HAMSTR_concatenate_orthologs.sh.
#Fasta headers must be in the following format: >orthology_group_ID|species_name_abbreviation|annotation_or_sequence_ID_information
#Example: >0001|LGIG|Contig1234
#Fasta headers may not include spaces or non-alphanumeric characters except for underscores (pipes are OK as field delimiters only).


########################################################################
################## Change these before your first use ##################
########################################################################
#Set values for variables.
MIN_SEQUENCE_LENGTH=50 #Deletes original seqs shorter than this length
MIN_ALIGNMENT_LENGTH=50 #Minimum length of a trimmed alignment in amino acids
MIN_TAXA=50 #Minimum number of OTUs to keep an OG
CLASS_PATH=/usr/local/bin/ #Location of PhyloTreePruner .class files
########################################################################


#Backup all sequences
echo "Making a backup of all sequences before beginning..."
mkdir unedited_sequences
cp *.fa ./unedited_sequences/
echo
echo Done
echo


#Delete sequences shorter than $MIN_SEQUENCE_LENGTH
echo "Deleting sequences shorter than $MIN_SEQUENCE_LENGTH AAs..."
for FILENAME in *.fa
do
grep -B 1 "[^>].\{$MIN_SEQUENCE_LENGTH,\}" $FILENAME > $FILENAME.out
sed -i 's/--//g' $FILENAME.out
sed -i '/^$/d' $FILENAME.out
rm -rf $FILENAME
mv $FILENAME.out $FILENAME
done
echo Done
echo


#If fewer than $MIN_TAXA different species are represented in the file, move that file to the "rejected_few_taxa" directory. 
echo "Removing groups with fewer than $MIN_TAXA taxa..."
mkdir -p rejected_few_taxa_1
for FILENAME in *.fa
do
awk -F"|" '/^>/{ taxon[$2]++ } END{for(o in  taxon){print o,taxon[o]}}' $FILENAME > $FILENAME\.taxon_count #Creates temporary file with taxon abbreviation and number of sequences for that taxon in $FILENAME
taxon_count=`grep -v 0 $FILENAME\.taxon_count | wc -l` #Counts the number of lines with an integer >0 (= the number of taxa with at least 1 sequence)
if [ "$taxon_count" -lt "$MIN_TAXA" ] ; then 
echo $FILENAME
mv $FILENAME ./rejected_few_taxa_1/
fi
done
rm -rf *.taxon_count
echo Done
echo


#List the remaining OTUs.
echo "List of OTUs:"
cat *.fa | awk -F"|" '{print $2}' | sed 's/^$//g' | sort | uniq
echo


#Remove redundant identical sequences in an alignment using uniqHaplo
echo "Removing redundant sequences using uniqHaplo..."
for FILENAME in *.fa
do
uniqHaplo.pl -a $FILENAME > $FILENAME".uniq"
rm -rf $FILENAME
mv $FILENAME".uniq" $FILENAME
done
echo Done
echo


#If one of the first 20 characters of a sequence is an X, that X and all characters before it are removed
echo "Trimming 5' ends..."
for FILENAME in *.fa
do
sed 's/^[^>]\{,19\}X//' $FILENAME > $FILENAME.trim
done
rename -f 's/.fa.trim/.fa/g' *.fa.trim
echo Done
echo


#If one of the last 20 characters of a sequence is an X, that X and all characters after it are removed
echo "Trimming 3' ends..."
for FILENAME in *.fa
do
sed '/>/! s/X.\{,19\}$//' $FILENAME > $FILENAME.trim
done
rename -f 's/.fa.trim/.fa/g' *.fa.trim
echo Done
echo


#Align the remaining sequences using Mafft.
echo "Aligning sequences using Mafft (auto)..."
mkdir backup_alignments
for FILENAME in *.fa
do
mafft --auto --localpair --maxiterate 1000 $FILENAME > $FILENAME.aln
done
rm -rf *.fa
rename 's/.fa.aln/.fa/g' *.fa.aln
cp *.fa ./backup_alignments/
echo Done
echo


#Use nentferner.pl to remove newlines.
echo "Removing linebreaks in sequences using nentferner.pl"
for FILENAME in *.fa
do
nentferner.pl -in=$FILENAME -out=$FILENAME.nent
done
rename -f 's/.fa.nent/.fa/' *.fa.nent
echo Done
echo


#Trim alignments using aliscore and alicut
echo "Trimming alignments in aliscore and alicut..."
sed -i 's/|/@/g' *.fa
for FILENAME in *.fa
do
perl /usr/local/bin/Aliscore.02.2.pl -i $FILENAME
perl /usr/local/bin/ALICUT_V2.3.pl -s
rm -rf $FILENAME
mv "ALICUT_"$FILENAME $FILENAME.out
done
rename 's/.out//g' *.out
mkdir alicut_files
mv *.txt ./alicut_files/
rm -rf *.svg
rm -rf *.xls
sed -i 's/@/|/g' *.fa


#Use nentferner.pl to remove newlines.
echo "Removing linebreaks in sequences using nentferner.pl"
for FILENAME in *.fa
do
nentferner.pl -in=$FILENAME -out=$FILENAME.nent
done
rename -f 's/.fa.nent/.fa/' *.fa.nent
echo Done
echo


#Remove highly divergent sequences using the EMBOSS program infoalign
#NOTE: headers must not contain pipe symbols!
echo "Removing highly divergent sequences using the EMBOSS program infoalign..."
sed -i 's/|/@/g' *.fa
for FILENAME in *.fa
do
ORTHOLOGY_GROUP=`echo $FILENAME | cut -d . -f 1`
infoalign -only -name -change -sequence $FILENAME -outfile $ORTHOLOGY_GROUP".info"
awk '{if ($NF <= 75) print $1;}' $ORTHOLOGY_GROUP".info" >> $ORTHOLOGY_GROUP".list" #Change the value from 75 here to whatever is needed
select_contigs.pl -e -p -l 9999 -n $ORTHOLOGY_GROUP".list" $FILENAME $FILENAME.reduced
done
mkdir infoalign_files
mv *.info ./infoalign_files/
mv *.list ./infoalign_files/
mv *.fa ./infoalign_files/
rename 's/.fa.reduced/.fa/g' *.fa.reduced
sed -i 's/@/|/g' *.fa
echo Done
echo


#Delete any stretches of 20 or fewer A-Z (but not X) characters surrounded by 10 or more gaps (-) on either side and replaces them with gaps.
echo "Deleting any incorrectly aligned misaligned sequence ends..."
tag=`cat /proc/sys/kernel/random/uuid`
for FILENAME in *.fa
do
sed  -r '/-{10,}[A-WY-Z]{,20}-{10,}/ { s//\n&\n/;s/^/'$tag'/ }' $FILENAME | sed -r '/^'$tag/' { s///;h;N;s/.*\n//;s/./-/g;H;N;s/.*\n//;H;g;s/\n//g; }' > $FILENAME.fa.cleaned
rm -rdf $FILENAME
mv $FILENAME.fa.cleaned $FILENAME
done
echo Done
echo


#Remove spaces in sequences and delete gap-only columns and columns with four or fewer non-gap characters.
echo "Removing gap-only columns..."
for FILENAME in *.fa
do
awk 'BEGIN { FS = "" }
!/^>/ { \
  sequence[NR] = $0 
  for ( i = 1; i <= NF; i++ ) \
    position[i] += ($i ~ /[A-WY-Z]/) \
} \
/^>/ { \
  header[NR] = $0 \
} \
END { \
  for ( j = 1; j <= NR; j++ ) { \
    if ( j in header) print header[j]
    if ( j in sequence ) { 
    $0 = sequence[j]
    for ( i = 1; i <= NF; i++)
    if ( position[i] > 4 ) printf "%s", $i
    printf "\n"
    } \
  } \
}' $FILENAME > $FILENAME.nogaps
done
rename -f 's/.fa.nogaps/.fa/' *.fa.nogaps
echo Done
echo


#Move trimmed alignments shorter than $MIN_ALIGNMENT_LENGTH AAs to the short_alignment folder.
echo "Moving alignments shorter than $MIN_ALIGNMENT_LENGTH AAs to the rejected_short_alignment folder..."
mkdir -p rejected_short_alignment
for FILENAME in *.fa
do
length=`awk '!/^>/{ lines++; total+= length($1) } END { average=total/(lines); printf(average);}' $FILENAME`
if [ "$length" -lt $MIN_ALIGNMENT_LENGTH ] ; then 
mv $FILENAME ./rejected_short_alignment/
fi
done
echo Done
echo


#Remove any sequences that don't overlap with all other sequences by at least 20 amino acids. This check runs until all sequences overlap with all other sequences by at least 20 amino acids.
echo "Removing short sequences that don't overalp with all other sequences by at least 20 AAs..."
for FILENAME in *.fa
do
java -cp /usr/local/bin AlignmentCompare $FILENAME
done
echo Done
echo
rm -rf myTempFile.txt


#If fewer than $MIN_TAXA different species are represented in the file, move that file to the "rejected_few_taxa_2" directory. 
echo "Removing groups with fewer than $MIN_TAXA taxa..."
mkdir -p rejected_few_taxa_2
for FILENAME in *.fa
do
awk -F "|" '/^>/{ taxon[$2]++ } END{for(o in taxon){print o,taxon[o]}}' $FILENAME > $FILENAME\.taxon_count #Creates temporary file with taxon abbreviation and number of sequences for that taxon in $FILENAME
taxon_count=`grep -v 0 $FILENAME\.taxon_count | wc -l` #Counts the number of lines with an integer >0 (= the number of taxa with at least 1 sequence)
if [ "$taxon_count" -lt "$MIN_TAXA" ] ; then 
echo $FILENAME
mv $FILENAME ./rejected_few_taxa_2/
fi
done
rm *.fa.taxon_count
echo Done
echo


#Backup the output up to this point and re-format headers for PhyloTreePruner.
echo "Making a backup of remaining .fa files and changing headers to remove OG number..."
mkdir backup_pre-phylotreepruner
cp *.fa ./backup_pre-phylotreepruner/
sed -i 's/>[0-9][0-9][0-9][0-9][0-9]|/>/g' *.fa
sed -i 's/>[0-9][0-9][0-9][0-9]|/>/g' *.fa
sed -i 's/|/@/g' *.fa
echo Done
echo


#Makes a tree for each OG with FastTree
echo "Making a tree for each OG using FastTreeMP..."
for FILENAME in *.fa
do
FastTreeMP -slow -gamma $FILENAME > $FILENAME.tre
done
rename 's/.fa.tre/.tre/g' *.fa.tre
echo Done


#Run PhyloTreePruner (TREES SHOULD BE CORRECTLY ROOTED!)
echo "Running PhyloTreePruner to remove paralogous sequences..."
	#Usage: java PhyloTreePruner input_tree_file min_number_of_taxa input_fasta bootstrap_cutoff r/u
	#r = redundant, let SCaFoS pick the best sequence for a given OTU
	#u = unique, pick the longest sequence for a given OTU
	#ex: java PhyloTreePruner 0001.tre 10 0001.fa 0.5 r

for FILENAME in *.tre
do
ORTHOLOGY_GROUP=`echo $FILENAME | cut -d . -f 1 | sed 's/.\+\///g'`
echo $ORTHOLOGY_GROUP
###################################################################################
####You may want to change some of the PhyloTreePruner settings in the next line###
###################################################################################
java -cp /usr/local/bin/PhyloTreePruner PhyloTreePruner $ORTHOLOGY_GROUP".tre" 50 $ORTHOLOGY_GROUP".fa" 0.95 u
done
echo Done
echo
echo End of script
