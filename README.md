# Lophotrochozoa_phylogenomics_scripts_06-2015
Script used for lophotrochozoan phylogenomics--This script was tested on Ubuntu Linux. Absolute paths (where used) may need to be edited.

Lophotrochozoa_orthologyScript_HaMSTR13.sh takes the output of HAMSTR and performs several steps to remove groups and sequences that are not suitable for phylogenomic analysis.

The final product of this script is a set of trimmed amino acid alignments representing putatively orthologous groups suitable for phylogenomic analysis.

A number of programs must also be in the path including Aliscore, Alicut, MAFFT, PhyloTreePruner, HaMStR (for nentferner.pl), and the uniqHaplo.pl script.

The path for AlignmentCompare (see basal_metazoan_phylogenomics_scripts_01-2015) must be specified inside the script, and AlignmentCompare should be compiled beofer running this script

A number of variables must be modified for your purposes within the bash script. We suggest you examine the entire script carefully and modify it as needed.

Input fasta file headers must be in the following format: >orthology_group_ID|species_name_abbreviation|annotation_or_sequence_ID_information

Example: >0001|LGIG|Contig1234

Fasta headers may not include spaces or non-alphanumeric characters except for underscores (pipes are OK as field delimiters only).

