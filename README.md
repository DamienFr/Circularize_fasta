# DESCRIPTION
Produces circular sequence from linear nucleotide sequence with repeated ends.
	Compare the beginning and the end of a nucleotide sequence from a fasta file in order to circularize it by cutting out the repeated segment. Comparison made by blastn. In case of multiple matches, sequence is not modified. Check log file. Multifasta files are accepted.

# USAGE

	Usage: /home/bacterio/Copy/Dropbox/xx.scripts/circularize_fasta/circularize_fasta.sh -f fasta_file [-m 0] [-l 25] [-d 2] [-k] [-h]
[] means argument is optional

available arguments:

	-f	input file, fasta format
	-m	number of mismatches allowed, default zero.
	-l	minimum length of identity between the two edges of the linear contig. Default 25.
	-d	distance of the repeated segment used to circularize from the edges of the sequence. Use this option carefully. Default 2.
	-k	activate log mode. Keep the two subfolders created for the analysis instead of removing them 
	-h	display this help

# OUTPUT

	- Circularized fasta file
  - 6 columns log file: Contig_name Circularized(yes/no) Old_length New_length Circularized_from Circularized_to

# AUTHOR

	Damien Richard, 2019
