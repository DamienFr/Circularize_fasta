#!/bin/bash


usage() { printf "\nCompare the beginning and the end of a nucleotide sequence from a fasta file in order to circularize it by cutting out the repeated segment. Comparison made by blastn. In case of multiple matches, sequence is not modified. Check log file.\nMultifasta files are accepted.\n\nUsage: $0 -f fasta_file [-m 0] [-l 25] [-d 2] [-k] [-h]
[] means argument is optional
\navailable arguments:
\n\t-f\tinput file, fasta format
\t-m\tnumber of mismatches allowed, default zero.
\t-l\tminimum length of identity between the two edges of the linear contig. Default 25.
\t-d\tdistance of the repeated segment used to circularize from the edges of the sequence. Use this option carefully. Default 2.
\t-k\tactivate log mode. Keep the two subfolders created for the analysis instead of removing them 
\t-h\tdisplay this help\n" 1>&2; exit 1; }

log_mode=false

while getopts f:m:l:d:kh? opts; do
   case ${opts} in
      f) input=${OPTARG} ;;
      m) mismatches=${OPTARG} ;;
      l) min_length=${OPTARG} ;;
      d) distance=${OPTARG} ;;
      k) log_mode=true ;;
      h|\?) usage;  exit 0 ;;
   esac
done

[ -z "$input" ]  &&  usage
[ -z "$mismatches" ]  && mismatches=0  && echo "-m (allowed mismatches) not specified, using default value 0" 
[ -z "$min_length" ]  &&  min_length=25  && echo "-l (identical segment min length) not specified, using default value 25" 
[ -z "$distance" ]  && distance=2 && echo "-d (distance of fragment from edges) not specified, using default value 2."


RED='\E[31;40m'
GREEN='\E[32;40m'

command -v makeblastdb >/dev/null 2>&1 || { printf "${RED}I require makeblastdb but it's not installed.\n to install it type into the terminal\n \"sudo apt-get install makeblastdb\" \n ${GREEN}" >&2; exit 1; }
command -v blastn >/dev/null 2>&1 || { printf "${RED}I require blastn but it's not installed.\n to install it type into the terminal\n \"sudo apt-get install blastn\" \n ${GREEN}" >&2; exit 1; }

fasta_file=$(basename "$input")
fasta_name=${fasta_file%.*}
FASTAPATH="$( cd "$(dirname "$input")" ; pwd -P )"
cd $FASTAPATH

DATE=`date +%Y-%m-%d:%H:%M:%S`

mkdir ./tmp_circularize_${fasta_name}_${DATE}

printf "\n#####################\nSanitiser of fasta file :\n";
printf "Convertion of commas to \"_\"\n";
printf "Conversion of spaces in fasta seq titles to \"_\"\n";
printf "Convertion to one lined sequences\n";
printf "Check for duplicated names\n#####################\n";
mv ${FASTAPATH}/${fasta_file} ${FASTAPATH}/${fasta_file}.original


# duplicated names !
env fasta_name=$fasta_name date=$DATE perl -ne 's/\r?\n//g; if(/^>/){ s/,|\s|\t/_/g ; $a++; if(exists($h{$a})){$error++} ; $h{$a} = $_ . "\n"}else{ $h{$a} .= $_ };END{ foreach $key (sort { $a <=> $b } keys %h){print $h{$key} . "\n"} };END{if($error){
open (ERO, ">", "./tmp_circularize_" . $ENV{fasta_name} . "_" . $ENV{date} . "/error.txt"); print ERO "I have found at least one duplicated file name in your fasta file input. commas and spaces are replaced with \"_\" so differences cannot be based on these characters.\nExiting ..." }}' ${FASTAPATH}/${fasta_file}.original > ${FASTAPATH}/${fasta_file}

if [ -f ./tmp_circularize_${fasta_name}_${DATE}/error.txt ]; then
    printf "\n## WARNING !!! ##\nI have found at least one duplicated file name in your fasta file input. commas and spaces are replaced with \"_\" so differences cannot be based on these characters.\nExiting ...\n"
    exit;
fi


printf "\n#####################\nMaking blast database ...\n#####################\n";
makeblastdb -in ${FASTAPATH}/${fasta_file} -dbtype nucl -parse_seqids -out db_circularize${DATE} -title "DB_"
printf "\n#####################\nBlasting ...\n#####################\n";
blastn  -num_alignments 1000000 -query ${FASTAPATH}/${fasta_file} -out resultat_blast.tmp -task blastn -db db_circularize${DATE} -outfmt '10 std qlen slen gaps'
echo "qseqid,sseqid,pcent_length_QUERY,pcent_length_SUBJECT,pident Percentage of identical matches,Alignment length,Number of mismatch,Number of gap openings,Start of alignment in query,End of alignment in query,Start of alignment in subject(ref),End of alignment in subject(ref),evalue,bitscore,Query_seq_len,Sbjct_seq_len, gaps" > ./tmp_circularize_${fasta_name}_${DATE}/resultat_blast
perl -F',' -ane '$qlength = 100*$F[3]/$F[12]; $slength = 100*$F[3]/$F[13]; $F[2] = $qlength . "," . $slength . "," . $F[2]; $toprint = join(",",@F);print $toprint' resultat_blast.tmp >> ./tmp_circularize_${fasta_name}_${DATE}/resultat_blast
#cat resultat_blast.tmp >> resultat_blast
rm -f resultat_blast.tmp

printf "\n#####################\nPrefiltering of blast output file ...\n#####################\n";
# gérer la présence de virgules dans les noms avec soit un parseur soit en modifiant le fasta en input !!!
# pas d'espace dans nom fasta car blast coupe au premier espace
env mism=$mismatches min_len=$min_length dist=$distance perl -F',' -ane 'if($. > 1 && $F[0] eq $F[1] && $F[6] <= $ENV{mism} && $F[2] != 100 && $F[3] != 100 && $F[5] >= $ENV{min_len} && $F[8] <= $ENV{dist} && $F[11] >= ($F[14] - $ENV{dist})){
$h1{$F[0]} = "$F[0],$F[8],$F[10],$F[14],"; $h2{$F[0]}++; push @{$h_multiple_matches{$F[0]}}, "$F[0],$F[8],$F[10],$F[14],multiple_matches\n";
};
END{
foreach $k (keys %h1){ if($h2{$k} == 1){ print $h1{$k} . "unique_match\n"}else{
# print $h1{$k} . "multiple_matches\n" ;
foreach $item (@{$h_multiple_matches{$k}}){
print $item; 
} 
}}
}' ./tmp_circularize_${fasta_name}_${DATE}/resultat_blast > ./tmp_circularize_${fasta_name}_${DATE}/resultat_blast_filtrated

printf "\n#####################\nReconstructing fasta file ...\n#####################\n";
# don't forget to print log file for those who have multiple matches ... can't have this special case happen in test seqs ...
env fasta_name=$fasta_name date=$DATE perl -ne 'BEGIN{open (GB, "<", "./tmp_circularize_" . $ENV{fasta_name} . "_" . $ENV{date} . "/resultat_blast_filtrated");while ($l = <GB>) { $l =~ s/\r?\n//g;  @fields = split /,/,$l;
$h{$fields[0]} = $l;
$fields[-1] =~ s/\r?\n// ; if($fields[-1] eq "multiple_matches"){ push @{$h_multiple_matches{$fields[0]}},$l }
}
open (LOG, ">", "log_" . $ENV{fasta_name} . "_" . $ENV{date} . "_log.txt");
print LOG "Contig_name\tCircularized\tOld_length\tNew_length\tCircularized_from\tCircularized_to\n";
};
if(/^>/){ $a = $_;$a =~ s/^>|\r?\n//g;
if(exists($h{$a})){
	@cases = split /,/,$h{$a};
	if($cases[-1] eq "multiple_matches"){
	print ">" . $a . "_uncircularized\n";
		print STDERR "Sequence $a has been labelled \"uncircularized\" multiple repeated segments that match your parameters, you have to manually choose which to consider to circularize.\nCheck log file\n\n";
		foreach $item (@{$h_multiple_matches{$a}}){
			@blast_line = split /,/,$item; 
			$debut_log= $blast_line[1] - 1; $fin_log = $blast_line[2]; $length_log = $blast_line[2]-$blast_line[1];
			print LOG $a . "\tMultiple_matches\t" . $blast_line[3] . "\t" . $length_log . "\t" . $debut_log . "\t" . $fin_log . "\n";
		}
		$circular = 0;
	}else{
	print ">" . $a . "_circularized\n"; $circular = 1}
}else{print ">" . $a . "_uncircularized\n"; $circular = 0}
}else{  
	if($circular){
		$_ =~ s/\r?\n//; $top = substr $_ ,($cases[1] - 1); $top2 = substr $top ,0,($cases[2]-$cases[1]); print $top2 . "\n";
		$debut_log= $cases[1] - 1; $fin_log = $cases[2]; $length_log = $cases[2]-$cases[1];
		print LOG $a . "\tYes\t" . $cases[3] . "\t" . $length_log . "\t" . $debut_log . "\t" . $fin_log . "\n";
	}else{
		print LOG $a . "\tNo\t" . $cases[3] . "\t\t\t\n";
		print};
;$a = ""; undef @cases}' ${FASTAPATH}/${fasta_file} > ${FASTAPATH}/${fasta_file}.circularized


# mv: cannot stat ‘/home/bacterio/Copy/Dropbox/xx.PUBLICATIONS/99.publis_done/xx.scripts/circularize_fasta/*.fasta.original’: No such file or directory

mv ${FASTAPATH}/${fasta_file} ${FASTAPATH}/${fasta_file}.tmp
mv ${FASTAPATH}/${fasta_file}.original ${FASTAPATH}/${fasta_file}

rm -f ${fasta_file}.tmp
rm -f db_circularize${DATE}*

if ! $log_mode  ; then 
rm -rf tmp_circularize_${fasta_name}_${DATE}
fi

printf "\n#####################\nEnd of the pipeline\n#####################\n\nYour output files are:\n\n-${FASTAPATH}/${fasta_file}.circularized\n-${FASTAPATH}/log_${fasta_name}_${DATE}_log.txt\n\n"



