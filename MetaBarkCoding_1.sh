#!/usr/bin/env bash

#######################################################################################
#										      #
#		                          MetaBarkCoding		       	      #
#									       	      #
#######################################################################################

# By Remi Petrolli (modified from Marcin Jakalski's original script)
# Article: A fine-scale spatial analysis of fungal communities on tropical tree bark unveils the epiphytic rhizosphere in orchids (2021)
# Authors: REMI PETROLLI, CONRADO AUGUSTO VIEIRA, MARCIN JAKALSKI, MELISSA F. BOCAYUVA, CLEMENT VALLE, EVERALDO DA SILVA CRUZ, MARC-ANDRÉ SELOSSE, FLORENT MARTOS, MARIA CATARINA M. KASUYA


# Bioinformatics pipeline for processing paired-end Illumina sequencing of fungal ITS-2 with two sets of primers: ITS86-F/ITS4 and ITS86-F/ITS4-Tul
# This script uses the following algorithms:

# BBMerge. Bushnell B, Rood J, Singer E. 2017. BBMerge – Accurate paired shotgun read merging via overlap. PLoS ONE 12: 1–15.
# BBDuk (sourceforge.net/projects/bbmap/)
# CUTADAPT Martin M. 2011. Cutadapt removes adapter sequences from high-throuoghput sequencing reads. EMBnet.journal: 10–12.
# SWARM 3.0 Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. 2015. Swarmv2: highly-scalable and high-resolution amplicon clustering. PeerJ: 1–12.
# VSEARCH Rognes T, Flouri T, Nichols B, Quince C, Mahé F. 2016. VSEARCH: a versatile open source tool for metagenomics. PeerJ: 1–22.
# BLASTN Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. Journal of Molecular Biology 215: 403–410.
# QIIME 1.9.1 Caporaso G, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Pena AG, Goodrich JK, Gordon JI, et al. 2010. QIIME allows analysis of high-throughput community sequencing data. Nat Methods 7: 335–336.


#######################################################
# Basic parameters for BBTools and performance
CPU=4
RAM=4

export PATH="my_path_to_softwares:$PATH"
BBToolsDir=$PATH

########################################################
# Merge paired-end sequences

FQ1=$1 # R1_file.fastq
FQ2=$2 # R2_file.fastq
out=$3 # outputfile.fastq

bbmerge.sh in1="${FQ1}" in2="${FQ2}" out="${out}" \
				   outu1="${FQ1%.*}".UNMERGEDnoQC.fastq \
				   outu2="${FQ2%.*}".UNMERGEDnoQC.fastq \
				   usequality=f threads="${CPU}" -Xmx"${RAM}"g

########################################################
# Demultiplex sequences

FQ=$1 		# merged sequences
mapping=$2	# mapping.txt (each row must contain at least the sample ID (first column, "#SampleID"), the forward and reverse barcode sequences, the forward and reverse primer sequences.
outdir=$3	# output directory

hdist=0 #  maximum number of allowed mismatches.
kdiff=0 #  if set to 0, the searched motif should have exactly the same length as the input sequence
barcodeFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $mapping | grep -nx 'barcodeFw' | cut -d: -f1) - 1 ))" # getting the index of column with forward barcode etc. "-1" is applied as array below counts from 0
primerFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $mapping | grep -nx 'primerFw' | cut -d: -f1) -1 ))"
barcodeRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $mapping | grep -nx 'barcodeRev' | cut -d: -f1) -1 ))"
primerRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $mapping | grep -nx 'primerRev' | cut -d: -f1) -1 ))"

mkdir "$outdir"
cd "$outdir"


while read -r line; do
	if [[ ! "$line" =~ ^#.* && ! "$line" =~ ^Sample.* ]]; then
		IFS=$'\t' read -r -a array <<< "$line"
		frwd="${array[${barcodeFwColumnIdx}]}${array[${primerFwColumnIdx}]}" # demultiplex with barcode+primer sequences
		rev="${array[${barcodeRevColumnsIdx}]}${array[${primerRevColumnsIdx}]}"
		sample="${array[0]}"
		echo -e "\n***** Processing: $sample\t$frwd\t$rev"
		frwdRC=$(echo "${frwd}" | tr ACGTacgt TGCAtgca | rev) # reverse complementing barcode+primer
		revRC=$(echo "${rev}" | tr ACGTacgt TGCAtgca | rev)
		
		# Depending on the sequencing procedure, the merged sequences may be both forward and reverse sequences.
		echo -e "\n***** Looking for barcode+primer in forward orientation..."
		CMD1=( "${BBToolsDir}"/bbduk.sh in=../"${FQ}" outm=stdout.fq literal="$frwd" k="$((${#frwd}-$kdiff))" hdist="${hdist}" -Xmx"${RAM}"g threads="${CPU}" overwrite=t rcomp=f )
		CMD2=( "${BBToolsDir}"/bbduk.sh in=stdin.fq outm="${sample}"-matched_FR.fastq literal="$revRC" k="$((${#revRC}-$kdiff))" hdist="${hdist}" -Xmx"${RAM}"g threads="${CPU}" overwrite=t rcomp=f int=f )
		"${CMD1[@]}" | "${CMD2[@]}" &> "nohup.${sample}_FR.out"
		echo -e "\n***** Looking for barcode+primer in reverse orientation..."
		CMD3=( "${BBToolsDir}"/bbduk.sh in=../"${FQ}" outm=stdout.fq literal="$rev" k="$((${#rev}-$kdiff))" hdist="${hdist}" -Xmx"${RAM}"g threads="${CPU}" overwrite=t rcomp=f )
		CMD4=( "${BBToolsDir}"/bbduk.sh in=stdin.fq outm="${sample}"-matched_RF.fastq literal="$frwdRC" k="$((${#frwdRC}-$kdiff))" hdist="${hdist}" -Xmx"${RAM}"g threads="${CPU}" overwrite=t rcomp=f int=f )
		"${CMD3[@]}" | "${CMD4[@]}" &> "nohup.${sample}_RF.out"
	
		cat "${sample}"-matched_FR.fastq "${sample}"-matched_RF.fastq > "${sample}"-matched_both.fastq
    fi
done < "../${mapping}"
rm ./*-matched_FR.fastq ./*-matched_RF.fastq


# Setting environmental variables
uchimeref="$PATH/uchime_ref.fasta"
qiimefasta="$PATH/qiime_ref.fasta"
qiimetxt="$PATH/qiime_taxo.txt"
mapping="mapping.txt"

primerFwColumnIdx=$(sed -n $'1s/\t/\\\n/gp' $mapping | grep -nx 'primerFw' | cut -d: -f1)
primerRevColumnsIdx=$(sed -n $'1s/\t/\\\n/gp' $mapping | grep -nx 'primerRev' | cut -d: -f1)

##########################################################################################################
# Compile reads with both primers (forward and reverse) into a single file, and filter by size and quality

mkdir -p 02.combinedReads
mkdir -p cutadapt_tmp
		
cat 00.reads/*-matched_both.fastq > 02.combinedReads/ITS86both.all.fastq
	
p=0
	
# reverse complement of all the sequences
vsearch --quiet \
--fastx_revcomp 02.combinedReads/ITS86both.all.fastq  \
--fastqout 02.combinedReads/ITS86both.all_RC.fastq

for pair in $(grep -v "^#" "${mapping}" | cut -f $primerFwColumnIdx,$primerRevColumnsIdx | tr "\t" ":" | sort | uniq ); do
	p1=$( echo "${pair}" | cut -f 1 -d ":" )
	p1rc=$( echo "${p1}" | tr ACGTacgt TGCAtgca | rev )
	p2=$( echo "${pair}" | cut -f 2 -d ":" )
	p2rc=$( echo "${p2}" | tr ACGTacgt TGCAtgca | rev )
	
	CMD1=( cutadapt -g "${p1}...${p2rc}" -O 15 -m 200 -n 1 -q 25 -e 0 --trimmed-only -o 02.combinedReads/ITS86.fullQ25.pair"${p}".fwd.fastq 02.combinedReads/ITS86both.all.fastq ) # forward orientation
	CMD2=( cutadapt -g "${p1}...${p2rc}" -O 15 -m 200 -n 1 -q 25 -e 0 --trimmed-only -o 02.combinedReads/ITS86.fullQ25.pair"${p}".rev.fastq 02.combinedReads/ITS86both.all_RC.fastq ) # reverse orientation
	
	echo "CMD: ${CMD1[@]}" && "${CMD1[@]}" >> cutadapt_tmp/cutadapt_combine_reads.fwd.log
	echo "CMD: ${CMD2[@]}" && "${CMD2[@]}" >> cutadapt_tmp/cutadapt_combine_reads.rev.log
	p=$((p + 1))
	echo
done
cat 02.combinedReads/ITS86.fullQ25.pair*.fastq > 02.combinedReads/ITS86.fullQ25.fastq
echo "***************************************************************"


function fq2fasta {
	# Conversion from fastq to fasta
	cat $1 | awk '{if(NR%4 ~ "1|2"){print $0}}' | tr '@' '>'
}

function remove_N_records {
	awk '/^>/ { printf("\n%s\t",$0);next; } { printf("%s", $0); } END { printf("\n"); }' < $1 | awk -F "\t" '$2 !~ /N|n/ {print}' | gsed -e 's/\t/\n/g' | gsed '/^$/d' | gsed -e 's/\t/\n/g'
}


# Combine both datasets and convert into fasta
cat 02.combinedReads/ITS86.fullQ25.fastq > 02.combinedReads/ALL.ITS86.fastq
fq2fasta 02.combinedReads/ALL.ITS86.fastq > 02.combinedReads/ALL.ITS86.fna
# Removing sequences with Ns, as swarm doesn't accept such input
remove_N_records 02.combinedReads/ALL.ITS86.fna > 02.combinedReads/ALL.ITS86.noNs.fna


# FastQC of all this
fastqc -t 4 02.combinedReads/*.fastq -o 02.combinedReads


##########################################################################################################

# SWARM-2 needs a dereplication step that can be achieved by VSEARCH
vsearch \
--derep_fulllength 02.combinedReads/ALL.ITS86.noNs.fna \
--sizeout \
--relabel_sha1 \
--output 02.combinedReads/ALL.ITS86.noNs.drp.fna


swarmRes=5
mkdir 03.swarmRef
# Pick SWARM otus
swarm -d $swarmRes -t 4 -z -w 03.swarmRef/ALL.ITS86.noNs.drp_otus.fna -o 03.swarmRef/ALL.ITS86.noNs.drp_otus.swarm -s 03.swarmRef/ALL.ITS86.noNs.drp_otus.stat < 02.combinedReads/ALL.ITS86.noNs.drp.fna
	
# Sort representatives and remove singletons
vsearch --fasta_width 0 \
	--sortbysize 03.swarmRef/ALL.ITS86.noNs.drp_otus.fna \
	-minsize 2 \
	--output 03.swarmRef/ALL.ITS86.noNs_otus_sort.NS.fna


# De novo chimera checking
vsearch --uchime_denovo 03.swarmRef/ALL.ITS86.noNs_otus_sort.NS.fna \
	--uchimeout 03.swarmRef/ALL.ITS86.noNs_otus.NS.uchime \
	--nonchimeras 03.swarmRef/ALL.ITS86_otus_FINAL.fna

##########################################################################################################

# Trim all the initial sequences for blasting on SWARM representatives
# This step deals with sequences of different size, especially when using two sets of primers amplifying different fragment lengths.

FIELDS=( $(grep "SampleID" "${mapping}" | gsed 's/#//g' | gsed 's/\./_/g') )
while IFS=$'\t' read "${FIELDS[@]}"; do
	if [ "${SampleID}" != "#SampleID" ]; then
		echo "Sample ${SampleID} is being processed"
		p1="${primerFw}"
		p1rc=$( echo "${p1}" | tr ACGTacgt TGCAtgca | rev )
		p2="${primerRev}"
		p2rc=$( echo "${p2}" | tr ACGTacgt TGCAtgca | rev )
		CMD1=( cutadapt -g "${p1}...${p2rc}" 00.reads/"${OriginalFileName}" -e 0 -n 2 -m 200 -q 25 -O 15 --trimmed-only -o 00.reads/"${SampleID}".filtered.fwd.fastq )
		CMD2=( cutadapt -g "${p2}...${p1rc}" 00.reads/"${OriginalFileName}" -e 0 -n 2 -m 200 -q 25 -O 15 --trimmed-only -o 00.reads/"${SampleID}".filtered.rev.fastq )
		echo CMD: "${CMD1[@]}" && "${CMD1[@]}" >> cutadaptIT86fwd.log
		echo CMD: "${CMD2[@]}" && "${CMD2[@]}" >> cutadaptIT86rev.log
		cat 00.reads/"${SampleID}".filtered.fwd.fastq 00.reads/"${SampleID}".filtered.rev.fastq > 00.reads/"${SampleID}".filtered.fastq
		fq2fasta 00.reads/"${SampleID}".filtered.fastq > 00.reads/"${InputFileName}"
		echo
	fi
done < "${mapping}"

##########################################################################################################

add_qiime_labels.py -i 00.reads -m "${mapping}" -c InputFileName -n 0 -o 00.reads
mkdir -p 04.assignOTUs
mv 00.reads/combined_seqs.fna 04.assignOTUs/AllSamples.fna

# Blast via Qiime
REFERENCESWARM="03.swarmRef_SWARM5/ALL.ITS86_otus_FINAL.fna"
parallel_pick_otus_blast.py -i 04.assignOTUs/AllSamples.fna -r $REFERENCESWARM -O $CPU -o 04.assignOTUs --similarity 0.97 # 97% threshold
cat 04.assignOTUs/AllSamples_otus.txt | gawk '{ if($3) print $0 }' > 04.assignOTUs/AllSamples_otus.NS.txt


##########################################################################################################

cp 03.swarmRef_SWARM5/ALL.ITS86_otus_FINAL.fna 04.assignOTUs/AllSamples_repset.swarm.fna
	
pick_rep_set.py -i 04.assignOTUs/AllSamples_otus.NS.txt -f 04.assignOTUs/AllSamples.fna -m most_abundant -o 04.assignOTUs/AllSamples_repset.NSma.fna
$PATH/combinerepset 04.assignOTUs/AllSamples_repset.swarm.fna 04.assignOTUs/AllSamples_repset.NSma.fna > 04.assignOTUs/AllSamples_repset.NScombined.fna
	
# Search for chimeras against UCHIME database
uchime --input 04.assignOTUs/AllSamples_repset.NScombined.fna --db $uchimeref --uchimeout 04.assignOTUs/NS.uchime --uchimealns 04.assignOTUs/NS.align
cat 04.assignOTUs/NS.uchime | cut -f 2,17 | sort -k 1,1 | grep "Y$" | cut -f 1 -d " " > 04.assignOTUs/NS_sorted.uchime
sort -k 1,1 04.assignOTUs/AllSamples_otus.NS.txt | tr "\t" " " | join -j 1 -v 1 - 04.assignOTUs/NS_sorted.uchime | tr " " "\t" > 04.assignOTUs/AllSamples_otus.OK.txt


##########################################################################################################


# assigning taxonomy
assign_taxonomy.py -i 04.assignOTUs/AllSamples_repset.NScombined.fna -r $qiimefasta -t $qiimetxt -m blast -o 04.assignOTUs


# make OTU tables
mkdir -p 05.otuTables
make_otu_table.py -i 04.assignOTUs/AllSamples_otus.OK.txt -m $mapping -t 04.assignOTUs/AllSamples_repset.NScombined_tax_assignments.txt -o 05.otuTables/SWARM.biom


# summarize OTU tables
biom summarize-table -i 05.otuTables/SWARM.biom -o 05.otuTables/SWARM.biom-summary.txt
biom summarize-table --qualitative -i 05.otuTables/SWARM.biom -o 05.otuTables/SWARM.biom-summary_qualitative.txt
biom summarize-table --observations -i 05.otuTables/SWARM.biom -o 05.otuTables/SWARM.biom-summary_observations.txt
biom convert -i 05.otuTables/SWARM.biom -o 05.otuTables/SWARM.biom.toText.txt --to-tsv
cat 05.otuTables/SWARM.biom.toText.txt | ggrep -v "^#" | while read -r line; do a=$(echo "${line}"| cut -f 1); b=$(ggrep -P "$a\t" 04.assignOTUs/AllSamples_repset.NScombined_tax_assignments.txt | cut -f 2); echo -e "${b}\t${line}"; done | sort > 05.otuTables/SWARM.biom.toText.withTaxonomy.txt
ggrep "^#" 05.otuTables/SWARM.biom.toText.txt | sed -e 's/#OTU ID/#OTU Taxonomy	OTU ID/g' | cat - 05.otuTables/SWARM.biom.toText.withTaxonomy.txt > /tmp/out && mv /tmp/out 05.otuTables/SWARM.biom.toText.withTaxonomy.txt


# Generate summary tables for OTUs
mkdir -p 06.summaryTables-2
sort_otu_table.py -i 05.otuTables/SWARM.biom -o 06.summaryTables-2/SWARM_sorted.biom
summarize_taxa.py -L 2,3,4,5,6,7 -i 06.summaryTables-2/SWARM_sorted.biom -o 06.summaryTables-2
plot_taxa_summary.py -i 06.summaryTables-2/SWARM_sorted_L2.txt,06.summaryTables-2/SWARM_sorted_L3.txt,06.summaryTables-2/SWARM_sorted_L4.txt,06.summaryTables-2/SWARM_sorted_L5.txt,06.summaryTables-2/SWARM_sorted_L6.txt,06.summaryTables-2/SWARM_sorted_L7.txt -o 06.summaryTables-2/taxa_summary_plots/ -c bar -x 20
summarize_taxa.py -L 2,3,4,5,6,7 -i 06.summaryTables-2/SWARM_sorted.biom -o 06.summaryTables-3abso -a
