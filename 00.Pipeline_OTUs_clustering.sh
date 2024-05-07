
#####################################################################################################
######################  CLUSTERING AND OTU TABLES FROM  METABARCODING DATASETS  #####################
#####################################################################################################



## You first need to install PYTHON3, VSEARCH, CUTADAPT, FASTQC and SWARM (v3)
## see https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline for details
# PYTHON3 in is /usr/local/bin
# FASTQC, VSEARCH, CUTADAPT and SWARM are in $path/software
# see https://manpages.debian.org/stretch/vsearch/vsearch.1.en.html for details on VSEARCH

#####################################

# Modified from :
# Perez-Lamarque, B., Petrolli, R., Strullu-Derrien, C., Strasberg, D., Morlon, H., Selosse, M. A., & Martos, F. (2022). Structure and spe- cialization of mycorrhizal networks in phylogenetically diverse tropical communities. Environmental Microbiome, 17, 38.
# https://github.com/BPerezLamarque/Scripts

####################################

path="MyPath/"
export PATH="$path/software/bin:$PATH"
cd $path

#####################################################################################################
######################  			STEP 1: MERGE R1 & R2 sequences	           ######################
#####################################################################################################

Marker="ITS"
THREADS=2

function load_dataset {

	echo "**************** Loading data"

	########   ITS   ###########

	# Load R1 and R2 files from paired-end sequencing
	INPUT_R1="R1_file.fastq"
	INPUT_R2="R2_file.fastq"


	OUTPUT="Quality.encoding_"$Marker".log"


	# Check quality encoding (33 or 64?)
	echo "**************** Check quality encoding"
	vsearch \
		--fastq_chars ${INPUT_R1} 2> ${OUTPUT/.fastq/.log}
		
	echo "**************** Checking quality with FASTQC"
	# FastQC function is required in a file named software in the $path location
	mkdir FastQC_Unmerged/
	software/FastQC/fastqc -t 4 ${INPUT_R1} -o FastQC_Unmerged/
	software/FastQC/fastqc -t 4 ${INPUT_R2} -o FastQC_Unmerged/
		
}



################################# Merge paired-reads

function merge_pairs {

	ENCODING=33
	
	INPUT_R1="R1_file.fastq"
	INPUT_R2="R2_file.fastq"

	OUTPUT="MERGED_R1R2_"$Marker".fastq"
	OUTPUT_NOTMERGED="UNMERGED_R1R2_"$Marker


	##### Merge read pairs
	echo "**************** Merge paired-reads"
	
	MIN_OVERLAP=10 # Default = 10
	MAX_DIFF=10 # Default = 10
	# Read this for information on default parameters https://forum.qiime2.org/t/default-parameters-on-vsearch-join-pairs/6501/2
	# and https://vcru.wisc.edu/simonlab/bioinformatics/programs/vsearch/vsearch_manual.pdf
	# (1) "Decreased   default   value   forfastq_minovlen option from 16 to 10." (2) "The default value for the fastq_maxdiffs option is increasedfrom  5  to  10. There  are  nowother  more  important  restrictions  that  will  avoid  merging  reads  thatcannot be reliably aligned."

	vsearch \
		--threads ${THREADS} \
		--fastq_mergepairs ${INPUT_R1} \
		--reverse ${INPUT_R2} \
		--fastq_ascii ${ENCODING} \
		--fastqout ${OUTPUT} \
		--fastqout_notmerged_fwd ${OUTPUT_NOTMERGED}"_fwd.fastq"\
		--fastqout_notmerged_rev ${OUTPUT_NOTMERGED}"_rev.fastq"\
		--fastq_allowmergestagger \
		--fastq_minovlen ${MIN_OVERLAP} \
		--fastq_maxdiffs ${MAX_DIFF} \
		--quiet 2>> ${OUTPUT/.fastq/.log}
}


#### Checking quality with FASTQC

function check_FastQC {
	# FastQC is available on https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
	
	OUTPUT="MERGED_R1R2_"$Marker".fastq"

	echo "**************** Checking quality with FASTQC"
	# FastQC function is required in a file named software in the $path location
	mkdir FastQC_Merged/
	software/FastQC/fastqc -t 4 ${OUTPUT} -o FastQC_Merged/
}
    



#####################################################################################################
######################  			STEP 2: DEMULTIPLEX SEQUENCES	           ######################
#####################################################################################################


################################# Checking mapping file

function check_mapping {

	# Info :
	# The 'mapping' file should contain at least 5 columns, the first one being the sample names (starting with '#'), and the following being, not necessary in this order but with the...
	# ... exact names : "primerFw", "primerRev", "barcodeFw", "barcodeRev".
	# The last column should not be one of these four columns. Choose instead, e.g., "NumeroPool"
	# /!\ Please do not forget to check whether the last raw is read by the function !
	# --> If not : add a space (i.e., additional line)

	cd $path
	
	echo "**************** Checking mapping file"

	Marker="ITS"

	MAPPING="mapping_file.txt"


	barcodeFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeFw' | cut -d: -f1) - 1 ))" # getting the index of column with forward barcode etc. "-1" is applied as array below counts from 0
	primerFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerFw' | cut -d: -f1) -1 ))"
	barcodeRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeRev' | cut -d: -f1) -1 ))"
	primerRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerRev' | cut -d: -f1) -1 ))"


	# Check that all the samples are presents
	while read -r line; do
		if [[ ! "$line" =~ ^#.* && ! "$line" =~ ^Sample.* ]]; then
			IFS=$'\t' read -r -a array <<< "$line"
			
			echo $line
			
			# Get sequences
			FwBarcode="${array[${barcodeFwColumnIdx}]}"
			FwPrimer="${array[${primerFwColumnIdx}]}"
			RevBarcode="${array[${barcodeRevColumnsIdx}]}"
			RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgt TGCAtgca | rev )
			#RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev ) # if reverse is degenerated
			RevPrimer="${array[${primerRevColumnsIdx}]}"
			RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgt TGCAtgca | rev )
			#RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev ) # if reverse is degenerated
			SAMPLE_NAME="${array[0]}"
			
			# Some information
			echo "${SAMPLE_NAME} is being tested.."
			echo "Barcode Fw: ${FwBarcode}"
			echo "Primer Fw: ${FwPrimer}"
			echo "Primer Rev (RC): ${RevPrimerRC}"
			echo "Barcode Rev (RC): ${RevBarcodeRC}"
		fi
		
	done < "${MAPPING}"

}

################################# Demultiplex sequences

function demultiplex_sequences {

	export PATH="$path/software/bin:$PATH"
	echo 'export PATH="$path/software/bin:$PATH"' >> ~/.bash_profile


	cd $path


	# Discard sequences containing Ns, remove if >1 error on average, add ee to sequence header, convert to fasta
	# See https://drive5.com/usearch/manual/exp_errs.html for information about the expected errors.

	 
	 echo "******************** Remove 'Ns' and error containing sequences"
	 
	vsearch \
		--fastq_filter "MERGED_APGR20_R1R2_"$Marker".fastq" \
		--fastq_maxns 0 \
		--fastq_maxee 2 \
		--fastaout "MERGED_R1R2_"$Marker".fasta"
		# --eeout \ useless here as sequences with poor quality are removed with maxee
		
	 
	# Define binaries and output files
	INPUT="MERGED_R1R2_"$Marker".fasta"
	INPUT_REVCOMP="${INPUT/.fasta/_RC.fasta}"
	MAPPING="mapping_file.txt"

	# Reverse complement fastq file

	echo "**************** Reverse complement fasta file"
	vsearch --quiet \
		--fastx_revcomp "${INPUT}" \
		--fastaout "${INPUT_REVCOMP}"


	mkdir $path/"01.Demultiplexed_data_"$Marker/

	MIN_LENGTH=200

	# Allows detection of different primer pairs in a single dataset
	echo "******************** Demultiplex sequences"
	
	barcodeFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeFw' | cut -d: -f1) - 1 ))" # getting the index of column with forward barcode etc. "-1" is applied as array below counts from 0
	primerFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerFw' | cut -d: -f1) -1 ))"
	barcodeRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeRev' | cut -d: -f1) -1 ))"
	primerRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerRev' | cut -d: -f1) -1 ))"


	while read -r line; do
		if [[ ! "$line" =~ ^#.* && ! "$line" =~ ^Sample.* ]]; then
			IFS=$'\t' read -r -a array <<< "$line"

			# Get sequences
			FwBarcode="${array[${barcodeFwColumnIdx}]}"
			FwPrimer="${array[${primerFwColumnIdx}]}"
			RevBarcode="${array[${barcodeRevColumnsIdx}]}"
			RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgt TGCAtgca | rev )
			#RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev ) # if reverse is degenerated
			RevPrimer="${array[${primerRevColumnsIdx}]}"
			RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgt TGCAtgca | rev )
			#RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev ) # if reverse is degenerated
			SAMPLE_NAME="${array[0]}"
			
			# Output file names
			LOG="01.Demultiplexed_data_"$Marker"/${SAMPLE_NAME}.log"
			FINAL_FASTA="01.Demultiplexed_data_"$Marker"/${SAMPLE_NAME}.fas"
			
			# Some information
			echo "${SAMPLE_NAME} is being processed.."
			echo "Barcode Fw: ${FwBarcode}"
			echo "Primer Fw: ${FwPrimer}"
			echo "Primer Rev (RC): ${RevPrimerRC}"
			echo "Barcode Rev (RC): ${RevBarcodeRC}"
			
			function trim_without_ambiguity {
			
				SEQTOT="${FwBarcode}${FwPrimer}${RevPrimerRC}${RevBarcodeRC}"
				MIN_MATCHED=${#SEQTOT}
				echo "Minimum length chosen: ${MIN_MATCHED}"
				#MIN_MATCHED=55 # total length = 57 (8+21+20+8) ITS86F/ITS4) and 55 (5.8SOF/ITS4-Tul)
				ERROR_RATE=0
				cat "${INPUT}" "${INPUT_REVCOMP}" | cutadapt -g "${FwBarcode}${FwPrimer}...${RevPrimerRC}${RevBarcodeRC}" --discard-untrimmed --minimum-length "${MIN_LENGTH}" -O "${MIN_MATCHED}" -e "${ERROR_RATE}" - 2> "${LOG}" > "temp_"$Marker".fasta"
			}
			
			function trim_with_ambiguities {

				MIN_MATCHED=20
				ERROR_RATE_BARCODE=0
				ERROR_RATE_PRIMER=0 # You may change this with degenerate adapters however CUTADAPT does not take into account "N" for estimating the number of error allowed
		
				# You assume that both barcodes have the same length
				
				cat "${INPUT}" "${INPUT_REVCOMP}" | cutadapt -g "${FwBarcode}...${RevBarcodeRC}" --discard-untrimmed --minimum-length "${MIN_LENGTH}" -O "${#FwBarcode}" -e "${ERROR_RATE_BARCODE}" - 2> "${LOG}" | cutadapt -g "${FwPrimer}...${RevPrimerRC}" --discard-untrimmed --minimum-length "${MIN_LENGTH}" -O "${MIN_MATCHED}" -e "${ERROR_RATE_PRIMER}" - 2>> "${LOG}" > "temp_"$Marker".fasta"
			}


			#================> Function to run: SELECT ONE <=================#
		
			trim_without_ambiguity
			#trim_with_ambiguities

		
			# Dereplicate at the study level
			# relabel_sha1 gives a unique label for each sequence, see https://manpages.debian.org/jessie-backports/vsearch/vsearch.1.en.html for details
			vsearch --quiet \
				--derep_fulllength "temp_"$Marker".fasta" \
				--sizeout \
				--fasta_width 0 \
				--relabel_sha1 \
				--output "${FINAL_FASTA}" 2>> "${LOG}"

		fi
	done < "${MAPPING}"


	# Clean
	rm -f "${INPUT}" "${INPUT_REVCOMP}"
	rm -f "temp_"$Marker".fastq" "temp_"$Marker".fasta"


}


###############################################################################################
###############################       STEP 3:  MAKE OTU at 97%        #########################
###############################################################################################

function classical_clustering {

	cd $path
	
	mkdir 02.Clustered_sequences

	cat "01_Demultiplexed_data_ITS_Total/"*.fas > "02.Clustered_sequences/reads_"$Marker".fa"

	## Dereplication using VSEARCH
	echo "******************* Dereplicate sequences"
	vsearch \
		--derep_fulllength "02.Clustered_sequences/reads_"$Marker".fa" \
		--sizein \
		--sizeout \
		--relabel_sha1 \
		--fasta_width 0 \
		--output "02.Clustered_sequences/reads_"$Marker"_derep.fa"
	# sha1 (encoding system) is giving the same names to the identical amplicons across samples
	
	rm "02.Clustered_sequences/reads_"$Marker".fa"
	

	# Sort by size
	vsearch -sortbysize "02.Clustered_sequences/reads_"$Marker"_derep.fa" -output "02.Clustered_sequences/reads_"$Marker"_sorted.fa" -minsize 1

	## OTU clustering (at 97%)
	vsearch -cluster_size  "02.Clustered_sequences/reads_"$Marker"_sorted.fa" --id 0.97 --centroids "02.Clustered_sequences/reads_"$Marker"_OTU97.fa" --uc "02.Clustered_sequences/clusters_"$Marker"_OTU97.uc" --sizein --sizeout
	

	# python script present in https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/map2qiime.py
	python3 /Users/remi/Desktop/These/Methodes/Barcoding/05.Emiliane/software/map2qiime.py "02.Clustered_sequences/clusters_"$Marker"_OTU97.uc" > "02.Clustered_sequences/reads_"$Marker"_mapped_OTU97.txt"


	# One line per OTU sequence
	vsearch --fasta_width 0 \
	--sortbysize "02.Clustered_sequences/reads_"$Marker"_OTU97.fa" \
	--output "02.Clustered_sequences/reads_"$Marker"_OTU97_final.fa"

}

function chimera_denovo_2 {

	# Make stats file
	# python script present in https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/make_stats.py
	python3 /Users/remi/Desktop/These/Methodes/Barcoding/05.Emiliane/software/make_stats.py "02.Clustered_sequences/reads_"$Marker"_OTU97_final.fa" > "02.Clustered_sequences/stats_"$Marker"_OTU97.txt"


	# Chimera filtering
	echo "**************** De novo chimera checking"
	vsearch --uchime_denovo "02.Clustered_sequences/reads_"$Marker"_OTU97_final.fa" --uchimeout "02.Clustered_sequences/reads_"$Marker"_OTU97_final.uchime"
	
}


function assign_taxonomy_2 {


	database_taxonomy="path_to_my_taxonomy.fasta"

	# Assign taxonomy
	vsearch --usearch_global "02.Clustered_sequences/reads_"$Marker"_OTU97_final.fa" \
	--threads 2 \
	--dbmask none \
	--qmask none \
	--rowlen 0 \
	--notrunclabels \
	--userfields query+id1+target \
	--maxaccepts 0 \
	--maxrejects 32 \
	--top_hits_only \
	--output_no_hits \
	--db "database/"$database_taxonomy \
	--id 0.5 \
	--iddef 4 \
	--userout "02.Clustered_sequences/taxonomy_"$Marker"_OTU97.txt"
		
}


function make_OTU_table_2 {
	
	# Make OTU table
	STATS="02.Clustered_sequences/stats_"$Marker"_OTU97.txt"
	OTUS="02.Clustered_sequences/reads_"$Marker"_mapped_OTU97.txt"
	REPRESENTATIVES="02.Clustered_sequences/reads_"$Marker"_OTU97_final.fa"
	UCHIME="02.Clustered_sequences/reads_"$Marker"_OTU97_final.uchime"
	ASSIGNMENTS="02.Clustered_sequences/taxonomy_"$Marker"_OTU97.txt"
	OTU_TABLE="02.Clustered_sequences/OTU_table_"$Marker"_OTU97.txt"

	SCRIPT="my_path/OTU_contingency_table.py"


	python3 \
		"${SCRIPT}" \
		"${REPRESENTATIVES}" \
		"${STATS}" \
		"${OTUS}" \
		"${UCHIME}" \
		"${ASSIGNMENTS}" \
		$path"01_Demultiplexed_data_ITS_Total/"*.fas > "${OTU_TABLE}"
		



	# Filter per OTU size or spread, quality and chimeric status:
	TABLE="02.Clustered_sequences/OTU_table_"$Marker"_OTU97.txt"
	FILTERED="${TABLE/.txt/_filtered.txt}"
	
	head -n 1 "${TABLE}" > "${FILTERED}"
	cat "${TABLE}" | awk '$5 == "N" && $4 >= 200 && $2 >= 10 && $6 >= 1' >> "${FILTERED}"
	# remove chimera, keep OTU of more than 200 bp, with an abundance of at least 10 reads, and spread of >=1

}


##################################################################################### Functions to select

############### STEP 1: MERGE R1 & R2 sequences

load_dataset
#merge_pairs
#check_FastQC

############### STEP 2: DEMULTIPLEX SEQUENCES

#check_mapping
#demultiplex_sequences

############### STEP 3:  MAKE OTU at 97%

#classical_clustering
#chimera_denovo_2
#assign_taxonomy_2
#make_OTU_table_2
