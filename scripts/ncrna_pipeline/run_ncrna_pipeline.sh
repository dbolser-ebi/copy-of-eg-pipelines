#!/bin/sh
# Copyright [1999-2014] EMBL-European Bioinformatics Institute
# and Wellcome Trust Sanger Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Todo: Check that the ncRNA analysis entry exist already in the database
# Maybe when missing add the ncRNA one from ../sql/ncRNA_analysis.sql

# Todo: Make sure the Mitochondrion seq_region is flagged properly for tRNAscan

NCRNA_LOGIC_NAME="ncrna_eg"

########################
#
# UPDATE ACCORDINGLY
#
########################

LSF_QUEUE="production-rh6"

# Default is toplevel
COORD_SYSTEM="toplevel"

# Default is 'CHROMOSOME'
DUMPING_TYPE="CHROMOSOME"
# the other option is 'SCAFFOLD'

###

CLEANUP=1

if [ $# != 3 ]
then
    echo "Wrong number of command line arguments"
    echo "sh run_ncrna_pipeline.sh schizosaccharomyces_pombe_core_10_63_1 division dumping_type"
    echo "division => [EPl EF EM EB EPr]"
    echo "dumping_type => [CHROMOSOME SCAFFOLD]"
    exit 1
fi

DB_NAME=$1
DIVISION=$2
DUMPING_TYPE=$3

# specific env variables in a user specific file now

ENV_FILE="${HOME}/${USER}.env.sh"

echo "Using env file: $ENV_FILE"

if [ ! -f ${ENV_FILE} ] 
then
    echo "can't find an environment file, ${ENV_FILE}, that defines the MySQL server parameters"
    echo "create one first"
    exit 1
fi

source $ENV_FILE

echo "database server host, $DB_HOST"


########################

SPECIES=`echo $DB_NAME | perl -ne '$_ =~ /^([^_]+)_([^_]+)_core.+/; $a = $1; $b = $2; print $a . "_" . $b;'`
echo "SPECIES: $SPECIES"

# e.g. anidulans
SPECIES_SHORT_NAME=`echo $SPECIES | perl -ne '$_ =~ /^(\w)[^_]+_(\w+)/; $a = $1; $b = $2; print "$a$b";'`
# e.g. ANI
SPECIES_PREFIX=`echo $SPECIES_SHORT_NAME | perl -ne '$_ =~ /^(\w\w\w).+/; $a = uc($1); print "$a";'`

echo "SPECIES_SHORT_NAME: $SPECIES_SHORT_NAME"
echo "SPECIES_PREFIX: $SPECIES_PREFIX"
echo "$DIVISION"

echo ""

OUTPUT_DIR=/nfs/nobackup2/ensemblgenomes/${USER}/ncgenes_pipelines/data/${SPECIES_SHORT_NAME}
LSF_OUTPUT=${OUTPUT_DIR}/lsf_output

CODE_ROOT_DIR=/nfs/panda/ensemblgenomes/production/ncgenes_pipelines

NCGENES_MODULES_PATH=${CODE_ROOT_DIR}/eg-pipelines/scripts/ncrna_pipeline
NCGENES_SCRIPTS_PATH=${CODE_ROOT_DIR}/eg-pipelines/scripts/ncrna_pipeline
NCGENES_SQL_PATH=${CODE_ROOT_DIR}/eg-pipelines/sql
#ENSEMBL_PATH=/nfs/panda/ensemblgenomes/apis/ensembl/69/ensembl
ENSEMBL_PATH=${CODE_ROOT_DIR}/ensembl-head
ENSEMBL_ANALYSIS_PATH=/nfs/panda/ensemblgenomes/apis/ensembl/master/ensembl-analysis/
PERL_PATH=/nfs/panda/ensemblgenomes/perl/
BIOPERL_PATH=/nfs/panda/ensemblgenomes/apis/bioperl/stable/

RNAMMER_PATH=/nfs/panda/ensemblgenomes/external/rnammer/rnammer
# old style
TRNASCAN_PATH=/nfs/panda/ensemblgenomes/external/tRNAscan-SE-1.3.1/bin
# now:
#TRNASCAN_PATH=/sw/arch/
#TRNASCAN_PATH=/nfs/panda/ensemblgenomes/external/tRNAscan-SE-1.3.1/
TRNASCAN_BIN=tRNAscan-SE

RFAMSCAN_PATH=/nfs/panda/ensemblgenomes/external/rfam_scan/rfam_scan.pl

export PATH=/nfs/panda/ensemblgenomes/perl/perlbrew/perls/5.14.2/bin:${TRNASCAN_PATH}:$PATH

export PERL5LIB=${NCGENES_MODULES_PATH}:${ENSEMBL_PATH}/modules:${BIOPERL_PATH}:${TRNASCAN_PATH}


# Rfam 10.1
# RFAM_DB_PATH=/nas/seqdb/integr8/production/data/mirror/data/Rfam/
# Rfam 11.0 and later
RFAM_DB_PATH=/nfs/panda/ensemblgenomes/production/ncgenes_pipelines/data/Rfam

if [ ! -d "${LSF_OUTPUT}" ]
then
    echo "Creating directory ${LSF_OUTPUT}"
    mkdir -p ${LSF_OUTPUT}/trnascan
    mkdir ${LSF_OUTPUT}/rnammer
    mkdir ${LSF_OUTPUT}/rfamscan
fi

if [ ! -d "${OUTPUT_DIR}/unmasked_seq" ]
then
    echo "Creating directory ${OUTPUT_DIR}/unmasked_seq"
    mkdir ${OUTPUT_DIR}/unmasked_seq
    
    if [ ${DUMPING_TYPE} == "CHROMOSOME" ] 
    then

        # Dump the sequences - one file - one sequence

	cd ${ENSEMBL_ANALYSIS_PATH}/scripts/
	echo "Dumping genomic sequences"

	echo "perl sequence_dump.pl -dbhost $DB_HOST -dbport $DB_PORT -dbuser $DB_USER -dbpass $DB_PASS -dbname $DB_NAME -coord_system_name $COORD_SYSTEM -output_dir ${OUTPUT_DIR}/unmasked_seq"
    
	perl sequence_dump.pl -dbhost $DB_HOST -dbport $DB_PORT -dbuser $DB_USER -dbpass $DB_PASS -dbname $DB_NAME -coord_system_name $COORD_SYSTEM -output_dir ${OUTPUT_DIR}/unmasked_seq

    else

	# Dump all sequences into one file - then split it up, 100 sequences per file
    
	cd ${ENSEMBL_ANALYSIS_PATH}/scripts/
	echo "Dumping genomic sequences"

	echo "perl sequence_dump.pl -dbhost $DB_HOST -dbport $DB_PORT -dbuser $DB_USER -dbpass $DB_PASS -dbname $DB_NAME -coord_system_name $COORD_SYSTEM -output_dir ${OUTPUT_DIR} -onefile"

	perl sequence_dump.pl -dbhost $DB_HOST -dbport $DB_PORT -dbuser $DB_USER -dbpass $DB_PASS -dbname $DB_NAME -coord_system_name $COORD_SYSTEM -output_dir ${OUTPUT_DIR} -onefile

	/nfs/panda/ensemblgenomes/external/bin/FastaToTbl ${OUTPUT_DIR}/toplevel.fa > ${OUTPUT_DIR}/toplevel.tbl
	cd ${OUTPUT_DIR}/unmasked_seq
	split -l 100 ../toplevel.tbl
	for f in `ls`
	do
	    echo "mv $f $f".tbl""
	    mv $f $f".tbl"
	    /nfs/panda/ensemblgenomes/external/bin/TblToFasta $f".tbl" > $f".fa"
	    rm -f $f".tbl"
	done
    fi
fi

if [ ! -d "${OUTPUT_DIR}/temp/trnascan" ]
then
    mkdir -p ${OUTPUT_DIR}/temp/trnascan
    mkdir ${OUTPUT_DIR}/temp/rnammer
    mkdir ${OUTPUT_DIR}/temp/rfamscan
fi

cd ${OUTPUT_DIR}

gff3_output=${SPECIES_SHORT_NAME}.ncRNAs.gff3

if [ -f "$gff3_output" ]
then
    echo "Deleting gff3 file, $gff3_output"
    rm -f $gff3_output
fi

touch $gff3_output

echo "getting the list of fasta files"

INDEX=0
for f in `find ${OUTPUT_DIR}/unmasked_seq -name "*.fa"`
do
    echo "Processing file, $f"
    echo ""

    INDEX=`expr $INDEX + 1`

    sequence_name=`basename $f .fa`
    trnascan_out=$sequence_name".trnascan"
    trnascan_path="temp/trnascan/"$trnascan_out
    trnascan_gff3=$sequence_name".trnascan.gff3"
    trnascan_options=""

    if [ -f "$trnascan_path" ]
    then
	echo "deleting $trnascan_path"
	rm -f $trnascan_path
    fi

    # Todo: Update the seq_region_name for MT

    if [ `basename $f .fa` == "MT" ]
    then
	echo "trnascan organelle mode"
	trnascan_options="-O"
    fi
    
    echo "Running tRNAScan-SE and parsing its results"
    echo "$TRNASCAN_BIN -o $trnascan_path $trnascan_options $f"
    echo "perl ${NCGENES_SCRIPTS_PATH}/trnascan_to_gff3.pl $trnascan_path `basename $f .fa` > $trnascan_gff3"

    bsub -q $LSF_QUEUE -J "GENEPRED"$INDEX -o ${LSF_OUTPUT}"/trnascan/"$sequence_name".trnascan.lsf.out" "$TRNASCAN_BIN -o $trnascan_path $trnascan_options $f; perl ${NCGENES_SCRIPTS_PATH}/trnascan_to_gff3.pl $trnascan_path `basename $f .fa` > $trnascan_gff3"

    INDEX=`expr $INDEX + 1`

	
    # Running Rfamscan
	
    rfamscan_out=`basename $f .fa`.rfamscan
    rfamscan_path="temp/rfamscan/"$rfamscan_out
    rfamscan_gff3=`basename $f .fa`.rfamscan.gff3
    rfamscan_options=""

    echo "Running Rfamscan and parsing its results"
    echo "perl $RFAMSCAN_PATH -o $rfamscan_path --nobig -v -filter wu --masking --blastdb ${RFAM_DB_PATH}/Rfam.fasta ${RFAM_DB_PATH}/Rfam.cm $f; perl ${NCGENES_SCRIPTS_PATH}/rfamscan10_to_gff3.pl $rfamscan_path `basename $f .fa` > $rfamscan_gff3"
    # Tell LSF we will use 4 CPUs (because wublast will) and more memory requirements

    bsub -M 8192 -R "rusage[mem=8192]" -n 4 -q $LSF_QUEUE -J "GENEPRED"$INDEX -o ${LSF_OUTPUT}"/rfamscan/"$sequence_name".rfamscan.lsf.out" "perl $RFAMSCAN_PATH -o $rfamscan_path --nobig -v -filter wu --masking --blastdb ${RFAM_DB_PATH}/Rfam.fasta ${RFAM_DB_PATH}/Rfam.cm $f; perl ${NCGENES_SCRIPTS_PATH}/rfamscan10_to_gff3.pl $rfamscan_path `basename $f .fa` > $rfamscan_gff3"
    

    rnammer_out=`basename $f .fa`.rnammer
    rnammer_path="temp/rnammer/"$rnammer_out
    rnammer_gff3=`basename $f .fa`.rnammer.gff3

    echo "Running RNAmmer and parsing its results"
    echo "$RNAMMER_PATH -T /tmp/ -S euk -m lsu,ssu,tsu -gff $rnammer_path -h "$rnammer_path".hmmreport "$f
    echo "perl ${NCGENES_SCRIPTS_PATH}/rnammer_to_gff3.pl $rnammer_path > $rnammer_gff3"

    # in rnammer wrapper script, hmmsearch configured to use 2 CPUs so tell LSF to allocate 2 CPUs

    bsub -n 2 -q $LSF_QUEUE -J "GENEPRED"$INDEX -o ${LSF_OUTPUT}"/rnammer/"$sequence_name".rnammer.lsf.out" "$RNAMMER_PATH -T /tmp/ -S euk -m lsu,ssu,tsu -gff $rnammer_path -h $rnammer_path".hmmreport" $f; perl ${NCGENES_SCRIPTS_PATH}/rnammer_to_gff3.pl $rnammer_path `basename $f .fa` > $rnammer_gff3"
	

    INDEX=`expr $INDEX + 1`

    echo ""
    
done

##
# Check that the jobs have finished
# while they haven't, don't carry on!
##

# Need to redirect STDERR into STDOUT, that's why we pipe the STDERR |& !

JOBS=`bjobs -w 2>&1| perl -ne 'if ($_ eq "No unfinished job found\n") {print "$_"; exit 0;} if ($_ =~ /^[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+(GENEPRED\d+)\s+.+/) { $job = $1; print "GENEPRED still running"; exit 0;}'`

echo "JOBS: $JOBS"

#if [ "$JOBS" == "" ]
#then
#    echo "no jobs matching 'GENEPRED'"
#fi

while [ "$JOBS" == "GENEPRED still running" ]
do    
    # echo "waiting 20 seconds"
    # wait 20 seconds
    
    sleep 20
    
    JOBS=`bjobs -w 2>&1| perl -ne 'if ($_ eq "No unfinished job found\n") {print "$_"; exit 0;} if ($_ =~ /^[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+[^\s]+\s+(GENEPRED\d+)\s+.+/) { $job = $1; print "GENEPRED still running"; exit 0;}'` 

    # echo "JOBS: $JOBS"
done

echo "JOBS: $JOBS"
echo "Gene prediction jobs execution done"

echo ""

# cleanup

if [ "$CLEANUP" == "1" ]
then
   
   echo "performing cleanup"
   
   bsub -q production-rh6 -o ${LSF_OUTPUT}/cleanup.lsf.out 'find temp -name "*.*scan" -exec rm {} \;'
   find . -name "*.rnammer" -exec rm {} \;
fi

# Generate a unique file

echo "Generating a unique with all gene predictions"

for f in `find ${OUTPUT_DIR}/unmasked_seq -name "*.fa"` 
do
    rnammer_gff3=`basename $f .fa`.rnammer.gff3
    cat $rnammer_gff3 >> $gff3_output
    trnascan_gff3=`basename $f .fa`.trnascan.gff3
    cat $trnascan_gff3 >> $gff3_output
    rfamscan_gff3=`basename $f .fa`.rfamscan.gff3
    cat $rfamscan_gff3 >> $gff3_output
done

# Clean the sequence line to only have the chromosome name

perl -i.bak -pe 'if ($_ =~ /^[^:]+:[^:]*:([^:]+):/) { my $chr = $1; $_ =~ s/^[^\t]+/$chr/; }' $gff3_output

echo ""

# Loading

echo "Loading gene predictions into db, $DB_NAME"
echo "perl ${NCGENES_SCRIPTS_PATH}/gene_store.pl -host $DB_HOST -port $DB_PORT -user $DB_USER -pass $DB_PASS -dbname $DB_NAME -gff3 $gff3_output -ncrna_gene_analysis=${NCRNA_LOGIC_NAME} -coord_system=${COORD_SYSTEM}"

perl ${NCGENES_SCRIPTS_PATH}/gene_store.pl -host $DB_HOST -port $DB_PORT -user $DB_USER -pass $DB_PASS -dbname $DB_NAME -gff3 ${gff3_output} -ncrna_gene_analysis=${NCRNA_LOGIC_NAME} -coord_system=${COORD_SYSTEM}

# Set these genes as 'NOVEL'

cd $NCGENES_SQL_PATH
cat set_genes_as_novel.sql | mysql -h $DB_HOST -P $DB_PORT -u $DB_USER -p$DB_PASS $DB_NAME

# generate_stable_ids

echo "perl ${NCGENES_SCRIPTS_PATH}/generate_ncrna_stable_ids.pl -dbuser $DB_USER -dbhost $DB_HOST -dbport $DB_PORT -dbpass $DB_PASS -dbname $DB_NAME -start ${DIVISION}${SPECIES_PREFIX}00000000000"

perl ${NCGENES_SCRIPTS_PATH}/generate_ncrna_stable_ids.pl -dbuser $DB_USER -dbhost $DB_HOST -dbport $DB_PORT -dbpass $DB_PASS -dbname $DB_NAME -start "${DIVISION}""${SPECIES_PREFIX}"00000000000

# Adding RFAM xrefs to RNAmmer rRNA genes

echo "perl ${NCGENES_SCRIPTS_PATH}/add_rfam_xrefs_to_rRNAs.pl -dbuser $DB_USER -dbhost $DB_HOST -dbport $DB_PORT -dbpass $DB_PASS -species $SPECIES"

perl ${NCGENES_SCRIPTS_PATH}/add_rfam_xrefs_to_rRNAs.pl -dbuser $DB_USER -dbhost $DB_HOST -dbport $DB_PORT -dbpass $DB_PASS -species $SPECIES

