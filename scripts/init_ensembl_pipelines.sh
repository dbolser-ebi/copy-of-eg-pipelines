#!/bin/sh

# Setup the dna and protein pipelines

# Todo: Problem, needs the right species name for the repeatmasker command line!

if [ $# != 2 ]
then
    echo "Wrong number of command line arguments"
    echo "sh init_ensembl_pipelines.sh saccharomyces_cerevisiae_core_9_62_3 chromosome|supercontig|scaffold|toplevel"
    exit 1
fi

DB_NAME=$1
COORD_SYSTEM=$2

echo "DB_NAME: $DB_NAME"
echo "COORD_SYSTEM: $COORD_SYSTEM"

if [ "$COORD_SYSTEM" != "chromosome" -a "$COORD_SYSTEM" != "supercontig" -a "$COORD_SYSTEM" != "scaffold" -a "$COORD_SYSTEM" != "toplevel" ]
then
    echo "wrong coord_system argument"
    echo "sh init_ensembl_pipelines.sh saccharomyces_cerevisiae_core_9_62_3 chromosome|supercontig|scaffold|toplevel"
    exit 1
fi

# One file specific for each user, based on ${USER}

ENV_FILE="${USER}.env.sh"

echo "Using env file: $ENV_FILE"

if [ ! -f ${ENV_FILE} ] 
then
    echo "can't find an environment file, ${ENV_FILE}, that defines the MySQL server parameters"
    echo "create one first"
    exit 1
fi

source $ENV_FILE

echo "database server host, $DB_HOST"

# Overwrite CONFIG_DIR to the generic location

CONFIG_DIR="/nfs/panda/ensemblgenomes/production/ensembl_pipelines_init/config/"

echo "CONFIG_DIR: $CONFIG_DIR


SPECIES=`echo $DB_NAME | perl -ne '$_ =~ /^([^_]+)_([^_]+)_core.+/; $a = $1; $b = $2; print $a . "_" . $b;'`
echo "SPECIES: $SPECIES"

# e.g. anidulans
SPECIES_SHORT_NAME=`echo $SPECIES | perl -ne '$_ =~ /^(\w)[^_]+_(\w+)/; $a = $1; $b = $2; print "$a$b";'`

echo "SPECIES_SHORT_NAME: $SPECIES_SHORT_NAME"

echo ""

PROTEIN_OUTPUT_DIR=/nfs/panda/ensemblgenomes/production/protein_pipelines/data/${SPECIES_SHORT_NAME}

PERL_PATH=/nfs/panda/ensemblgenomes/perl/
#ENS_PATH=/nfs/panda/ensemblgenomes/apis/ensembl/64/ensembl
ENS_PATH=/nfs/panda/ensemblgenomes/production/ensembl_pipelines_init/ensembl-head
ENS_ANALYSIS_PATH=/nfs/panda/ensemblgenomes/apis/ensembl/analysis/head
ENS_PIPELINE_PATH=/nfs/panda/ensemblgenomes/apis/ensembl/pipeline/head
BIOPERL_PATH=/nfs/panda/ensemblgenomes/apis/bioperl/stable/

export PERL5LIB=${CONFIG_DIR}:${ENS_PATH}/modules:${ENS_ANALYSIS_PATH}/modules:${ENS_PIPELINE_PATH}/modules:${ENS_PIPELINE_PATH}/scripts/:${BIOPERL_PATH}

# Add ensgen perl binary path
export PATH=${PERL_PATH}/perlbrew/perls/5.14.2/bin:/nfs/panda/ensemblgenomes/external/bin:$PATH


# 1/ cat the tables.sql

echo "cat ensembl-pipeline table.sql"

cd ${ENS_PIPELINE_PATH}/sql 

cat table.sql | mysql -h $DB_HOST -u $DB_USER -p$DB_PASS -P $DB_PORT $DB_NAME


# 2/ analysis_setup.pl

cd ${ENS_PIPELINE_PATH}/scripts

echo "perl analysis_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/rawcomputes.analysis"

perl analysis_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/rawcomputes.analysis

echo "perl analysis_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/protein_pipelines.analysis"

perl analysis_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/protein_pipelines.analysis

# 3/ rules_setup.pl

echo "perl rule_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/rawcomputes.rules"

perl rule_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/rawcomputes.rules

echo "perl rule_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/protein_pipelines.rules"

perl rule_setup.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -read -file ${CONFIG_DIR}/protein_pipelines.rules


# 4/ Dump the protein sequences and chunk them

echo "Dump the protein sequences"

cd ${ENS_PIPELINE_PATH}/scripts/protein_pipeline

PROTEIN_FILE=${SPECIES_SHORT_NAME}.pep
PROTEIN_FILE_PATH=${PROTEIN_OUTPUT_DIR}/proteins/${PROTEIN_FILE}

if [ ! -d "${PROTEIN_OUTPUT_DIR}/proteins" ]
then
    echo "mkdir -p ${PROTEIN_OUTPUT_DIR}/proteins/"
    mkdir -p ${PROTEIN_OUTPUT_DIR}/proteins
fi

if [ ! -f "$PROTEIN_FILE_PATH" ]
then
    echo "perl ./dump_translations.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -db_id > $PROTEIN_FILE_PATH"
    perl ./dump_translations.pl -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -db_id > $PROTEIN_FILE_PATH
fi

echo "Chunk the protein file"

if [ ! -d "${PROTEIN_OUTPUT_DIR}/protein_chunks" ]
then
    echo "mkdir -p ${PROTEIN_OUTPUT_DIR}/protein_chunks"
    mkdir ${PROTEIN_OUTPUT_DIR}/protein_chunks
fi

# Define how many chunks to produce, with 100 proteins per chunk

NB_PROTEIN=`grep -c '^>' $PROTEIN_FILE_PATH`
NB_CHUNKS=`expr $NB_PROTEIN / 100`

echo "Splitting protein file in $NB_CHUNKS chunks"
echo "fastasplit $PROTEIN_FILE_PATH $NB_CHUNKS ${PROTEIN_OUTPUT_DIR}/protein_chunks"

fastasplit $PROTEIN_FILE_PATH $NB_CHUNKS ${PROTEIN_OUTPUT_DIR}/protein_chunks

echo ""


# 5/ Make input_ids for the raw compute analysis in the pipeline

echo "Run make_input_ids"

cd ${ENS_PIPELINE_PATH}/scripts

# DNA ones

echo "perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitSlice -slice -coord_system $COORD_SYSTEM -slice_size 300000"

perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitSlice -slice -coord_system $COORD_SYSTEM -slice_size 300000
 perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitChromosome -slice -coord_system $COORD_SYSTEM
 perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name Submit30kSlice -slice -coord_system $COORD_SYSTEM -slice_size 30000

# Protein ones

echo "perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitProteome -single -single_name $PROTEIN_FILE"

perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitProteome -single -single_name $PROTEIN_FILE

echo "perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitTranslation -translation_ids"

perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitTranslation -translation_ids

echo "perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitChunk -file -dir ${PROTEIN_OUTPUT_DIR}/protein_chunks"

perl make_input_ids -dbuser $DB_USER -dbpass $DB_PASS -dbhost $DB_HOST -dbport $DB_PORT -dbname $DB_NAME -logic_name SubmitChunk -file -dir ${PROTEIN_OUTPUT_DIR}/protein_chunks

