#!/bin/sh

if [ $# != 1 ]
then
    echo "Wrong number of command line arguments"
    echo "sh run_xrefs_pipeline_generic.sh saccharomyces_cerevisiae_core_10_63_3"
    exit 1
fi

source /homes/oracle/ora920setup.sh

CORE_DB_NAME=$1

XREF_DB_NAME=`echo "$CORE_DB_NAME" | sed s/core/xref/`

echo "XREF_DB_NAME: $XREF_DB_NAME"

# specific env variables in a user specific file now

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

# Test the CONFIG_FILE

if [ ! -d ${CONFIG_DIR} ] 
then
    echo "CONFIG_DIR is not a valid directory, ${CONFIG_DIR}"
    echo "check you have specified it correctly if your ${USER}.env.sh environment file"
    exit 1
fi

echo "Config Dir: $CONFIG_DIR"

#DB_HOST=mysql-cluster-eg-prod-1.ebi.ac.uk
#DB_PORT=4238
#DB_USER=ensrw
#DB_PASS=writ3rp1

SPECIES=`echo $CORE_DB_NAME | perl -ne '$_ =~ /^([^_]+)_([^_]+)_core.+/; $a = $1; $b = $2; print $a . "_" . $b;'`
echo "SPECIES: $SPECIES"

# e.g. anidulans
SPECIES_SHORT_NAME=`echo $SPECIES | perl -ne '$_ =~ /^(\w)[^_]+_(\w+)/; $a = $1; $b = $2; print "$a$b";'`

echo "SPECIES_SHORT_NAME: $SPECIES_SHORT_NAME"

PERL_PATH=/nfs/panda/ensemblgenomes/perl/
ENSEMBL_ROOT_DIR=/nfs/panda/ensemblgenomes

BIOPERL_PATH=/nfs/panda/ensemblgenomes/apis/bioperl/stable/
DATA_DOWNLOAD_DIR=/nfs/nobackup/ensemblgenomes/production/xrefs_pipeline/data/${SPECIES_SHORT_NAME}/input_data
# Right now, data can't go to nobackup (excepts download) because of nfs issue
DATA_OUTPUT_DIR=/nfs/panda/ensemblgenomes/production/xrefs_pipeline/data/${SPECIES_SHORT_NAME}

if [ ! -d "${DATA_DOWNLOAD_DIR}" ] 
then
    echo "Creating directory structure, ${DATA_DOWNLOAD_DIR}"
    mkdir -p ${DATA_DOWNLOAD_DIR}
fi

if [ ! -d "${DATA_OUTPUT_DIR}" ] 
then
    echo "reating directory structure, ${DATA_OUTPUT_DIR}"
    mkdir -p ${DATA_OUTPUT_DIR}
fi

export PERL5LIB=${ENSEMBL_ROOT_DIR}/production/xrefs_pipeline/ensembl-head/modules:${ENSEMBL_ROOT_DIR}/apis/ensembl/pipeline/head/modules:${BIOPERL_PATH}

# Add ensgen perl binary path
export PATH=${PERL_PATH}/perlbrew/perls/5.14.2/bin:$PATH

cd ${ENSEMBL_ROOT_DIR}/production/xrefs_pipeline/ensembl-head/misc-scripts/xref_mapping

# Parsing stage

echo ""
echo "Running xref_parser.pl"

echo "y" | perl xref_parser.pl -user $DB_USER -pass $DB_PASS -host $DB_HOST -port $DB_PORT -species $SPECIES -create -dbname $XREF_DB_NAME -checkdownload -download_dir ${DATA_DOWNLOAD_DIR} -drop_db

# Mapping stage

if [ ! -f ${CONFIG_DIR}"/"${SPECIES_SHORT_NAME}"_xref_mapper.input" ]
then
    echo "Config file, "${CONFIG_DIR}"/"${SPECIES_SHORT_NAME}"_xref_mapper.input, not found!"
    exit 1
fi

echo ""
echo "Running xref_mapper.pl"
echo "perl xref_mapper.pl -file "${CONFIG_DIR}"/"${SPECIES_SHORT_NAME}"_xref_mapper.input -upload"

# 55 and further

perl xref_mapper.pl -file ${CONFIG_DIR}"/"${SPECIES_SHORT_NAME}"_xref_mapper.input" -upload

