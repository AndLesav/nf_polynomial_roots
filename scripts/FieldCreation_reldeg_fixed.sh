#!/bin/bash

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 arguments";
    echo "    $0 DIM_EXT DIM_G_m DIM_G_M and opt. argument:";
    echo "    NUMBER_FIELDS DEGREE_EQ";
    echo "    ....DIM_EXT: ext. dim. of relativ ext [L:K] considered";
    echo "    ....DIM_G_m: min. dim. of ground field [K:QQ] considered";
    echo "    ....DIM_G_M: max. dim. of ground field [K:QQ] considered";
    echo "    ....NUMBER_FIELDS: number of fields to be created for each set of
    	      			 param., default is 50";	 
    echo "    ....DEGREE_EQ: degree of potential equation f(T), default is 10";
	
    exit;
fi

# mandatory parameters
DIM_EXT=$1
DIM_G_m=$2
DIM_G_M=$3

# DEFAULT parameters
DEFAULT_NUMBER_FIELDS=50;
DEFAULT_DEGREE_EQ=10;

# assign variables to default values
NUMBER_FIELDS=${4-${DEFAULT_NUMBER_FIELDS}}
DEGREE_EQ=${5-${DEFAULT_DEGREE_EQ}}


PARAMS=(DIM_EXT DIM_G_m DIM_G_M
	NUMBER_FIELDS DEGREE_EQ);

STR_TAIL=""
for PARAM in ${PARAMS[@]}; do
    STR_TAIL="${STR_TAIL}_${!PARAM}";
done

# Script folder
EXE_DIR=$(dirname $(readlink -f $0));
ROOT_DIR=$(dirname ${EXE_DIR});
DATA_DIR="${ROOT_DIR}/data";
LOGS_DIR="${ROOT_DIR}/logs";
HEAD_DIR="${ROOT_DIR}/heads";

HEAD_FILE="${HEAD_DIR}/head_field_creation_reldeg_fixed${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/field_creation_reldeg_fixed${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/field_creation_reldeg_fixed${STR_TAIL}";

# Just check that parent folders are indeed where they should be
[[ ! -d ${DATA_DIR} ]] && {
    echo -e "\x1b[31m[Err]\x1b[0m Data directory ${DATA_DIR} does not exist.";
    exit 1;
};

[[ ! -d ${LOGS_DIR} ]] && {
    echo -e "\x1b[31m[Err]\x1b[0m Logs directory ${LOGS_DIR} does not exist.";
    exit 1;
};

for PARAM in ${PARAMS[@]}; do
    echo "$PARAM := ${!PARAM};" >> ${HEAD_FILE};
done

cat ${HEAD_FILE} "field_creation_reldeg_fixed.m" > ${CODE_FILE};

magma ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
