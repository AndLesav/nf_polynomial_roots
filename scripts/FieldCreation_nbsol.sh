#!/bin/bash

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 arguments";
    echo "    $0 DIM_ABS DIM_EXT_m DIM_EXT_M and opt. argument:";
    echo "    NUMBER_FIELDS";
    echo "    ....DIM_A: abs. dim. of field [L:\QQ] considered";
    echo "    ....DIM_E_m: min. dim. of relative ext [L:K] considered";
    echo "    ....DIM_E_M: max. dim. of relative ext [L:K] considered";
    echo "    ....NUMBER_FIELDS: number of fields to be created for each set of param., default is 50";
    exit;
fi

# mandatory parameters
DIM_ABS=$1
DIM_EXT_m=$2
DIM_EXT_M=$3

# DEFAULT parameters
DEFAULT_NUMBER_FIELDS=50;

# assign variables to default values
NUMBER_FIELDS=${4-${DEFAULT_NUMBER_FIELDS}}

PARAMS=(DIM_ABS DIM_EXT_m DIM_EXT_M
	 NUMBER_FIELDS);

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

HEAD_FILE="${HEAD_DIR}/head_field_creation_nbsol${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/field_creation_nbsol${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/field_creation_nbsol${STR_TAIL}";

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

cat ${HEAD_FILE} "field_creation_nbsol.m" > ${CODE_FILE};

magma ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
