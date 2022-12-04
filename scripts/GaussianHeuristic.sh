#!/bin/bash

if (( $# < 1 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 1 argument";
    echo "    $0 DIM_m TYPE_FIELD and opt. argument:";
    echo "    NUMBER_DIM NUMBER_TESTS";
    echo "    .... DIM_m: min. dim. of field [K:\QQ] considered";
    echo "    .... TYPE_FIELD: type of field, 1 for real (r1>0) or 0 for complex (r1=0)";
    echo "    .... NUMBER_DIM: number of dim. to be computed, default is 8 (thought to start from DIM_G = 15)";
    echo "    .... NUMBER_TESTS: number of fields to be tested for each set of param., default is 50";
    exit;
fi

# mandatory parameters
DIM_m=$1
TYPE_FIELD=$2;

# DEFAULT parameters
DEFAULT_NUMBER_DIM=8;
DEFAULT_NUMBER_TESTS=50;

# assign variables to default values
NUMBER_DIM=${3-${DEFAULT_NUMBER_DIM}}
NUMBER_TESTS=${4-${DEFAULT_NUMBER_TESTS}}

PARAMS=(DIM_m TYPE_FIELD
	NUMBER_DIM NUMBER_TESTS);

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

HEAD_FILE="${HEAD_DIR}/head_gaussian_heuristic${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/gaussian_heuristic${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/gaussian_heuristic${STR_TAIL}";

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

cat ${HEAD_FILE} "gaussian_heuristic.m" > ${CODE_FILE};

magma ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
