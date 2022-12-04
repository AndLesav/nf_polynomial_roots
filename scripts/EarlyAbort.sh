#!/bin/bash

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 argument";
    echo " $0 DIM_START SIZE_POL_FIELD TYPE_FIELD";
    echo "    with opt. arguments NUMBER_DIM NUMBER_TESTS";
    echo "    .... DIM_START: smallest dim [K:\QQ] considered";
    echo "    .... SIZE_POL_FIELD: size of defining pol P_K(X)";
    echo "    .... TYPE_FIELD: real or complex";
    echo "    .... NUMBER_DIM: number of dim considered, default is 5 (25+4*25=150)";
    echo "    .... NUMBER_TESTS: number of tests for each fields, default is 100";
    exit;
fi


# mandatory parameters
DIM_START=$1;
SIZE_POL_FIELD=$2;
TYPE_FIELD=$3;

# DEFAULT parameters
DEFAULT_NUMBER_DIM=5;
DEFAULT_NUMBER_TESTS=100;

# assign variables to default values
NUMBER_DIM=${4-${DEFAULT_NUMBER_DIM}}
NUMBER_TESTS=${5-${DEFAULT_NUMBER_TESTS}}

PARAMS=(DIM_START SIZE_POL_FIELD TYPE_FIELD
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

HEAD_FILE="${HEAD_DIR}/head_early_abort${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/early_abort${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/early_abort${STR_TAIL}";

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
    echo "$PARAM = ${!PARAM};" >> ${HEAD_FILE};
done

cat ${HEAD_FILE} "skel_early_abort.c" > ${CODE_FILE};

gp $CODE_FILE 1>$LOG_FILE 2>&1  &

exit 0;
