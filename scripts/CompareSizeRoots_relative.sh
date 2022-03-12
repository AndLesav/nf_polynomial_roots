#!/bin/bash

# DISCLAIMER : SHOULD NOT BE USED WITH SAME TYPE_FIELD AND 

if (( $# < 2 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 2 arguments";
    echo " Polfield_global.sh DIM_G DIM_E and opt. arguments:";
    echo "    DEGREE_EQ NUMBER_TESTS";
    echo "    ....DIM_G: dim of ground field [K:\QQ] considered";
    echo "    ....DIM_E: dim of relative ext [L:K] considered";
    echo "    ....DEGREE_EQ: degree of f(T), default is 25";
    echo "    ....NUMBER_TESTS: number of tests for each fields, default is 50";
    exit;
fi

# mandatory parameters
DIM_G=$1
DIM_E=$2

# DEFAULT parameters
DEFAULT_DEGREE_EQ=25;
DEFAULT_NUMBER_TESTS=50;

# assign variables to default values
DEGREE_EQ=${3-${DEFAULT_DEGREE_EQ}}
NUMBER_TESTS=${4-${DEFAULT_NUMBER_TESTS}}

PARAMS=(DIM_G DIM_E
	DEGREE_EQ NUMBER_TESTS);

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

HEAD_FILE="${HEAD_DIR}/head_compar_sizeroots_relative${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/compar_sizeroots_relative${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/compar_sizeroots_relative${STR_TAIL}";

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

cat ${HEAD_FILE} "skel_comp_sizeroots_relative.c" > ${CODE_FILE};

gp ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
