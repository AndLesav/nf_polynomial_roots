#!/bin/bash

# DISCLAIMER : SHOULD NOT BE USED WITH SAME TYPE_FIELD AND 

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 argument";
    echo " $0 DIM_START SIZE_POL_FIELD TYPE_FIELD";
    echo "    with opt. arguments NUMBER_DIM DEGREE_EQ  NUMBER_TESTS";
    echo "    .... DIM_START: smallest dim [K:\QQ] considered";
    echo "    .... SIZE_POL_FIELD: size of defining pol P_K(X)";
    echo "    .... TYPE_FIELD: real or complex";
    echo "    .... NUMBER_DIM: number of dim considered, default is 20 (15+19*5=110)";
    echo "    .... DEGREE_EQ: degree of eq. f(T), default is 25";
    echo "    .... NUMBER_TESTS: number of tests for each fields, default is 100";
    exit;
fi


# mandatory parameters
DIM_START=$1;
SIZE_POL_FIELD=$2;
TYPE_FIELD=$3;

# DEFAULT parameters
DEFAULT_NUMBER_DIM=20;
DEFAULT_DEGREE_EQ=25;
DEFAULT_NUMBER_TESTS=100;

# assign variables to default values
NUMBER_DIM=${4-${DEFAULT_NUMBER_DIM}}
DEGREE_EQ=${5-${DEFAULT_DEGREE_EQ}}
NUMBER_TESTS=${6-${DEFAULT_NUMBER_TESTS}}

PARAMS=(DIM_START SIZE_POL_FIELD TYPE_FIELD
	NUMBER_DIM DEGREE_EQ NUMBER_TESTS);

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

HEAD_FILE="${HEAD_DIR}/head_norm_evaluation${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/norm_evaluation${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/norm_evaluation${STR_TAIL}";

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

cat ${HEAD_FILE} "skeletons/skel_norm_eval.c" > ${CODE_FILE};

gp $CODE_FILE 1>$LOG_FILE 2>&1  &

exit 0;
