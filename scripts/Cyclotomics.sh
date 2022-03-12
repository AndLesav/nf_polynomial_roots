#!/bin/bash

# DISCLAIMER : SHOULD NOT BE USED WITH SAME TYPE_FIELD AND 

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 arguments";
    echo " $0 DIM_MAX DEGREE_EQ SIZE_ROOTS";
    echo "    with opt. arguments NUMBER_TESTS CONDS_MIN CONDS_MAX";
    echo "    ....DIM_MAX: largest dim [K:\QQ] considered";
    echo "    ....DEGREE_EQ: degree of eq. f(T)";
    echo "    ....SIZE_ROOTS: size of roots of f(T)";
    echo "    ....NUMBER_TESTS: number of tests for each fields, default is 50";
    echo "    ....CONDS_MIN: smallest conductor of K considered, default is 20";
    echo "    ....CONDS_MAX: largest conductor of K considered, default is 500";
    exit;
fi


# mandatory parameters
DIM_MAX=$1
DEGREE_EQ=$2
SIZE_ROOTS=$3

# DEFAULT parameters
DEFAULT_NUMBER_TESTS=50;
DEFAULT_CONDS_MIN=20;
DEFAULT_CONDS_MIN=500;

# assign variables to default values
NUMBER_TESTS=${4-${DEFAULT_NUMBER_TESTS}}
CONDS_MIN=${5-${DEFAULT_CONDS_MIN}}
CONDS_MAX=${6-${DEFAULT_CONDS_MAX}}

PARAMS=(DIM_MAX DEGREE_EQ SIZE_ROOTS
	 NUMBER_TESTS CONDS_MIN CONDS_MAX);

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

HEAD_FILE="${HEAD_DIR}/head_cyclotomics${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/cyclotomics${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/cyclotomics${STR_TAIL}";

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

cat ${HEAD_FILE} "skel_cyclotomics.c" > ${CODE_FILE};

gp $CODE_FILE 1>$LOG_FILE 2>&1  &

exit 0;
