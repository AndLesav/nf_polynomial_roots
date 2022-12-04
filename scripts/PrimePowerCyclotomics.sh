#!/bin/bash

# DISCLAIMER : SHOULD NOT BE USED WITH SAME TYPE_FIELD AND 

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 arguments";
    echo " $0 PRIME DEGREE_EQ SIZE_ROOTS";
    echo "    with opt. arguments NUMBER_TESTS DIM_MAX";
    echo "    .... PRIME: prime p defining conductor of L";
    echo "    .... DEGREE_EQ: degree of eq. f(T)";
    echo "    .... SIZE_ROOTS: size of roots of f(T)";
    echo "    .... NUMBER_TESTS: number of tests for each fields, default is 15";
    echo "    .... DIM_MAX: largest dim [L:\QQ] considered, default is 256";
    exit;
fi


# mandatory parameters
PRIME=$1
DEGREE_EQ=$2
SIZE_ROOTS=$3

# DEFAULT parameters
DEFAULT_NUMBER_TESTS=15;
DEFAULT_DIM_MAX=250;

# assign variables to default values
NUMBER_TESTS=${4-${DEFAULT_NUMBER_TESTS}}
DIM_MAX=${5-${DEFAULT_DIM_MAX}}

PARAMS=(PRIME DEGREE_EQ SIZE_ROOTS
	 NUMBER_TESTS DIM_MAX);

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

HEAD_FILE="${HEAD_DIR}/head_primepower_cyclotomics${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/primepower_cyclotomics${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/primepower_cyclotomics${STR_TAIL}";

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

cat ${HEAD_FILE} "skel_primepower_cyclotomics.c" > ${CODE_FILE};

gp $CODE_FILE 1>$LOG_FILE 2>&1  &

exit 0;
