#!/bin/bash

# DISCLAIMER : SHOULD NOT BE USED WITH SAME TYPE_FIELD AND 

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 argument";
    echo " Polfield_global.sh DIM SIZE_POL_FIELD DEGREE_EQ";
    echo "    with opt. arguments TYPE_FIELD NUMBER_TESTS TYPE_EQ";
    echo "    ....DIM: dim [K:\QQ] considered";
    echo "    ....SIZE_POL_FIELD: size of defining pol P_K(X)";
    echo "    ....DEGREE_EQ: degree of eq. f(T)";
    echo "    ....TYPE_FIELD: real or complex, default is real";
    echo "    ....NUMBER_TESTS: number of tests for each fields, default is 50";
    echo "    ....TYPE_EQ: split or not, default is split";
    exit;
fi


# mandatory parameters
DIM=$1
SIZE_POL_FIELD=$2
DEGREE_EQ=$3

# DEFAULT parameters
DEFAULT_TYPE_FIELD="real"
DEFAULT_NUMBER_TESTS=50;
DEFAULT_TYPE_EQ="split";

# assign variables to default values
TYPE_FIELD=${4-${DEFAULT_TYPE_FIELD}}
NUMBER_TESTS=${5-${DEFAULT_NUMBER_TESTS}}
TYPE_EQ=${6-${DEFAULT_TYPE_EQ}}

PARAMS=(DIM SIZE_POL_FIELD DEGREE_EQ
	TYPE_FIELD NUMBER_TESTS TYPE_EQ);

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

HEAD_FILE="${HEAD_DIR}/head_roots_absolute${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/sizeroots_absolute${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/sizeroots_absolute${STR_TAIL}";

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

cat ${HEAD_FILE} "skel_roots_absolute.c" > ${CODE_FILE};

gp ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
