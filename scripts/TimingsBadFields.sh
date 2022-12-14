#!/bin/bash

# DISCLAIMER : SHOULD NOT BE USED WITH SAME TYPE_FIELD AND 

if (( $# < 2 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 2 argument";
    echo " PropBadFields.sh SIZE_POL_FIELD DEGREE_EQ";
    echo "    with opt. arguments TYPE_FIELD  NUMBER_TESTS";
    echo "    ....SIZE_POL_FIELD: size of defining pol P_K(X)";
    echo "    ....DEGREE_EQ: degree of potential eq. f(T)";
    echo "    ....TYPE_FIELD: real or complex, default is real";
    echo "    ....NUMBER_TESTS: number of tests for each fields, default is 20";
    exit;
fi


# mandatory parameters
SIZE_POL_FIELD=$1
DEGREE_EQ=$2

# DEFAULT parameters
DEFAULT_TYPE_FIELD="real"
DEFAULT_NUMBER_TESTS=20

# assign variables to default values
TYPE_FIELD=${3-${DEFAULT_TYPE_FIELD}}
NUMBER_TESTS=${4-${DEFAULT_NUMBER_TESTS}}


PARAMS=(SIZE_POL_FIELD DEGREE_EQ
	TYPE_FIELD  NUMBER_TESTS);

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

HEAD_FILE="${HEAD_DIR}/head_timingsbadfields${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/timingsbadfields${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/timingsbadfields${STR_TAIL}";

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

cat ${HEAD_FILE} "skeletons/skel_timings_badfields.c" > ${CODE_FILE};

gp $CODE_FILE 1>$LOG_FILE 2>&1  &

exit 0;
