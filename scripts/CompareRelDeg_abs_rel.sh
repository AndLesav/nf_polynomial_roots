#!/bin/bash


# CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS 
if (( $# < 1 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 1 argument";
    echo "    $0 DIM_G and opt. arguments:";
    echo "    DEGREE_EQ SIZE_ROOTS NUMBER_TESTS";
    echo "    ....DIM_G: dim. of ground field [K:\QQ] considered";
    echo "    ....DEGREE_EQ: degree of f(T), default is 50";
    echo "    ....SIZE_ROOTS: size of roots of f(T), default is 10";
    echo "    ....NUMBER_TESTS: number of tests for each field, default is 50";
    exit;
fi


# mandatory parameters
DIM_G=$1

# DEFAULT parameters
DEFAULT_DEGREE_EQ=50;
DEFAULT_SIZE_ROOTS=10;
DEFAULT_NUMBER_TESTS=50;

# assign variables to default values
DEGREE_EQ=${2-${DEFAULT_DEGREE_EQ}}
SIZE_ROOTS=${3-${DEFAULT_SIZE_ROOTS}}
NUMBER_TESTS=${4-${DEFAULT_NUMBER_TESTS}}

PARAMS=(DIM_G
	DEGREE_EQ SIZE_ROOTS NUMBER_TESTS);

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

HEAD_FILE="${HEAD_DIR}/head_compar_reldeg_abs_rel${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/compar_reldeg_abs_rel${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/compar_reldeg_abs_rel${STR_TAIL}";

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

cat ${HEAD_FILE} "skel_reldeg_abs_rel.c" > ${CODE_FILE};

gp ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
