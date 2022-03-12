#!/bin/bash

# CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS
if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 arguments";
    echo "    $0 EXPONENT LENGTH DEGREE_EQ and opt. arguments:";
    echo "    SIZE_ROOTS NUMBER_TESTS VERSION";
    echo "    .... EXPONENT: exp. of Kummer ext. L considered";
    echo "    .... LENGTH: length of sequence defining L, i.e. exp^length = [L:QQ]";
    echo "    .... DEGREE_EQ: degree of f(T)";
    echo "    .... SIZE_ROOTS: size of roots of f(T), default is 10";
    echo "    .... NUMBER_TESTS: number of tests for each field, default is 25";
    echo "    .... VERSION: shape of f(T) (split or single), default is split";
    exit;
fi


# mandatory parameters
EXPONENT=$1
LENGTH=$2
DEGREE_EQ=$3;

# DEFAULT parameters
DEFAULT_SIZE_ROOTS=10;
DEFAULT_NUMBER_TESTS=25;
DEFAULT_VERSION="split";

# assign variables to default values
SIZE_ROOTS=${4-${DEFAULT_SIZE_ROOTS}}
NUMBER_TESTS=${5-${DEFAULT_NUMBER_TESTS}}
VERSION=${6-${DEFAULT_VERSION}}

PARAMS=(EXPONENT LENGTH DEGREE_EQ
	SIZE_ROOTS NUMBER_TESTS VERSION)


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

# files - temp and logs
HEAD_FILE="${HEAD_DIR}/head_kummer_rel${STR_TAIL}";
CODE_FILE_GP="${HEAD_DIR}/kummer_gp${STR_TAIL}";
CODE_FILE_LLL_ABS="${HEAD_DIR}/kummer_lll_abs${STR_TAIL}";
CODE_FILE_LLL_REL="${HEAD_DIR}/kummer_lll_rel${STR_TAIL}";
LOG_FILE_GP="${LOGS_DIR}/kummer_gp${STR_TAIL}";
LOG_FILE_LLL_ABS="${LOGS_DIR}/kummer_lll_abs${STR_TAIL}";
LOG_FILE_LLL_REL="${LOGS_DIR}/kummer_lll_rel${STR_TAIL}";

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

# create code files
cat ${HEAD_FILE} "skel_kummer_absolute.c" > ${CODE_FILE_LLL_ABS};
cat ${HEAD_FILE} "skel_kummer_relative.c" > ${CODE_FILE_LLL_REL};
cat ${HEAD_FILE} "skel_kummer_gp.c" > ${CODE_FILE_GP};

# launch computations
gp ${CODE_FILE_LLL_ABS} 1>${LOG_FILE_LLL_ABS} 2>&1  &
gp ${CODE_FILE_LLL_REL} 1>${LOG_FILE_LLL_REL} 2>&1  &
gp ${CODE_FILE_GP} 1>${LOG_FILE_GP} 2>&1  &

exit 0;
