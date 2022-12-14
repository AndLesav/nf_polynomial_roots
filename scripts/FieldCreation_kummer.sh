#!/bin/bash

if (( $# < 2 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 2 arguments";
    echo "    $0 EXPONENT LENGTH and opt. arguments:";
    echo "    PRIMES_RANGE NUMBER_FIELDS";
    echo "    ....EXPONENT: exp. of Kummer ext. L considered";
    echo "    ....LENGTH: length of sequence defining L, i.e. exp^length = [L:QQ]";
    echo "    ....PRIMES_RANGE: range for primes pi defining L, default is 40";
    echo "    ....NUMBER_FIELDS: number of tests, default is 25";
    exit;
fi

# mandatory parameters
EXPONENT=$1
LENGTH=$2

# DEFAULT parameters
DEFAULT_PRIMES_RANGE=40;
DEFAULT_NUMBER_FIELDS=25;

# assign variables to default values
PRIMES_RANGE=${3-${DEFAULT_PRIMES_RANGE}}
NUMBER_FIELDS=${4-${DEFAULT_NUMBER_FIELDS}}

PARAMS=(EXPONENT LENGTH
	PRIMES_RANGE NUMBER_FIELDS);

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

HEAD_FILE="${HEAD_DIR}/head_field_creation_kummer${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/field_creation_kummer${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/field_creation_kummer${STR_TAIL}";

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

cat ${HEAD_FILE} "skeletons/field_creation_kummer.m" > ${CODE_FILE};

magma ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
