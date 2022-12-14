#!/bin/bash


if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 argument";
    echo " Polfield_global.sh DIM TYPE_FIELD FRAC_ROOTS and opt. arguments:";
    echo "    SIZE_POL_FIELD SIZE_ROOTS NUMBER_TESTS TYPE_EQ";
    echo "    ....DIM: dim [K:\QQ] considered";
    echo "    ....TYPE_FIELD: real or complex";
    echo "    ....FRAC_ROOTS: deg(g)/deg(f) where f(T)=g(T)h(T) w. g split, Z(g)=Z(h)"; 
    echo "    ....SIZE_POL_FIELD: size of defining pol P_K(X), default is 1";
    echo "    ....SIZE_ROOTS: size of roots of f(T), default is 10";
    echo "    ....NUMBER_TESTS: number of tests for each fields, default is 50";
    echo "    ....TYPE_EQ: uniform or not, default is not (enforced anyway)";
    exit;
fi

# mandatory parameters
DIM=$1
TYPE_FIELD=$2;
FRAC_ROOTS=$3;

# DEFAULT parameters
DEFAULT_SIZE_POL_FIELD=2;
DEFAULT_SIZE_ROOTS=10;
DEFAULT_NUMBER_TESTS=50;
DEFAULT_TYPE_EQ="notuni";

# assign variables to default values
SIZE_POL_FIELD=${4-${DEFAULT_SIZE_POL_FIELD}}
SIZE_ROOTS=${5-${DEFAULT_SIZE_ROOTS}}
NUMBER_TESTS=${6-${DEFAULT_NUMBER_TESTS}}
TYPE_EQ=${7-${DEFAULT_TYPE_EQ}}

PARAMS=(DIM TYPE_FIELD FRAC_ROOTS
	SIZE_POL_FIELD SIZE_ROOTS NUMBER_TESTS TYPE_EQ);

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

HEAD_FILE="${HEAD_DIR}/head_nbsol_absolute${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/nbsol_absolute${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/nbsol_absolute${STR_TAIL}";

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

cat ${HEAD_FILE} "skeletons/skel_nbsol_absolute.c" > ${CODE_FILE};

gp ${CODE_FILE} 1>$LOG_FILE 2>&1  &

exit 0;
