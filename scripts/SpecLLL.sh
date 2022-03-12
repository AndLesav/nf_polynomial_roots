#!/bin/bash

# DISCLAIMER : SHOULD NOT BE USED WITH SAME TYPE_FIELD AND 

if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 argument";
    echo " $0 DIM_START SIZE_POL_FIELD TYPE_FIELD";
    echo "    with opt. arguments NUMBER_DIM";
    echo "    .... DIM_START: smallest dim [K:\QQ] considered";
    echo "    .... SIZE_POL_FIELD: size of defining pol P_K(X)";
    echo "    .... TYPE_FIELD: real or complex";
    echo "    .... NUMBER_DIM: number of dim considered, default is 4 (50+3*50=200)";
    exit;
fi


# mandatory parameters
DIM_START=$1;
SIZE_POL_FIELD=$2;
TYPE_FIELD=$3;


# need diff type field for magma
if [ $TYPE_FIELD == "real" ]
then
    TYPE_FIELD_MAGMA=1;
elif [ $TYPE_FIELD == "complex" ]
then
    TYPE_FIELD_MAGMA=0;
fi

# DEFAULT parameters
DEFAULT_NUMBER_DIM=4;

# assign variables to default values
NUMBER_DIM=${4-${DEFAULT_NUMBER_DIM}}

PARAMS=(DIM_START SIZE_POL_FIELD # TYPE_FIELD will be handled separately
	NUMBER_DIM);

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

HEAD_FILE="${HEAD_DIR}/head_spec_lll${STR_TAIL}";
CODE_FILE="${HEAD_DIR}/spec_lll${STR_TAIL}";
LOG_FILE="${LOGS_DIR}/spec_lll${STR_TAIL}";

HEAD_FILE_GP="${HEAD_DIR}/head_spec_lll_gp${STR_TAIL}";
CODE_FILE_GP="${HEAD_DIR}/spec_lll_gp${STR_TAIL}";
LOG_FILE_GP="${LOGS_DIR}/spec_lll_gp${STR_TAIL}";

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
    echo "$PARAM = ${!PARAM};" >> ${HEAD_FILE_GP};
done

echo "TYPE_FIELD := ${TYPE_FIELD_MAGMA};" >> ${HEAD_FILE};
echo "TYPE_FIELD = ${TYPE_FIELD};" >> ${HEAD_FILE_GP};

cat ${HEAD_FILE} "skel_spec_lll.m" > ${CODE_FILE};
cat ${HEAD_FILE_GP} "skel_spec_lll.c" > ${CODE_FILE_GP};

magma $CODE_FILE 1>$LOG_FILE 2>&1  &
gp $CODE_FILE_GP 1>$LOG_FILE_GP 2>&1  &

exit 0;
