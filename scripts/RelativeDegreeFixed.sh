#!/bin/bash

# CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS
if (( $# < 3 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 3 arguments";
    echo "    $0 DIM_EXT DIM_G_m DIM_G_M and opt. arguments:";
    echo "    DEGREE_EQ NUMBER_TESTS VERSION";
    echo "    ....DIM_EXT: ext. dim. of relativ ext [L:K] considered";
    echo "    ....DIM_G_m: min. dim. of ground field [K:QQ] considered";
    echo "    ....DIM_G_M: max. dim. of ground field [K:QQ] considered";
    echo "    .... DEGREE_EQ: degree of f(T), default is 10";
    echo "    .... NUMBER_TESTS: number of tests for each field, default is 25";
    echo "    .... SIZE_ROOTS: max size of roots, default is 100";
    echo "    .... VERSION: shape of f(T) (split or single), default is single";
    exit;
fi


# mandatory parameters
DIM_EXT=$1
DIM_G_m=$2
DIM_G_M=$3

# DEFAULT parameters
DEFAULT_DEGREE_EQ=10
DEFAULT_NUMBER_TESTS=25;
DEFAULT_SIZE_ROOTS=100;
DEFAULT_VERSION="single";

# assign variables to default values
DEGREE_EQ=${4-${DEFAULT_DEGREE_EQ}}
NUMBER_TESTS=${5-${DEFAULT_NUMBER_TESTS}}
SIZE_ROOTS=${6-${DEFAULT_SIZE_ROOTS}}
VERSION=${7-${DEFAULT_VERSION}}


# need diff VERSION for magma
if [ $VERSION == "single" ]
then
    VERSION_MAGMA=0;
elif [ $VERSION == "split" ]
then
    VERSION_MAGMA=1;
fi

PARAMS=(DIM_EXT DIM_G_m DIM_G_M
	DEGREE_EQ NUMBER_TESTS SIZE_ROOTS) # VERSION will be handled separately


STR_TAIL=""
for PARAM in ${PARAMS[@]}; do
    STR_TAIL="${STR_TAIL}_${!PARAM}";
done
STR_TAIL="${STR_TAIL}_${VERSION}";

# Script folder
EXE_DIR=$(dirname $(readlink -f $0));
ROOT_DIR=$(dirname ${EXE_DIR});
DATA_DIR="${ROOT_DIR}/data";
LOGS_DIR="${ROOT_DIR}/logs";
HEAD_DIR="${ROOT_DIR}/heads";

# files - temp and logs
HEAD_FILE_MAGMA="${HEAD_DIR}/head_reldeg_fixed_mgm${STR_TAIL}";
HEAD_FILE="${HEAD_DIR}/head_reldeg_fixed_gp${STR_TAIL}";

CODE_FILE_MAGMA="${HEAD_DIR}/reldeg_fixed_mgm${STR_TAIL}";
CODE_FILE_GP="${HEAD_DIR}/reldeg_fixed_gp${STR_TAIL}";
CODE_FILE_LLL_ABS="${HEAD_DIR}/reldeg_fixed_lll_abs${STR_TAIL}";
CODE_FILE_LLL_REL="${HEAD_DIR}/reldeg_fixed_lll_rel${STR_TAIL}";

LOG_FILE_MAGMA="${LOGS_DIR}/reldeg_fixed_mgm${STR_TAIL}";
LOG_FILE_GP="${LOGS_DIR}/reldeg_fixed_gp${STR_TAIL}";
LOG_FILE_LLL_ABS="${LOGS_DIR}/reldeg_fixed_lll_abs${STR_TAIL}";
LOG_FILE_LLL_REL="${LOGS_DIR}/reldeg_fixed_lll_rel${STR_TAIL}";

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
    echo "$PARAM := ${!PARAM};" >> ${HEAD_FILE_MAGMA};
done

echo "VERSION = ${VERSION};" >> ${HEAD_FILE};
echo "VERSION := ${VERSION_MAGMA};" >> ${HEAD_FILE_MAGMA};


# create code files
cat ${HEAD_FILE} "skeletons/skel_reldegree_fixed_absolute.c" > ${CODE_FILE_LLL_ABS};
cat ${HEAD_FILE} "skeletons/skel_reldegree_fixed_relative.c" > ${CODE_FILE_LLL_REL};
cat ${HEAD_FILE} "skeletons/skel_reldegree_fixed_gp.c" > ${CODE_FILE_GP};
cat ${HEAD_FILE_MAGMA} "skeletons/skel_reldegree_fixed_mgm.m" > ${CODE_FILE_MAGMA};

# launch computations
gp ${CODE_FILE_LLL_ABS} 1>${LOG_FILE_LLL_ABS} 2>&1  &
gp ${CODE_FILE_LLL_REL} 1>${LOG_FILE_LLL_REL} 2>&1  &
gp ${CODE_FILE_GP} 1>${LOG_FILE_GP} 2>&1  &
magma ${CODE_FILE_MAGMA} 1>${LOG_FILE_MAGMA} 2>&1 &

exit 0;
