#!/bin/bash

timestamp() {
  date +"%Y-%m-%d_%H-%M-%S-%N"_$RANDOM
}

function format_fplll2magma {
    # reformat the output of fplll to suit the magma syntax
    # remove unnecessary spaces
    local FPLLL_OUT=$1;
    tr -s " " < ${FPLLL_OUT} > ${FPLLL_OUT}_1
    sed ${FPLLL_OUT}_1 -e "s/\[ /[/g" > ${FPLLL_OUT}
    sed ${FPLLL_OUT} -e "s/ \]/]/g" > ${FPLLL_OUT}_1
    #replace leftover space by commas
    tr " " "," < ${FPLLL_OUT}_1 > ${FPLLL_OUT}
    #separate in two parts to replace "] [" by "],[" except the end
    head -n-2 ${FPLLL_OUT} > ${FPLLL_OUT}_h
    tail -n2 ${FPLLL_OUT} > ${FPLLL_OUT}_t
    sed ${FPLLL_OUT}_h -e "s/\]$/\],/g" > ${FPLLL_OUT}_1
    #get the final result, remove intermediates
    cat ${FPLLL_OUT}_1 ${FPLLL_OUT}_t > ${FPLLL_OUT}
    rm -f ${FPLLL_OUT}_1 ${FPLLL_OUT}_h ${FPLLL_OUT}_t
}

function split_file {
    local INPUT_FILE=$1
    local RANGE=$2
    # echo ${RANGE}
    local OUTPUT_FILE_1=$3
    local OUTPUT_FILE_2=$4
    # echo "local variables created"
    head -n-${RANGE} ${INPUT_FILE}  > ${OUTPUT_FILE_1}
    tail -n${RANGE} ${INPUT_FILE} > ${OUTPUT_FILE_2}
    # echo "split finished"
}

#writes the load script
# echo "print \"\";" 1> ${LOAD_SCRIPT}
# echo "print \"Loading matrix from fplll\";" 1>> ${LOAD_SCRIPT}
function load_file_creation {
    local PRINT_FILE=$1
    local INPUT_FILE=$2
    local MATRIX_NAME=$3
    local LOAD_SCRIPT_TMP=LOAD_SCRIPT_TMP_$(timestamp)
    echo "${MATRIX_NAME}:=Matrix(" 1> ${LOAD_SCRIPT_TMP}
    cat ${LOAD_SCRIPT_TMP} ${INPUT_FILE}  1> ${PRINT_FILE}
    echo ");" 1>> ${PRINT_FILE}
    echo "return ${MATRIX_NAME};" 1>> ${PRINT_FILE}
    rm -f ${LOAD_SCRIPT_TMP}
}



# echo "usage : bash ${0} PRE_LLL_MAGMA POST_LLL_MAGMA MATRIX_NAME OUTPUT_NAME"
# echo "PRE_LLL_MAGMA: magma code which creates the matrix to reduce, named MATRIX_NAME"
# echo "MATRIX_NAME: name of the matrix to reduce in both magma codes"
# echo "POST_LLL_MAGMA: magma code to run after the matrix is reduced, must be terminated by \"exit;\""
# echo "OUTPUT_NAME: where the magma/fplll output/verbose is printed"

if [ "$#" -ne 6 ]; then
    echo "Illegal number of parameters: $# parameters given instead of 6"
    exit
fi

PRE_LLL_MAGMA=${1}
# POST_LLL_MAGMA=${2}
MATRIX_NAME=${2}
OUTPUT_NAME=${3}
PRINT_FILE=${4}
DIMENSION=${5}
DELTA=${6}

#commands to call magma and fplll, depends of your configuration
MAGMA_COMMAND="magma"
# FPLLL_COMMAND="fplll -a lll -v"
if (($DIMENSION==0)); then
    FPLLL_COMMAND="fplll -a lll -d ${DELTA}"
elif (($DIMENSION > 0)); then
    FPLLL_COMMAND="fplll -a lll -of bu -d ${DELTA}"
else
    echo "wrong argument for dimension"
    exit
fi

#name of the temporary magma programs
LOAD_SCRIPT=load_$(timestamp)
SAVE_SCRIPT=save_$(timestamp)

#name of the temporary matrices
FPLLL_ARG=${MATRIX_NAME}_$(timestamp)
FPLLL_OUT=${FPLLL_ARG}_reduced

#writes the save_script
echo "print \"\";" 1> ${SAVE_SCRIPT}
echo "print \"Printing matrix for fplll\";" 1>> ${SAVE_SCRIPT}
echo "PrintFile(\"${FPLLL_ARG}\",\"[\": Overwrite := true);" 1>> ${SAVE_SCRIPT}
echo "FILE:=Open(\"${FPLLL_ARG}\", \"a\");"  1>> ${SAVE_SCRIPT}
echo "fprintf FILE, \"%o\", ${MATRIX_NAME};"   1>> ${SAVE_SCRIPT};
echo "delete FILE;"  1>> ${SAVE_SCRIPT}
# echo "PrintFile(\"${FPLLL_ARG}\",${MATRIX_NAME});" 1>> ${SAVE_SCRIPT}
echo "PrintFile(\"${FPLLL_ARG}\",\"]\");" 1>> ${SAVE_SCRIPT}
echo "exit;" 1>> ${SAVE_SCRIPT}

#creates part1 of the magma run
PART1=magma_$(timestamp)
cat ${PRE_LLL_MAGMA} ${SAVE_SCRIPT} > ${PART1}
# cat $PART1;

# execute the first batch and save
# ${MAGMA_COMMAND} ${PRE_LLL_MAGMA} ${SAVE_SCRIPT} 1> ${OUTPUT_NAME}
${MAGMA_COMMAND} ${PART1} 1> ${OUTPUT_NAME}
rm -f ${SAVE_SCRIPT} ${PART1}
# cat ${OUTPUT_NAME}

# execute fplll
echo "" >> ${OUTPUT_NAME}
echo "Running ${FPLLL_COMMAND}:" >> ${OUTPUT_NAME}
${FPLLL_COMMAND} ${FPLLL_ARG} 1> ${FPLLL_OUT} 2>> ${OUTPUT_NAME}
echo "End of ${FPLLL_COMMAND} command" >> ${OUTPUT_NAME}
rm -f ${FPLLL_ARG}
# cat ${FPLLL_OUT}

# reformat the output of fplll to suit the magma syntax



#writes the load script
# echo "print \"\";" 1> ${LOAD_SCRIPT}
# echo "print \"Loading matrix from fplll\";" 1>> ${LOAD_SCRIPT}
# echo "${MATRIX_NAME}:=Matrix(" 1> ${LOAD_SCRIPT}_tmp
# cat ${LOAD_SCRIPT}_tmp ${FPLLL_OUT}  1>> ${LOAD_SCRIPT}
# echo ");" 1>> ${LOAD_SCRIPT}
# echo "return ${MATRIX_NAME};" 1>> ${LOAD_SCRIPT}
# rm -f ${LOAD_SCRIPT}_tmp ${FPLLL_OUT}

if (($DIMENSION==0)); then
    format_fplll2magma ${FPLLL_OUT}
    load_file_creation ${LOAD_SCRIPT} ${FPLLL_OUT} ${MATRIX_NAME}
    rm -f ${FPLLL_OUT}
    cat ${LOAD_SCRIPT} > ${PRINT_FILE}
elif (($DIMENSION > 0)); then
    OUTPUT1=${FPLLL_OUT}_1
    OUTPUT2=${FPLLL_OUT}_2
    # echo "before split"
    split_file ${FPLLL_OUT} $((${DIMENSION}+1)) ${OUTPUT1} ${OUTPUT2}
    # echo "end split"
    format_fplll2magma ${OUTPUT1}
    format_fplll2magma ${OUTPUT2}
    load_file_creation ${PRINT_FILE} ${OUTPUT1} ${MATRIX_NAME}
    load_file_creation ${PRINT_FILE}_uni ${OUTPUT2} ${MATRIX_NAME}_uni
    rm -f ${OUTPUT1} ${OUTPUT2} ${FPLLL_OUT}
else
    echo "wrong argument for dimension"
    exit
fi

# cat $PRINT_FILE
# #creates part2 of the magma run
PART2=magma_$(timestamp)
# cat ${LOAD_SCRIPT} ${POST_LLL_MAGMA} > ${PART2}

# # execute the load script and run the second part of the program
# ${MAGMA_COMMAND} ${PART2} 1>> ${OUTPUT_NAME}
rm -f ${LOAD_SCRIPT} ${PART2}
# rm -f LOAD_SCRIPT_TMP*
