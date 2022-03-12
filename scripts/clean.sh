#!/bin/bash

if (( $# > 1 )); then
    echo -e  "\033[31merror in $0:\033[0m need at most 1 argument";
    echo " clean.sh DIR"
    echo "    ....DIR: logs heads or data";
    exit;
elif (( $# < 1 )); then
    echo -e  "\033[31merror in $0:\033[0m need at least 1 argument";
    echo " clean.sh DIR"
    echo "    ....DIR: logs heads or data";
    exit;
elif [ $1 != "logs" ] && [ $1 != "heads" ] && [ $1 != "data" ] && [ $1 != "inputs" ];
then
    echo -e  "\033[31merror in $0:\033[0m argument given is uncorrect";
    echo " clean.sh DIR"
    echo "    ....DIR: logs heads data or inputs";
    exit;
fi

DIR=$1;

EXE_DIR=$(dirname $(readlink -f $0));
ROOT_DIR=$(dirname ${EXE_DIR});
STR="${ROOT_DIR}/${DIR}/*";
rm $STR

exit 0;
