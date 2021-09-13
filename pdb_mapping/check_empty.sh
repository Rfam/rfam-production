#!/bin/bash

FILENAME="relX_text_search/families/error.log"

if [ -f ${FILENAME} ]
then
    if [ -s ${FILENAME} ]
    then
        echo "File exists and not empty"
        exit 1
    else
        echo "File exists but empty"
        exit 0
    fi
else
    echo "File does not exist"
    exit 1
fi