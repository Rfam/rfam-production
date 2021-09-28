#!/bin/bash
FILENAME=$1
if [ -f ${FILENAME} ]
then
    if [ -s ${FILENAME} ]
    then
        echo "File exists and not empty"
        exit 0
    else
        echo "File exists but empty"
        exit 1
    fi
else
    echo "File does not exist"
    exit 1
fi