#!/bin/bash

######################################################################
# Safe Mode Parameters
######################################################################
set -e
set -u
set -o pipefail

######################################################################
# Set user-inputted arguments
######################################################################
ID=$1
OUTPUT_DIR=$2
TARGET_INDIVIDUAL_ID=$3
TARGET_INDIVIDUAL_OUTPUT_DIR=$4

######################################################################
# Making the directories to do the project work in
######################################################################

if [ ! -d ${OUTPUT_DIR}/${ID} ]
then
        mkdir -p ${OUTPUT_DIR}/${ID}
fi

if [ ! -d ${OUTPUT_DIR}/${ID}/wgs ]
then
        mkdir -p ${OUTPUT_DIR}/${ID}/wgs \
                 ${OUTPUT_DIR}/${ID}/stats
fi

######################################################################
# Making symlinks
######################################################################
# wgs directory
for TARGET_FILE in ${TARGET_INDIVIDUAL_OUTPUT_DIR}/${TARGET_INDIVIDUAL_ID}/wgs/${TARGET_INDIVIDUAL_ID}.*
do
ln -s \
$TARGET_FILE \
${OUTPUT_DIR}/${ID}/wgs/${ID}.$(realpath $TARGET_FILE | rev | cut -d"/" -f 1 | rev | cut -d"." -f1 --complement)
done

# wgs stats
for TARGET_FILE in ${TARGET_INDIVIDUAL_OUTPUT_DIR}/${TARGET_INDIVIDUAL_ID}/stats/${TARGET_INDIVIDUAL_ID}.wgs*
do
ln -s \
$TARGET_FILE \
${OUTPUT_DIR}/${ID}/stats/${ID}.$(realpath $TARGET_FILE | rev | cut -d"/" -f 1 | rev | cut -d"." -f1 --complement)
done

# Other files that don't start with ${TARGET_INDIVIDUAL_ID}.wgs as literal file name
for TARGET_FILE in ${TARGET_INDIVIDUAL_OUTPUT_DIR}/${TARGET_INDIVIDUAL_ID}/stats/*.{html,zip,json,log}
do
ln -s \
$TARGET_FILE \
${OUTPUT_DIR}/${ID}/stats/$(basename $TARGET_FILE)
done

######################################################################
# End message
######################################################################

echo "Done with $(basename ${OUTPUT_DIR})/${ID} Module 1b: Setup, Calling, & Phasing Symlinks"
