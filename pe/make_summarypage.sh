#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <username> <output_path> <input_file>"
    exit 1
fi

# Assign input arguments to variables
USERNAME=$1
OUTPUT_PATH=$2
INPUT_FILE=$3

# Run the command with dynamic inputs
summarypages --v \
    --webdir "/home/${USERNAME}/public_html/${OUTPUT_PATH}" \
    --samples "${INPUT_FILE}" \
    --labels pTEOB \
    --gw \
    --no_ligo_skymap
