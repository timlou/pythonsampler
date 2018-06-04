#!/bin/bash

if [[ $# < 3 ]]; then
    echo "Usage: ./sampler.sh input.txt output.txt n_pts_generated"
    exit 0
fi

python3 PointsGenerator.py ${1} ${2} ${3}

