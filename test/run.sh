#!/bin/bash

FILEPATH="./MOLLUSC_Chr10_small.maf"
SP1="Octopusvulgaris6645"
SP2="Octopusbimaculoides37653"
SP3="Octopussinensis2607531"

python ../conservfinder.py -f ${FILEPATH} -t 0.7 -s "${SP1}" "${SP2}" "${SP3}" -o conserved_regions.rbh2
