#!/bin/bash

cat $1 | while read structure
do
	rosetta_scripts.linuxgccrelease -s $structure -beta -ignore_zero_occupancy false -nstruct 20 -ex1 -ex2aro -parser:protocol  design_lcyc.xml
done
