#!/bin/bash

folder=`pwd`

$folder/../setup.sh $folder/smis_test $folder/ecoli_data/draft_assembly.fasta $folder/ecoli_data/ont.fastq &> /dev/null

cd $folder/smis_test
#./mysmissv.sh
