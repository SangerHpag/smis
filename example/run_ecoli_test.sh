#!/bin/bash

folder=`pwd`

$folder/../setup.sh $folder/smis_test $folder/ecoli_data/draft_assembly.fasta $folder/ecoli_data/ont.fastq &> exmaple.log

cd $folder/smis_test
./mysmissv.sh  &>> example.log


echo; echo Scaffolds are in  smis_test/spinner_scaffolds.fasta
echo example log output: exmaple.log
echo smis logs in smis_test/logs/
echo
