#!/bin/bash

folder=`pwd`

$folder/../setup.sh $folder/smis_test $folder/ecoli_data/draft_assembly.fasta $folder/ecoli_data/ont.fastq &> log.out

cd $folder/smis_test
./mysmissv.sh  &>> log.out


echo Scaffolds are in  smis_test/spinner_scaffolds.fasta
