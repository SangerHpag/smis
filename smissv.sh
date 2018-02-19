#!/bin/bash
#set -o errexit
#set -o pipefail



#settings file
sfile=ODIR/mysettings.sh 

source $sfile

######################### A few Checks: ##########################
#mybwa if aligner, samtools,


### check if aligner defined and located:
echo "${aligner:?Need to choose an aligner in mysettings.sh}" > /dev/null
if [ $aligner == "bwa" ]; then
    echo "${mybwa:?bwa is not in your PATH, please set mybwa variable to bwa location in mysettings.sh}" > /dev/null
else
    echo "${mysmalt:?smalt is not in your PATH, please set mysmalt variable to smalt location mysettings.sh}" > /dev/null
fi

## check is samtools found
echo "${mysamtools:?samtools is not in your PATH, please set mysamtools variable to samtools location in mysettings.sh}" > /dev/null


######################### Check if input files exist: ##########################
if [ ! -f $fafile ]; then
  echo '   !!!!!!!!!!  Error: fasta file $fafile not found or not readable!' 
  echo '                 Please define the fafile variable in ./mysettings.sh file !' 
  exit
else  
  check=`head -1 $fafile`
    if [ -z "${check}" ]; then
        echo; echo " Error! file " $fafile " appears to be empty"
        exit
    fi
fi

if [ ! -f $fqfile ]; then
  echo '   !!!!!!!!!!  Error: fastq file $fqfile not found or not readable!'
  echo '                 Please define the fqfile variable in ./mysettings.sh file !' 
  exit
else  
  check=`head -1 $fqfile`
    if [ -z "${check}" ]; then
        echo; echo " Error! file " $fqfile " appears to be empty"
        exit
    fi
fi


mkdir -p $odir/logs
thisrun=`date +%s`
las=$odir/logs/launchedas_$thisrun.txt
outp=$odir/logs/output_$thisrun.txt


workdir=$odir/tempWork
mkdir -p $workdir
cd $workdir

############################ Save settings #####################################

echo '################## Settings Used for Running ##################' > $las 
cat $sfile >> $las


echo '################## Running Output ##################' > $outp

rerun=no

### rerun everything or keep existing??
if [ -f "fakemates_1.fastq" ];  then
  echo
  echo "Warning: It appears you already run the pipeline in this folder: do you want to recreate everything from scratch? [ENTER] "
  read rerun
  
  if [ $rerun == "y" ] || [ $rerun == "yes" ]; then
   echo -e "\nRerunning from scratch \n" >> $outp
  else
    echo -e "\nKeeping existing files from previous run \n" >> $outp
  fi
fi



if [ $rerun == "y" ] || [ $rerun == "yes" ] || [ ! -f "fakemates_1.fastq" ] || [ ! -f "fakemates_2.fastq" ] ; then 
  echo -e "\nCreating fake-mate pairs from long read fastq file \n" >> $outp
  $bindir/smis_shred -rlength $fakelen -step $step -minlen $minlen $fqfile fakemates_1.fastq fakemates_2.fastq >> $outp 
else
 echo -e "\nFake-mate pairs exist already \n" >> $outp  
fi



if [ $rerun == "y" ] || [ $rerun == "yes" ] || [ ! -f "genome.fasta" ] ; then 
 echo -e "\nPreparing Reference Fasta file\n" >> $outp 
 $bindir/smis_rename $fafile genome.fasta
else
 echo -e "\nReference Fasta file prepared already \n" >> $outp             
fi


if [ $aligner == "smalt" ]; then
  if [ -z $mysmalt ]; then 
      echo "Cannot find smalt: please define \"mysmalt\" in the settings file \@ line 33" 
      exit 
  fi

  if [ $rerun == "y" ] || [ $rerun == "yes" ] || [ ! -f "hash_target.smi" ] ; then  
    echo -e "\nAligning Fake-Mate Pairs to References with SMALT\n"  >> $outp 

    $mysmalt index -k 15 -s 6 hash_target genome.fasta 
    smaltpars="match=$match,subst=$subst,gapopen=$gapopen,gapext=$gapext"
    $mysmalt map -m $mapscore -n $nodes -O -f ssaha -o align.dat -i 40000 -j 1000 -S $smaltpars hash_target fakemates_1.fastq fakemates_2.fastq > /dev/null
  fi

 if [ $rerun == "y" ] || [ $rerun == "yes" ] || [ ! -f "mates_id60_mscore10.out" ] ; then
   echo -e "\n
     Sort Aligned Fake-Mate Pairs by their Position in Reference\n"  >> $outp
   $bindir/smis_sort align.dat
 fi 
elif  [ $aligner == "bwa" ]; then  
  if [ -z $mybwa ]; then
      echo "Cannot find smalt: please define \"mybwa\" in the settings file \@ line 33" 
      exit
  fi

  if [ $rerun == "y" ] || [ $rerun == "yes" ] || [ ! -f "bwa_sorted.bam" ] ; then
    echo -e "\nAligning Fake-Mate Pairs to References with BWA\n"  >> $outp
    
    $mybwa index genome.fasta &>> $outp 
    $mybwa mem -t $nodes -T $mapscore  -A $match -O $gapopen -E $gapext -B $subst  genome.fasta fakemates_1.fastq fakemates_2.fastq | samtools view -Sb - | samtools sort -o bwa_sorted.bam - &>> $outp
  fi

  if [ $rerun == "y" ] || [ $rerun == "yes" ] || [ ! -f "mates_id60_mscore10.out" ] ; then
    echo -e "\n
     Sort Aligned Fake-Mate Pairs by their Position in Reference\n"  >> $outp
    $bindir/smis_sort bwa_sorted.bam 
  fi
fi ##aligner selection


echo -e "\n
   Scaffolding Reference using Fake-Mate Pairs\n"  >> $outp

echo '################## Scaffolding ... ##################' > $outp
sed "s#MIN_EDGE#$min_edge#g" $bindir/../scafsettings.txt | sed  "s#NUM_SIGMA#$num_sigma#g" | sed  "s#MIN_LEN#$minlen#g" > settings.txt
printf "1\nTesting string\nfasta-in\ngenome.fasta\n" > files.txt

l1=(3000 5000 7000 9000 12000 20000)
l2=(500  1000 1500 2000 3000  4000)
ii=0
for gfile in genome-matepair-*; do
 echo $aligner >> files.txt
 
 nls=`wc -l $workdir/$gfile | awk '{print($1)}'`
 if [ $nls > 1000 ]; then
   echo $workdir/$gfile ${l1[$ii]} ${l2[$ii]} 151 70 in >> files.txt 
 fi
 ii=$(($ii+1))
done
printf "fasta-out\nspinner-sp2b.fasta\ncontigs-out\nspinner-contigs.dat\n">> files.txt 

$MYSMISDIR/smissv-bin/spinner_pg  > spinner.out

mv spinner-sp2b.fasta  ../spinner_scaffolds.fasta
cd ../

echo; echo " Scaffolds are in spinner_scaffolds.fasta"
echo " Summary of parameters used are in " $las
echo " Log is in " $outp



if [ ! $debug ]; then
 rm -rf $workdir/
fi


