
Scaffolding pipeline using data from long reads technologies (PacBio, ONT)
to scaffold an initial draft assembly. The long reads are shred in smaller segments 
(f.i. 1000 bp) to create fake mate-pairs. The fake mates are
then aligned against the draft assembly and the spinner scaffolder looks for
links between contigs and creates scaffolds. 

### Download and Compile:
Requirements for compiling: Cmake > = 2.6.4

	$ git clone https://github.com/SangerHpag/smis.git
	$ cd smis 
	$ ./makeall.sh

(Tested with gcc-4.9.2, bamtools-2.4.0) 

#### External packages
The smis pipeline downloads and installs the bamtools for reading bam files (https://github.com/pezmaster31/bamtools) 


### Run 
#### Setup 

	$MYSMISDIR/setup.sh </full/path/to/destdir> <draft_assembly> <long_reads>

	where:
   	   /full/path/to/destdir: folder where to run the pipeline (Please provide full path)
   	   draft assembly: fasta file of the assembly to be scaffolded  (Please provide full path)
	   long reads: fastq file of long reads for scaffolding (Please provide full path)

#### Parameters
The pipeline parameters can be modified in the /full/path/to/destdir/mysettings.sh .
The default aligner is bwa. Change to smalt by changing the 'aligner' variable in settings.sh
   
#### Run:
Requirements for running: samtools, bwa (or smalt) in PATH.

	cd /full/path/to/destdir
   	./mysmissv.sh

(Tested with samtools-1.3.1, bwa-0.7.12, smalt-0.7.4)

#### Results

Scaffolds will be in /full/path/to/destdir/spinner_scaffolds.fasta


