# SRTk

## SUMMARY
SRTk is a collection of scripts to find the split-read alignments of sequences.

## REQUIREMENTS
STRBait should work on any standard 64 bit Linux environment with 

- GCC
- Python version >= 2.7
- Cython

The following two python libraries should be installed on the system

- pysam 
- numpy

## ACKNOWLEDGEMENTS
SRTk uses LASTZ (http://www.bx.psu.edu/~rsharris/lastz/newer) >= 1.03.66 to
align the sequences to the reference genome. 

## INSTALLATION
python, gcc are normally pre-installed on most Linux systems. If not,
please go ahead and install them. Cython is required as well, and might need to
be installed on some systems.

Please do the following in the top level directory of the distribution:
```
make
```
This should create a `bin` folder and copy all the required binaries to it.
Please make sure that following files are in the `bin` folder:
```
align.sh
align_reads_with_lz
combine_alignments
combine_alignments.so
find_best_split
find_best_split.so
sr_sam
sr_sam.so
```

## DESCRIPTION
SRTk is a collection of scripts to output a SAM file with the split-read
alignments of sequences in a sample. The script of interest is called sr_sam

```
Find the best split alignments for reads that carry signatures of SVs.
    usage:
      sr_sam [options] reference.fa [splitters.bam] unmapped.fq.gz

    where the options are:
        -h,--help      : print usage and quit
        -d,--debug     : print debug information
        -m,--maxsplits : allow up to maxsplits primary and supplemental
                         alignments [2]
        -c,--coverage  : require at least this fraction of the read to be aligned
                         in the primary alignment [0.5]
        -l,--lastz     : the path to the LASTZ (32 bit version) binary [lastz_32]
        -t,--threads   : number of threads to use [1]
        -1,--onlylz    : only use LASTZ. default is to use both BWA and LASTZ

    1. reference.fa is the reference sequence in fasta format.
    2. splitters.bam is a indexed BAM file with the clipped reads
    3. unmapped.fq.gz is a zipped file of sequences that did not align
       to the reference
```
        
This script collects all the split read alignments from BWA or some other
aligner. It also uses LASTZ to align the unmapped reads to the reference, and
outputs the best set of split alignments that cover the query. The output can
then be given as input to LUMPY to find the SVs.

## TEST-DATASET
A test dataset is provided with the distribution in the `tests` folder.

Run
```
make
```

and that should run some simple tests to make sure that the program is working
as expected.
