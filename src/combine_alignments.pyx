#!/usr/bin/env python

"""
    usage:
        combine_alignments [options] name2col numthreads reads.fq

    where the options are:
        -h,--help : print usage and quit
        -d,--debug: print debug information
        -c,--comment: lines in alignment files beginning with this character
                      are treated as comments    

    Combine several alignment files produced by LASTZ as a result of running
    align_reads_with_lz script. This assumes (and does not check for) that 
    all alignment files have the same format (which is accurate if this is
    processing the output of align_reads_with_lz).
"""

from sys import argv, stderr, stdin, exit, stdout
from getopt import getopt, GetoptError

__author__ = "Aakrosh Ratan"
__email__  = "ratan@bx.psu.edu"

# do we want the debug information to be printed?
debug_flag = False

def NextFqSequence(reads_name):
    with open(reads_name, "r") as f:
        line = f.readline()
        assert line.startswith("@"), "%s"% line
        while line:
            name = line.strip()[1:]
            line = f.readline()
            line = f.readline()
            assert line.startswith("+")
            line = f.readline()     
            yield name
            line = f.readline()
            if not line: break
            assert line.startswith("@")

def Combine(name2col, numthreads, reads_name, tmpdir, comment):
    # keep the file handles for all the alignment sequences
    files=[open("%s/alignments.%d.lz" % (tmpdir,indx), "r") for indx in xrange(1,numthreads+1)]

    cache = []
    for indx,file in enumerate(files):
        line = file.readline()
        while line.startswith(comment):
            #if indx == 0: print line,
            line = file.readline()
        cache.append(line)   

    # iterate through the sequences in the reads file and output the alignments
    # in order
    indx = 1
    for name in NextFqSequence(reads_name):
        line = cache[indx-1]
        while line and line.strip().split("\t")[name2col-1] == name:
            print line,
            line = files[indx-1].readline()
        cache[indx-1] = line
        indx += 1
        if indx > numthreads:
            indx = 1
