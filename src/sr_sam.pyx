#!/usr/bin/env python

from sys import argv, stderr, stdin, exit, stdout
from tempfile import mkdtemp
from shutil import rmtree
from string import maketrans
from subprocess import Popen, PIPE, call
from os import path

import gzip 

import pysam

class Alignment:
    '''Minimal representation of a single SAM alignment.
    '''
    def __init__(self, alignment):
        self.query_name = alignment.query_name
        assert alignment.query_sequence != None
        self.query_sequence = alignment.query_sequence
        self.query_qualities = alignment.query_qualities
        self.flag = alignment.flag
        self.reference_id = alignment.reference_id
        self.reference_start = alignment.reference_start
        self.mapping_quality = alignment.mapping_quality
        self.cigar = alignment.cigar
        self.next_reference_id = alignment.next_reference_id
        self.next_reference_start= alignment.next_reference_start
        self.template_length= alignment.template_length
        self.tags = alignment.tags
        self.query_length = alignment.infer_query_length()
        self.is_read1 = alignment.is_read1
        self.is_read2 = alignment.is_read2
        self.is_reverse = alignment.is_reverse
    
    def PrintAsSam(self, outfile):
        a = pysam.AlignedSegment()
        a.query_name = self.query_name
        a.query_sequence = self.query_sequence
        a.flag = self.flag
        a.reference_id = self.reference_id
        a.reference_start = self.reference_start
        a.mapping_quality = self.mapping_quality
        a.cigar = self.cigar
        a.next_reference_id = self.next_reference_id
        a.next_reference_start = self.next_reference_start
        a.template_length = self.template_length
        a.query_qualities = self.query_qualities
        a.tags = self.tags
        outfile.write(a)
    
def ReverseComplement(sequence):
    complement = maketrans('atcgnATCGN', 'tagcnTAGCN')
    return sequence.translate(complement)[::-1]

def PrintFastqSequence(alignment, fqfile):
    qname = alignment.query_name
    if alignment.is_read2:
        qname += "/2"
    else:
        assert alignment.is_read1 or (not alignment.is_paired)
        qname += "/1"

    sequence = alignment.query_sequence
    if alignment.is_reverse:
        sequence = ReverseComplement(sequence)
    
    quality = alignment.query_qualities
    if alignment.is_reverse:
        quality = quality[::-1]
    quality = [chr(x+33) for x in quality]

    print >> fqfile, qname
    print >> fqfile, sequence
    print >> fqfile, "+"
    print >> fqfile, "".join(quality)

def AddSplitReadsToFq(split_name, working_directory):
    bamfile = pysam.AlignmentFile(split_name, "rb")
    fqfile  = open("%s/reads.fq" % working_directory, "w")

    numclipped = 0
    for alignment in bamfile:
        if alignment.is_qcfail: continue
        if alignment.is_secondary: continue
        if alignment.is_supplementary: continue
        if alignment.is_duplicate: continue
        if alignment.is_unmapped: continue

        PrintFastqSequence(alignment, fqfile)
        numclipped += 1
    bamfile.close()

    print >> stderr, "Written %d clipped sequences to input fq" % numclipped

def RunCommand(command):
    try:
        call(command)
    except OSError as e:
        print >> stderr, "Execution of %s failed: %s", (command,e)
        exit(6)

def AlignReadsUsingLASTZ(cwd, lastz, threads, ref_name):
    DIR = path.dirname(path.realpath(__file__))

    command = ["%s/align_reads_with_lz" % DIR,
               "-t", str(threads),
               lastz,
               ref_name,
               "%s/reads.fq" % cwd,
               cwd
              ]
    RunCommand(command)

def DecideOnBwaAlignments(alignments, samfile, fqfile, maxsplits):
    alength = len(alignments)
    is_candidate = False

    if alength == 1:
        PrintFastqSequence(alignments[0], fqfile)
        is_candidate = True
    elif alength >= 2 and alength <= maxsplits:
        for aln in alignments:
            aln.PrintAsSam(samfile)
    else:
        PrintFastqSequence(alignments[0], fqfile)
        is_candidate = True

    return is_candidate

def ProcessBwaAlignments(split_name, working_directory, outfile, maxsplits):
    # sort the split file by name
    namesorted = "%s/split.namesorted" % working_directory
    pysam.sort("-n", split_name, namesorted)
    bam_name = namesorted + ".bam"

    bamfile = pysam.AlignmentFile(bam_name, "rb")
    fqfile  = open("%s/reads.fq" % working_directory, "w")
    query_name = None
    alignments = []        
    numclipped = 0

    for alignment in bamfile:
        aln = Alignment(alignment)

        if query_name == None:
            query_name = aln.query_name
            alignments.append(aln)
        elif aln.query_name == query_name:
            alignments.append(aln)
        else:
            is_candidate = DecideOnBwaAlignments(alignments,outfile,fqfile,maxsplits)
            if is_candidate:
                numclipped += 1

            query_name = aln.query_name
            alignments = [aln]

    if query_name != None:
        is_candidate = DecideOnBwaAlignments(alignments,outfile,fqfile,maxsplits)
        if is_candidate:
            numclipped += 1

    bamfile.close()
    fqfile.close()
    print >> stderr, "Written %d clipped sequence to input fq" % numclipped

def PrintSamHeader(faname):
    print "@HD\tVN:1.3"
    reference = pysam.FastaFile(faname)
    for chrom,chromlength in zip(reference.references,reference.lengths):
        print "@SQ\tSN:%s\tLN:%d" % (chrom,chromlength)

def FindSplitAlignments(maxsplits, coverage, lastz, threads, onlylz, ref_name, unmap_name, split_name = None):
    # create a temporary working directory
    working_directory = mkdtemp(dir = '.')

    # the output SAM records
    if split_name != None:
        bamfile = pysam.AlignmentFile(split_name, "r")
        outfile = pysam.AlignmentFile("-", "wh", template = bamfile)
        bamfile.close()
   
        if onlylz:
            # Add the split reads to the input to LASTZ
            AddSplitReadsToFq(split_name, working_directory)     
        else:
            # Select the clipped reads that need to be realigned using LASTZ
            ProcessBwaAlignments(split_name,working_directory,outfile,maxsplits)
    else:
        PrintSamHeader(ref_name)
    stdout.flush()

    # add the unmapped sequences to the file
    numunmap = 0
    fqfile  = open("%s/reads.fq" % working_directory, "a")
    f = gzip.open(unmap_name, "r")
    for line in f:
        print >> fqfile, line,
        numunmap += 1
    f.close()
    print >> stderr, "Written %d unmapped sequences to input fq" % (numunmap / 4)
    fqfile.close() 

    # make sure that we have at least $threads sequences. If not we just use one
    # threads
    numlines = 0
    f  = open("%s/reads.fq" % working_directory, "r")
    for line in f:
        numlines += 1
    f.close()
    if (numlines / 4) < threads:   
        threads = 1        

    # now use LASTZ to align the sequences in the input fq 
    AlignReadsUsingLASTZ(working_directory, lastz, threads, ref_name)

    # clear the temporary working area
    rmtree(working_directory)
