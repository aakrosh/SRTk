#!/usr/bin/env python

"""
    usage:
        find_best_split [options] reference.fa < alignments.sam

    where the options are:
        -h,--help         : print usage and quit
        -d,--debug        : print debug information
        -s,--significant  : only output alignments where the best possible split
                            is better than the second best by at least this 
                            score [10]
        -c,--coverage     : require at least this % of the sequence to be 
                            aligned [50]                       
        -m,--maxalignments: ignore the sequence if more than this number of 
                            qualified alignments are found [2000]
        -i,--ignoreq      : ignore the mapping qualities from the aligner
"""

from sys import argv, stderr, stdin, exit, stdout
from getopt import getopt, GetoptError
from time import clock
from string import maketrans

import numpy as np

import pysam 

__author__ = "Aakrosh Ratan"
__email__  = "ratan@bx.psu.edu"

# do we want the debug information to be printed?
debug_flag = False

NEGATIVEINFINITY = np.iinfo('i').min

MATCH = 1
MISMATCH = -1
GAPOPEN = -4
GAPEXTEND = -3
JUMP = -15

MQTHRESHOLD = 1
NEIGHBORHOOD = 5 

class Alignment:
    '''Minimal representation of a single SAM alignment.
    '''
    def __init__(self, alignment, reference, infile):
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

        # these are some of the attributes that I need, and their definition
        # might not be consistent with the definition in pysam
        if alignment.is_reverse:
            self.query_start = self.query_length - alignment.query_alignment_end
            self.query_end = self.query_length - alignment.query_alignment_start
        else:
            self.query_start = alignment.query_alignment_start
            self.query_end = alignment.query_alignment_end   
        self.reference_name = infile.getrname(self.reference_id)
        self.reference_end = alignment.reference_end
        self.is_reverse = alignment.is_reverse
        
        # some other details from the alignments, such as whether a particular
        # base from the query was part of a mismatch, match or a gap.
        self.matches    = np.zeros(self.query_length, dtype=np.int8)
        self.mismatches = np.zeros(self.query_length, dtype=np.int8)
        self.rgaps      = np.zeros(self.query_length, dtype=np.int8)
        self.qgaps      = np.zeros(self.query_length, dtype=np.int8)
        self.pairs = []

        refseq = reference.fetch(self.reference_name, 
                                 alignment.reference_start,
                                 alignment.reference_end)

        qi = self.query_start - 1
        fseg = True
        for q,r in alignment.get_aligned_pairs():
            if fseg == True:
                if r == None:
                    continue
                else:
                    fseg = False
            if qi == (self.query_end - 1): break

            self.pairs.append((q,r))
            qbase = "-" if q == None else self.query_sequence[q].upper()
            rbase = "-" if r == None else refseq[r - self.reference_start].upper()
            
            if   rbase == "-" and qbase != "-":
                qi += 1
                self.rgaps[qi] = 1
            elif rbase != "-" and qbase == "-":
                self.qgaps[qi] += 1
            elif rbase != "-" and qbase != "-":
                qi += 1
                if rbase == 'N' or qbase == 'N' or rbase == qbase:
                    self.matches[qi] = 1
                else:
                    self.mismatches[qi] = 1
            
    
        #offset = self.reference_start 
        #qlast  = None
        #for q,r in alignment.get_aligned_pairs():
        #    if q == None:
        #        x = qlast
        #    else:
        #        qlast = q
        #        x = q

        #    if alignment.is_reverse:
        #        pq = self.query_length - x - 1
        #    else:
        #        pq = x

        #    if pq != None and (pq < self.query_start or pq >= self.query_end): 
        #        continue

        #    assert pq >= self.query_start and \
        #           pq < self.query_end, "%d : %d : %d" \
        #        % (pq, self.query_start, self.query_end)
        #    self.pairs.append((q,r))
        #    rbase = "-" if r == None else refseq[r - offset].upper()
        #    qbase = "-" if q == None else self.query_sequence[q].upper()

        #    if alignment.is_reverse:q = pq
    
        #    if rbase == "-" and qbase != "-":
        #        self.rgaps[q] = 1
        #    elif rbase != "-" and qbase == "-":
        #        self.qgaps[qlast+1] += 1
        #    elif rbase != "-" and qbase != "-":
        #        if rbase == 'N' or qbase == 'N' or rbase == qbase:
        #            self.matches[q] = 1
        #        else:
        #            self.mismatches[q] = 1


    def __str__(self):
        string  = "%s\n" % self.query_name
        string += "\t%s:%d-%d\t%d-%d\n" % (self.reference_name, self.reference_start,self.reference_end,self.query_start,self.query_end)
        string += "\tmatches:\t%s\n" % "".join([str(x) for x in self.matches])
        string += "\tmismatches:\t%s\n" % "".join([str(x) for x in self.mismatches])
        string += "\trgaps:\t\t%s\n" % "".join([str(x) for x in self.rgaps])
        string += "\tqgaps:\t\t%s\n" % "".join([str(x) for x in self.qgaps])
        return string

def PrintSamOutput(alignments, outfile):
    '''Print the alignments as SAM records.
    '''
    for aln in alignments:
        a = pysam.AlignedSegment()
        a.query_name = aln.query_name
        a.query_sequence = aln.query_sequence
        a.flag = aln.flag
        a.reference_id = aln.reference_id
        a.reference_start = aln.reference_start
        a.mapping_quality = aln.mapping_quality
        a.cigar = aln.cigar
        a.next_reference_id = aln.next_reference_id
        a.next_reference_start = aln.next_reference_start
        a.template_length = aln.template_length
        a.query_qualities = aln.query_qualities
        a.tags = aln.tags
        outfile.write(a)

def CreateCigar(aln, start, end):
    r = aln.reference_start
    q = 0

    mstart = aln.query_length - end if aln.is_reverse else start
    mend   = aln.query_length - start if aln.is_reverse else end

    cigar = []
    if mstart > 0: cigar.append([4,mstart])

    for q,r in aln.pairs:
        if q != None and q < mstart: continue
        if q != None and q >= mend: continue

        if   r == None and q != None:
            if cigar and cigar[-1][0] == 1:
                cigar[-1][1] += 1
            else:
                cigar.append([1,1])
        elif r != None and q == None:
            if cigar and cigar[-1][0] == 2:
                cigar[-1][1] += 1
            else:
                cigar.append([2,1])
        elif r != None and q != None:
            if cigar and cigar[-1][0] == 0:
                cigar[-1][1] += 1
            else:
                cigar.append([0,1])
        else:
            print >> stderr, \
            "Should not happen. All conditions are covered for pairs."  

    if mend < aln.query_length:
        cigar.append((4,aln.query_length - mend))

    return cigar

def FindRefStart(aln, start, end):
    search = aln.query_length - end if aln.is_reverse else start

    for q,r in aln.pairs:
        if q == search: 
            return r
        
    return None

def CreateAlignment(alignments, row_indx, col_start, col_end):
    aln = alignments[row_indx]

    a = pysam.AlignedSegment()
    a.query_name = aln.query_name
    a.query_sequence = aln.query_sequence
    a.flag = aln.flag
    a.reference_id = aln.reference_id
    a.mapping_quality = aln.mapping_quality
    a.next_reference_id = aln.next_reference_id
    a.next_reference_start = aln.next_reference_start
    a.template_length = aln.template_length
    a.query_qualities = aln.query_qualities

    tmp = FindRefStart(aln, col_start, col_end)
    if tmp == None:
        return None   
    else:
        a.reference_start = tmp

    a.cigar = CreateCigar(aln, col_start, col_end)
    assert a.cigar != None
    return a

class Node:
    def __init__(self, row_indx = -1, col_indx = -1):
        self.row = row_indx
        self.col = col_indx
        self.score = NEGATIVEINFINITY
        self.best_jump = None
        self.is_used = False

def PrintScoringMatrix(network, m, n):
    print >> stderr, "Scoring matrix:"
    print >> stderr, "\t".join([str(x) for x in xrange(0,n)])
    print >> stderr, "-"* n
    for row_indx in xrange(0,m):
        row = []
        for col_indx in xrange(0,n):
            val = network[row_indx][col_indx].score
            if val == NEGATIVEINFINITY:
                row.append("NI")   
            else:
                row.append(str(val))
        print >> stderr, "\t".join(row)  
    print >> stderr, ""

def BetterToSkip(col_best_score_nodes, col_indx, nd):
    nodes = col_best_score_nodes[col_indx - NEIGHBORHOOD:col_indx]
    for node in nodes:
        if node.score > nd.score:
            return True

def RunSanityChecks(alignments):
    '''Run some simple checks to ensure that this works in realistic ways.
        
    The following checks are currently employed:
    1. If there are more than 2 splits, then all local alignments that are 
       less than 30M are removed.
    '''
    best_alignments = []
    if len(alignments) <= 2:
        return alignments
    else:   
        for aln in alignments:
            alnlength = AlnCoverage(aln)
            if alnlength > 30:
                best_alignments.append(aln)

    return best_alignments

def FindBestSplit(infile, alignments, significant):
    '''Return the best alignments that cover the query.
    '''
    global NEGATIVEINFINITY

    m = len(alignments)
    n = alignments[0].query_length 

    # the score matrix
    network = [[Node(row_indx,col_indx) for col_indx in xrange(0,n)] 
                                        for row_indx in xrange(0,m)]

    col_best_score_nodes = []
    col_second_best_score_nodes = []

    for col_indx in xrange(0,n):
        col_best_score_node = Node()
        col_second_best_score_node = Node()

        for row_indx in xrange(0,m):
            aln  = alignments[row_indx]
            node = network[row_indx][col_indx]
 
            # this alignment is only useful if we are within the aligned segment
            if col_indx < aln.query_start or col_indx >= aln.query_end:
                node.score = NEGATIVEINFINITY           
            else:
                # what if I want to continue on the same alignment?
                path1 = 0
                if col_indx != aln.query_start:
                    path1 = network[row_indx][col_indx-1].score
                    if aln.qgaps[col_indx] != 0:
                        path1 += (GAPOPEN + aln.qgaps[col_indx]*GAPEXTEND)
               
                if aln.matches[col_indx] == 1:
                    path1 += MATCH
                elif aln.mismatches[col_indx] == 1:
                    #assert col_indx != aln.query_start
                    path1 += MISMATCH
                elif aln.rgaps[col_indx] == 1:
                    assert col_indx != aln.query_start
                    if aln.rgaps[col_indx-1] == 1:
                        path1 += GAPEXTEND
                    else:
                        path1 += (GAPOPEN + GAPEXTEND)
                else:
                    print >> stderr, "path1: unhandled situation"

                # what if I jump from another alignment to this alignment?
                if col_indx == 0:
                    path2 =  NEGATIVEINFINITY
                else:
                    nd = col_best_score_nodes[-1]
                    if nd == None:
                        path2 =  NEGATIVEINFINITY
                    else:
                        path2 = col_best_score_nodes[-1].score
                        path2 += JUMP
                # any gap that I see here should be treated as a new gap start 
                if aln.matches[col_indx] == 1:
                    path2 += MATCH
                elif aln.mismatches[col_indx] == 1:
                    #assert col_indx != aln.query_start
                    path2 += MISMATCH
                elif aln.rgaps[col_indx] == 1:
                    assert col_indx != aln.query_start
                    path2 += (GAPOPEN + GAPEXTEND)
                else:
                    print >> stderr, "path2: unhandled situation"
            
                # if the jump is profitable and this is the beginning of the
                # alignment, then I should jump
                if path2 > path1 and col_indx == aln.query_start:
                    node.best_jump = nd.row
                    if debug_flag:
                        print >> stderr, "Jumping from %d to %d at column %d" \
                               % (node.row, nd.row, col_indx)
                node.score = max(path1,path2)

            if node.score >= col_best_score_node.score: 
                col_second_best_score_node = col_best_score_node
                col_best_score_node = node
            elif node.score >= col_second_best_score_node.score:
                col_second_best_score_node = node

        col_best_score_nodes.append(col_best_score_node)
        col_second_best_score_nodes.append(col_second_best_score_node)

    # lets print the score matrix
    if debug_flag: PrintScoringMatrix(network, m, n)

    # now lets find the set of alignments that best covers the query
    best_alignments = []
    col_indx = n-1
    last_col_score = NEGATIVEINFINITY
    is_good = True

    while col_indx >= 0:
        # which alignment is the best at this query base
        nd = col_best_score_nodes[col_indx]
        if nd.is_used:
            # we have used this alignment once, the best and second best
            # alignment are too close for comfort
            is_good = False
            break 
        if nd.score == NEGATIVEINFINITY:    
            # none of the alignments work at this query base
            col_indx -= 1
            continue
        elif col_indx > NEIGHBORHOOD:
            if BetterToSkip(col_best_score_nodes,col_indx,nd):
                col_indx -= 1
                continue

        # this could be an interesting segment of the alignment
        row_indx = nd.row
        end = col_indx + 1
    
        # lets ensure that the best score here is significantly better than the
        # second best score
        if (nd.score - col_second_best_score_nodes[col_indx].score)< significant:
            is_good = False   
            break    

        nd.is_used = True
        while nd.score != NEGATIVEINFINITY and col_indx >= 0:
            if nd.best_jump == None:
                col_indx -= 1
                nd = network[row_indx][col_indx]
            else:
                start = col_indx
                if debug_flag:
                    print >> stderr, "Segment: %d - %d" % (col_indx+1,end)
                aln = CreateAlignment(alignments, row_indx, start, end)
                if aln == None:
                    is_good = False
                    break
                best_alignments.append(aln)
                col_indx -= 1
                row_indx = nd.best_jump
                end = col_indx + 1
                nd = network[row_indx][col_indx]
                if nd.is_used:
                    is_good = False
                    break
                if (col_best_score_nodes[col_indx].score - \
                    col_second_best_score_nodes[col_indx].score) < significant:
                    is_good = False
                    break

                nd.is_used = True
        
        if nd.score == NEGATIVEINFINITY:
            if debug_flag: print >> stderr, "Segment: %d - %d" % (col_indx+1,end)
            aln = CreateAlignment(alignments, row_indx, col_indx+1, end)
            if aln == None:
                is_good = False
                break
            best_alignments.append(aln)
        if col_indx == 0:
            if nd.score != NEGATIVEINFINITY:
                if debug_flag: print >> stderr, "Segment: %d - %d" % (0,end)
                aln = CreateAlignment(alignments, row_indx, 0, end)
                if aln == None:
                    is_good = False
                    break
                best_alignments.append(aln)
            break   
               
    if not is_good:
        best_alignments = []
    best_alignments = RunSanityChecks(best_alignments)

    return best_alignments

def AlnCoverage(aln):
    qaln_length = 0
    
    for op,oplen in aln.cigartuples:
        if op in [0,1,7,8]:
            qaln_length += oplen
    
    return qaln_length


def Coverage(alignments):
    qaln_length  = 0

    for aln in alignments:
        for op,oplen in aln.cigartuples:
            if op in [0,1,7,8]:
                qaln_length += oplen

    return qaln_length * 100.0 / alignments[0].query_length

def BothEndsMap(alignments):
    '''
    Most of the false positives occur due the problems with the aligner.
    Specifically if I see cases where xSyMzS or in other words where some
    subsequence in the middle is mapped, but the edges aren't that seems to be
    a definite error. Lets ignore those for now
    '''
    lmapped = False
    rmapped = False
    for aln in alignments:
        operations = list(aln.cigartuples)
        if aln.is_reverse: operations = operations[::-1]
        if operations[0][0] in [0,7]:
            lmapped = True
        if operations[-1][0] in [0,7]:
            rmapped = True

    return lmapped & rmapped 
    
def ProcessAlignments(alignments, infile, outfile, significant, coverage, maxalignments, maxsplits):
    '''Process alignments for this query.
    '''
    start_time = clock()
    if debug_flag:
        print >> stderr, \
        "Processing %s. m = %d, n = %d" % (alignments[0].query_name, len(alignments), alignments[0].query_length)

    if len(alignments) == 1:
        # if this is the only alignment for this query, then we should just 
        # print it out in the SAM format, i.e if the alignment is of good
        # quality
        # assert alignments[0].flag & 256 != 256
        if alignments[0].mapping_quality >= MQTHRESHOLD:
            PrintSamOutput(alignments, outfile)
    elif len(alignments) < maxalignments:
        # are all the alignments of good quality?
        do_process = True
        for aln in alignments:
            if aln.mapping_quality < MQTHRESHOLD:
                do_process = False
                break
        if do_process == True:
            # we have at least two alignments. So lets pick ones that best covers
            # the query as per the algorithm.
            best_alignments = FindBestSplit(infile, alignments, significant)
            if best_alignments != None and len(best_alignments) > 0:
                if maxsplits >= 1 and len(best_alignments) == 1:
                    PrintSamOutput(best_alignments[::-1], outfile)
                elif len(best_alignments) >= 2 and \
                     Coverage(best_alignments) > coverage and \
                     BothEndsMap(best_alignments) and \
                     len(best_alignments) <= maxsplits:
                    PrintSamOutput(best_alignments[::-1], outfile)

    end_time = clock()
    if debug_flag:
        print >> stderr, "Processing %s took %0.3f sec" \
           % (alignments[0].query_name, end_time - start_time)

def FindBestSplitFromAlignments(faname, significant, coverage, maxalignments, maxsplits, d_flag):
    debug_flag = d_flag

    query_name = None
    alignments = []

    # the reference sequence. I need it so that I can differentiate between 
    # matches and mismatches
    reference = pysam.FastaFile(faname)

    # print out the SAM header
    infile  = pysam.AlignmentFile("-", "r")
    outfile = pysam.AlignmentFile("-", "wh", template = infile)

    for alignment in infile:
        if alignment.is_unmapped: continue
        if alignment.is_duplicate: continue
        if alignment.is_qcfail: continue

        aln = Alignment(alignment, reference, infile)
        if debug_flag:
            print >> stderr, aln
    
        if query_name == None:
            query_name = aln.query_name
            alignments.append(aln)
        elif aln.query_name == query_name:
            alignments.append(aln)
        else:
            ProcessAlignments(alignments, infile, outfile, significant, coverage, maxalignments, maxsplits)
            query_name = aln.query_name
            alignments = [aln]

    if query_name != None:
        ProcessAlignments(alignments, infile, outfile, significant, coverage, maxalignments, maxsplits)

    reference.close()
    stdout.flush() # to get rid of warn "close failed in file object destructor:"
