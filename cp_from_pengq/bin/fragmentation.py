#!/bin/env python
"""Calculate reads distribution along trancripts by mapping result on genome.

Input: gpe gene structure; bam file (MUST BE INDEXED);
normalized length for output, default 25

Output: a 2 column matrix, with 1st column as factor, 2nd column as frequence,
thus can be used to draw bargraph.CI on R using package sciplot

Method:
    for each gpe, fetch unique reads with same location,
      for each read, check if read in the gpe, if yes, note the position
    calculate distribution on normalized length

    calculate mean distribution and output

Author: yansy
Revised: 130624
"""

import sys
import pysam
import argparse

class GenePredExt(object):
    ' a line of gpe file'
    def __init__(self, line, bincolumn=True):
        'initialize each column'
        self.record = line.strip().split()
        if bincolumn == True:
            start = 1
        else:
            start = 0
        self.name = self.record[start]
        self.chr = self.record[start + 1]
        self.strand = self.record[start + 2]
        self.txStart = int(self.record[start + 3])
        self.txEnd = int(self.record[start + 4])
        self.cdsStart = int(self.record[start + 5])
        self.cdsEnd = int(self.record[start + 6])
        self.exonCount = int(self.record[start + 7])
        self.exonStarts = [int(i) for i in \
                           self.record[start + 8].strip(',').split(',')]
        self.exonEnds = [int(i) for i in \
                         self.record[start + 9].strip(',').split(',')]
        self.name2 = self.record[start + 11]

    def exon_pos(self):
        'returns a list of exon position tuples like (a, b),0-based'
        pos_list = []
        for i in range(self.exonCount):
            pos_list.append((self.exonStarts[i],self.exonEnds[i]))
        return pos_list

    def transcript_length(self):
        'transcript length'
        length = 0
        for (st, ed) in self.exon_pos():
            length += (ed - st)
        return length
        
    def contain(self, pos):
        """if pos is in trancript region, return True,
           ignoring chromosome and strand"""
        for (st, ed) in self.exon_pos():
            if (pos >= st and pos < ed):
                return True
        return False
        
    def dist2left(self, pos):
        """distence from pos to left (not including pos),
pos must be contained in transcript"""
        dist = 0
        for (st, ed) in self.exon_pos():
            if pos > ed:
                dist += (ed - st)
            elif (pos >= st and pos < ed):
                dist += (pos - st)
            else:
                continue
        return dist


def norm_freq(numlist):
    'normalize frequence in list so that total freq will be 1'
    s = sum(numlist)
    for i,num in enumerate(numlist):
        numlist[i] = float(num)/s
    return numlist
    
def mean_freq(mtx):
    """for mtx having n lists, each m elements, return 1 list with m element,
each a mean for n identical position"""
    m = len(mtx)
    n = len(mtx[0])
    mean = []
    for j in range(n):
        s = 0
        for i in range(m):
            s += mtx[i][j]
        mean.append(s/m)
    return mean

def sample_evenly(lst, size=10):
    "take even sample from lst"
    lgth = len(lst)
    interval = lgth/size
    st = interval / 2
    samp = []
    for i in range(size):
        samp.append(lst[st + i*interval])
    return samp
    
def read_pos(read, gpe, n):
    """if sam read in gpe structure, return its start position and end position
on the transcript normalized by n; if not, return None"""
    THRESHOLD = 0.90
    try:
        if read.opt("XS") != gpe.strand:
            return None
    except KeyError:   # some reads don't have XS, treat as unstranded case.
        pass

    match_num = 0
    sample_pos = sample_evenly(read.positions, size=10)
    for rpos in sample_pos:
        if gpe.contain(rpos):
            match_num += 1
    if float(match_num)/len(sample_pos) < THRESHOLD:
        return None
    else:    ## add 1 to denominator to restrain index in range
        left_pos = int(gpe.dist2left(read.pos) /
                       float(gpe.transcript_length()+ 1)*n)
        right_pos = int(gpe.dist2left(read.pos + read.rlen -1) /
                        float(gpe.transcript_length() + 1)*n)
        if gpe.strand == "+":
            return (left_pos, right_pos)
        else:
            return (n - right_pos - 1, n - left_pos - 1)

        
def main(fgpe, fsam, has_bin, n=25, lenlimit=0):
    """print the 2 column n row mean frequece table to stdout
       fileter out gene less than 200 bp length or no reads mapped to"""
    allfreq = []
    for gpeline in fgpe:
        gpe = GenePredExt(gpeline, bincolumn=has_bin)
        
        if lenlimit and gpe.transcript_length() <= lenlimit:
            continue

        try:
            reads = fsam.fetch(reference=gpe.chr, start=gpe.txStart,
                               end=gpe.txEnd)
        except ValueError:
            continue

        pos = list("0"*n)
        for i in range(len(pos)):
            pos[i] = int(pos[i])

        for read in reads:
            nh_tophat = None
            try:
                nh_tophat = read.opt("NH")
            except KeyError:
                pass
            
            if nh_tophat and nh_tophat == 1:
                rpos = read_pos(read, gpe, n)
                if rpos:
                    st, ed = rpos[0], rpos[1]
                    while st <= ed:
                        pos[st] += 1
                        st += 1
            else:
                continue
            

        if sum(pos) == 0:
            continue

        normfreq = norm_freq(pos)
        allfreq.append(normfreq)
    
    try:
        assert len(allfreq) > 0
    except AssertionError:
        print "Is your bam indexed?"
        raise
        
    meanfreq = mean_freq(allfreq)

    fgpe.close()
    fsam.close()
    
    for i,freq in enumerate(meanfreq):
        print "%d\t%.8f" % (i+1, freq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--gpe",
                        help="gpe gene structure used to assess fragmentation")
    parser.add_argument("-b", action='store_true', 
                        help="specified if gpe has bin column")
    parser.add_argument("--bam", help="input bam file (must be indexed)")
    parser.add_argument("-n", type=int, default=25,
                        help="fraction length to be normalized to, default 25")
    parser.add_argument("--lenlimit", type=int, default=0,
                        help=("only take transcripts longer than this "
                              "into consideration. default inactivated."))
    args = parser.parse_args()

    fgpe = open(args.gpe, "r")
    fsam = pysam.Samfile(args.bam, "rb")
    main(fgpe, fsam, n=args.n, has_bin=args.b,
         lenlimit=args.lenlimit)
