import re
import os
import numpy as np
import sys
from math import log

class SeqFasta:
    def __init__(
            self,
            info=str(),
            seq=str()
        ):
        self.info = info
        self.seq = seq
    
    def __repr__(self):
        return f"info: {self.info}\nseq:\n{self.seq}"
    
    @staticmethod
    def Import(fname):
        fasta = dict()
        with open(fname, 'r') as f:
            for line in f:
                if line[0] == '>':
                    # assign the old new_read to dictionary fastq
                    try:
                        fasta[new_read.info] = new_read
                    except NameError:
                        pass
                    
                    # creat a SeqFasta object
                    new_read = SeqFasta()
                    # record fasta info
                    new_read.info = line.strip()[1:]
                else:
                    new_read.seq = line.strip()
            fasta[new_read.info] = new_read

        return fasta
                

class SeqFastq:
    def __init__(self,
                 info=str(),
                 seq=str(),
                 info2=str(),
                 qscore=str()
                 ):
        self.info = info
        self.info2 = info2
        self.adprimer = {
            'fp': None,
            'rp': None,
            'fp_rc': None,
            'rp_rc': None
        }
        self.bdprimer = {
            'fp': None,
            'rp': None,
            'fp_rc': None,
            'rp_rc': None
        }
        self.body_hit = None
        self.seq = seq
        self.qscore = qscore
        self.strand = None
        self.cutpoint = None
        self.cut_primer = False
        self.cut_ploy = False
        self.orientation = None

    def __repr__(self):
        return f"info: {self.info}\
                \ncutpoint: {self.cutpoint}\
                \nadaptor: {self.adprimer}\
                \nreadbody: {self.bdprimer}\
                \nbody_hit: {self.body_hit}\
                \nstrand: {self.strand}\
                \ncut_primer: {self.cut_primer}\
                \ncut_ployA: {self.cut_ploy}\
                \norientation: {self.orientation}\
                \nseq:\n{self.seq}\
                \nqscore:\n{self.qscore}"

    def mean_q(self):
        p = [10**(-(ord(q)-33)/10) for q in self.qscore]
        mq = -10 * log(np.mean(p),10)
        return mq 

    @staticmethod
    def Import(fname):
        fastq = dict()
        
        with open(fname) as f:
            n = 0
            
            for l in f:
                n += 1
                new_read = SeqFastq()
                
                # check name line start with @    
                if not '@' == l[0]:
                    print(f"The {n}th line of '{fname}' does not start with '@'")
                    break
                else:
                    new_read.info = re.split('@', l.strip())[1]
                
                # read seq line
                new_read.seq = f.readline().strip()
                n += 1
                
                # read second info line
                new_read.info2 = f.readline().strip()
                n += 1
                
                # read qscore line
                new_read.qscore = f.readline().strip()
            
                fastq[str(new_read.info)] = new_read
        
        return fastq
    
    @staticmethod
    def Export(fastq_dict, fname):
        with open(fname, 'w') as o:
            for i in list(fastq_dict.keys()):
                read = fastq_dict[i]
                o.write(read.info+'\n')
                o.write(read.seq+'\n')
                o.write(read.info2+'\n')
                o.write(read.qscore+'\n')
