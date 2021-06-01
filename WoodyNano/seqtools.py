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
                 info: str,
                 seq: str,
                 info2: str,
                 qscore: str
                 ) -> None:
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
    def read(fileObject):
        buffer = [fileObject.readline().strip() for _ in range(4)]

        if any([len(i) == 0 for i in buffer]):
            # print("EOF")
            return None

        else:
            if not buffer[0][0] == '@':
                raise Exception('Error: Missing lines')

            else:
                return SeqFastq(
                    info=buffer[0],
                    seq=buffer[1],
                    info2=buffer[2],
                    qscore=buffer[3]
                )

    def write(self, fileObject):
        fileObject.write(self.info+'\n')
        fileObject.write(self.seq+'\n')
        fileObject.write(self.info2+'\n')
        fileObject.write(self.qscore+'\n')
        return


    @staticmethod
    def Import(fname):
        fastq = dict()
        
        with open(fname) as f:
            new_read = SeqFastq.read(f)
            while new_read:
                fastq[new_read.info] = new_read
                new_read = SeqFastq.read(f)
        
        return fastq

    
    @staticmethod
    def Export(fastq_dict, fname):
        with open(fname, 'w') as o:
            for i in list(fastq_dict.keys()):
                fastq_dict[i].write(o)

    
    

