from WoodyNano import seqtools
from WoodyNano import utils


def read_primer_fasta(fasta_fname):
    f = seqtools.SeqFasta.Import(fname=fasta_fname)
    if 'fp' in f.keys() and 'rp' in f.keys():
        return f
    else:
        print('Unrecognized primer name')
        return


def primer_alignment(read, fp, rp, score, ap_length):

    primers = {
        'fp':fp,
        'rp':rp,
        'fp_rc': utils.reverse_complement(fp),
        'rp_rc': utils.reverse_complement(rp)
    }

    for pos,pnames in zip(['front', 'end'],[('fp','rp'),('fp_rc','rp_rc')]):
        if pos == 'front':
            ap_region = [0, ap_length[0]]
        if pos == 'end':
            ap_region = [len(read.seq)-ap_length[1], len(read.seq)]

        for p in pnames:
            # map to the adaptor region
            read.adprimer[p] = utils.aligner(
                primer=primers[p],
                sequence=read.seq,
                score=score,
                map_start=ap_region[0],
                map_end=ap_region[1]
            )
            # map to read body
            read.bdprimer[p] = utils.aligner(
                primer=primers[p],
                sequence=read.seq,
                score=score,
                map_start=ap_length[0],
                map_end=(len(read.seq) - ap_length[1])
            )
    
    return


def body_hit(read):
    try:
        read.body_hit = any([read.bdprimer[p]['editDistance'] !=
                              99 for p in read.bdprimer.keys()])
    except TypeError:
        pass

    return


def split_fusion(read, primer_cnfg):
    def correct_primer(read, pname1, pname2):
        correct = None
        err1 = read.adprimer[pname1]['errorRate']
        err2 = read.adprimer[pname2]['errorRate']

        if err2 > err1:
            correct = pname1
        elif err1 > err2:
            correct = pname2

        return correct

    correct_front = correct_primer(read, 'fp', 'rp')
    correct_end = correct_primer(read, 'fp_rc', 'rp_rc')
    
    paired_primer = {
        'fp':'rp_rc',
        'rp':'fp_rc',
        'fp_rc':'rp',
        'rp_rc':'fp'
    }
    split_read1 = None
    split_read2 = None

    if correct_front:
        if read.bdprimer[paired_primer[correct_front]]['editDistance'] != 99:
            middle_primer1 = paired_primer[correct_front]
            split_read1 = seqtools.SeqFastq(
                info=read.info + '_1',
                seq=read.seq,
                info2=read.info2,
                qscore=read.qscore
            )
            
            if correct_front == primer_cnfg['+'][0]:
                split_read1.strand = '+'
            
            elif correct_front == primer_cnfg['-'][0]:
                split_read1.strand = '-'

            split_read1.adprimer[middle_primer1] = read.bdprimer[middle_primer1]
            split_read1.adprimer[correct_front] = read.adprimer[correct_front]

    if correct_end:
        if read.bdprimer[paired_primer[correct_end]]['editDistance'] != 99:
            middle_primer2 = paired_primer[correct_end]
            split_read2 = seqtools.SeqFastq(
                info=read.info + '_2',
                seq=read.seq,
                info2=read.info2,
                qscore=read.qscore
            )

            if correct_end == primer_cnfg['+'][1]:
                split_read2.strand = '+'
            
            elif correct_end == primer_cnfg['-'][1]:
                split_read2.strand = '-'

            split_read2.adprimer[middle_primer2] = read.bdprimer[middle_primer2]
            split_read2.adprimer[correct_end] = read.adprimer[correct_end]

    return [split_read1, split_read2]

def is_full_length(read, primer_cnfg):
    
    def correct_primer(read, pname1, pname2):
        correct = None
        err1 = read.adprimer[pname1]['errorRate']
        err2 = read.adprimer[pname2]['errorRate']
        
        if err2 > err1:
            correct = pname1
        elif err1 > err2:
            correct = pname2
        
        return correct    

    if not read.body_hit:
        correct_front = correct_primer(read, 'fp', 'rp')
        correct_end = correct_primer(read, 'fp_rc', 'rp_rc')
        
        if correct_front and correct_end:
            if [correct_front, correct_end] == primer_cnfg['-']:
                read.strand = '-'

            elif [correct_front, correct_end] == primer_cnfg['+']:
                read.strand = '+'

    return


def generate_cutpoint(read):
    # cut 'GGG'/'CCC' at rp and rp_rc
    if read.strand:
        if read.strand == '-':
            read.cutpoint = (
                read.adprimer['fp']['locations'][-1][-1],
                read.adprimer['rp_rc']['locations'][0][0]-4
            )
        
        if read.strand == '+':
            read.cutpoint = (
                read.adprimer['rp']['locations'][-1][-1]+5,
                read.adprimer['fp_rc']['locations'][0][0]
            )
    return


def cut_primer(read):
    if not read.cut_primer:
        if read.cutpoint:
            read.seq = read.seq[read.cutpoint[0]:read.cutpoint[1]]
            read.qscore = read.qscore[read.cutpoint[0]:read.cutpoint[1]]
            read.cut_primer = True
    return


def read_orientation(read, direction='+'):
    if read.strand != direction:
        read.seq = utils.reverse_complement(read.seq)
        read.qscore = read.qscore[::-1]
        read.orientation = direction

    if read.strand == direction:
        read.orientation = direction

    return

def cut_ployA(read):
    if read.cut_primer:
        if not read.cut_ploy:
            if read.orientation == '+':
                map_end = -7
                tail = read.seq[map_end:map_end+6]
                while tail.count('A') >= 4:
                    map_end -= 1
                    tail = read.seq[map_end:map_end+6]
                map_end += 3
                
                if read.strand == '+':
                    f = read.cutpoint[0]
                    e = read.cutpoint[1] + map_end
                    read.cutpoint = (f, e)
                
                if read.strand == '-':
                    f = read.cutpoint[0] - map_end
                    e = read.cutpoint[1]
                    read.cutpoint = (f, e)

                read.seq = read.seq[:map_end]
                read.qscore = read.qscore[:map_end]

                read.cut_ploy = True
                read.info = '@%s:%s|%s strand=%s' % (
                    read.cutpoint[0], read.cutpoint[1], read.info, read.strand)
    return

