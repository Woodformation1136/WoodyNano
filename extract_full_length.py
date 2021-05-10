import argparse
import sys
from datetime import datetime
from WoodyNano import main
from WoodyNano import seqtools

parser = argparse.ArgumentParser(description='Arguments availible to WoodyNano')
parser.add_argument('-l', metavar='len_cutoff', type=int,
                    help='Reads length cutoff. Default = sum(ap_length)')
parser.add_argument('-q', metavar='q_cutoff', default=7, type=int,
                    help='Reads qscore cutoff.')
parser.add_argument('-e', metavar='error_rate_cutoff', default=0.3, type=float,
                    help='Primer error rate cutoff.')
parser.add_argument('-p', metavar='primer_fasta', help='Primer fasta.')
parser.add_argument('--ap_length', nargs=2, metavar=("front", "end"), default=[130, 60], type=int,
                    help='AP length, separated by space')
parser.add_argument('--primer_cnfg', metavar='primer_cnfg', required=False, help='Primer configuration.')
parser.add_argument('input_fastq', metavar='input_fastq', help='Input fastq.')
parser.add_argument('output_fastq', metavar='output_fastq',
                    help='Output fastq.')
args = parser.parse_args()

def print_to_stdout(*a):
    print(*a, file=sys.stdout)


def print_to_stderr(*a):
    print(*a, file=sys.stderr)

def classifier(
        input_dir,
        output_dir,
        fp, rp,
        len_cutoff,
        q_cutoff,
        ap_length,
        error_cutoff,
        primer_cnfg
    ):

    # report start time and variables
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print_to_stdout(f'Work start...\nCurrent time: {current_time}')
    print_to_stdout('Parameters:')
    print_to_stdout(f'--input file: {input_dir}')
    print_to_stdout(f'--output file: {output_dir}')
    print_to_stdout(f'--forward primer (fp): {fp}')
    print_to_stdout(f'--reverse primer (rp): {rp}')
    print_to_stdout(f'--read length cutoff: {len_cutoff}')
    print_to_stdout(f'--qscore cutoff: {q_cutoff}')
    print_to_stdout(f'--adaptor region: 1:{ap_length[0]}, -{ap_length[1]}:end')
    print_to_stdout(f'--primer error rate cutoff: {error_cutoff}')
    print_to_stdout(f'--primer configuration: {primer_cnfg}\n')
    
    # read fastq
    raw_fastq = seqtools.SeqFastq.Import(fname=input_dir)
    
    # record read names
    ttl_reads = list(raw_fastq.keys())
    fusion_reads = list()
    print_to_stdout(f"{len(ttl_reads)} raw reads.")

    # qc
    for rname in ttl_reads:
        read = raw_fastq[rname]
        if len_cutoff > 0:
            if len(read.seq) < len_cutoff:
                raw_fastq.pop(rname)
            else:
                if q_cutoff > 0:
                    if read.mean_q() < q_cutoff:
                        raw_fastq.pop(rname)

    print_to_stdout(f"{raw_fastq.__len__()} reads passed QC")

    qc_passed_reads = list(raw_fastq.keys())
    
    for rname in qc_passed_reads:
        read = raw_fastq[rname]
        # align primer
        main.primer_alignment(
            read=read,
            fp=fp,
            rp=rp,
            score=error_cutoff,
            ap_length=ap_length
        )
        
        # check if body primers exists
        main.is_fusion(read)
        main.is_full_length(read, primer_cnfg)
        main.generate_cutpoint(read)
        main.cut_primer(read)
        main.read_orientation(read,direction='+')
        main.cut_ployA(read)
        
        if read.is_fusion:
            fusion_reads.append(rname)
        
        if not read.cut_ploy:
            raw_fastq.pop(rname)

    print_to_stdout(f"{raw_fastq.__len__()} full-length reads.")

    # export full-length fastq
    seqtools.SeqFastq.Export(
        fastq_dict=raw_fastq,
        fname=output_dir
        )
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print_to_stdout(f'Finished...\nCurrent time: {current_time}')

    return

if args.p:
    primer_fasta = main.read_primer_fasta(fasta_fname=args.p)
    primers = {
        'fp':primer_fasta['fp'].seq,
        'rp':primer_fasta['rp'].seq
    }
else:
    primers = {
        'fp':'AAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTC',
        'rp':'AAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGC'
    }

if not args.l:
    len_cutoff = sum(args.ap_length)
else:
    len_cutoff = args.l

if not args.primer_cnfg:
    primer_cnfg = {
        '+':['rp','fp_rc'],
        '-':['fp','rp_rc']
    }
else:
    primer_cnfg = args.primer_cnfg

classifier(
    input_dir=args.input_fastq,
    output_dir=args.output_fastq,
    fp=primers['fp'],
    rp=primers['rp'],
    len_cutoff=len_cutoff,
    q_cutoff=args.q,
    ap_length=args.ap_length,
    error_cutoff=args.e,
    primer_cnfg=primer_cnfg
)
