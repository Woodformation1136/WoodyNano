import argparse
import sys
from datetime import datetime
from WoodyNano import main
from WoodyNano import seqtools

parser = argparse.ArgumentParser(description='Arguments availible to WoodyNano')
parser.add_argument('-l', metavar='len_cutoff', type=int,
        help='Reads length cutoff. (Default: sum(ap_length))')
parser.add_argument('-q', metavar='q_cutoff', default=7, type=int,
        help='Reads qscore cutoff. (Default: 7)')
parser.add_argument('-e', metavar='error_rate_cutoff', default=0.3, type=float,
        help='Primer error rate cutoff. (Default: 0.3)')
parser.add_argument('-p', metavar='primer_fasta', help='Primer fasta. (Not required)')
parser.add_argument('--ap_length', nargs=2, metavar='', default=[130, 60], type=int,
        help='AP length, separated by space (Default: 130 60)')
parser.add_argument('--primer_cnfg', metavar='', required=False, 
                    help='Primer configuration. (Not required)')
parser.add_argument('--fusion_read', metavar='', required=False,
                    help='Write splitted fusion reads to this file.')
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
        primer_cnfg,
        out_fusion
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
    print_to_stdout(f'--primer configuration: {primer_cnfg}')
    print_to_stdout(f'--fusion read fastq: {out_fusion}\n')


    if out_fusion:
        splitted_reads = dict()
    
    # read fastq
    raw_fastq = seqtools.SeqFastq.Import(fname=input_dir)
    
    # record read names
    ttl_reads = list(raw_fastq.keys())

    print_to_stdout(f"{len(ttl_reads)} raw reads.")

    # qc
    for rname in ttl_reads:
        read = raw_fastq[rname]
        if len_cutoff > 0:
            if len(read.seq) <= len_cutoff:
                raw_fastq.pop(rname)
            else:
                if q_cutoff > 0:
                    if read.mean_q() < q_cutoff:
                        raw_fastq.pop(rname)

    print_to_stdout(f"{raw_fastq.__len__()} reads passed QC")

    qc_passed_reads = list(raw_fastq.keys())
    
    
    body_primer = 0
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
        main.body_hit(read)
        if read.body_hit:
            body_primer += 1
            # if output were required
            if out_fusion:
                split_read1, split_read2 = main.split_fusion(read, primer_cnfg)
                for split_read in split_read1, split_read2:
                    if split_read:
                        main.generate_cutpoint(split_read)
                        main.cut_primer(split_read)
                        main.read_orientation(split_read, direction='+')
                        main.cut_ployA(split_read)
                        
                        if split_read.cut_ploy:
                            splitted_reads[split_read.info] = split_read

            raw_fastq.pop(rname)
        
        else:
            main.is_full_length(read, primer_cnfg)
            main.generate_cutpoint(read)
            main.cut_primer(read)
            main.read_orientation(read,direction='+')
            main.cut_ployA(read)
            
            if not read.cut_ploy:
                raw_fastq.pop(rname)
        
        # end of for-loop

    if out_fusion:
        print_to_stdout(f"{splitted_reads.__len__()} reads were generated from splitting fusion reads.")
    
    print_to_stdout(f"{body_primer} reads have body pimers.")
    print_to_stdout(f"{raw_fastq.__len__()} full-length reads.")

    # export full-length fastq
    seqtools.SeqFastq.Export(
        fastq_dict=raw_fastq,
        fname=output_dir
        )
    
    # export splitted fusion reads fastq
    if out_fusion:
        seqtools.SeqFastq.Export(
            fastq_dict=splitted_reads,
            fname=out_fusion
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

out_fusion = args.fusion_read if args.fusion_read else None

classifier(
    input_dir=args.input_fastq,
    output_dir=args.output_fastq,
    fp=primers['fp'],
    rp=primers['rp'],
    len_cutoff=len_cutoff,
    q_cutoff=args.q,
    ap_length=args.ap_length,
    error_cutoff=args.e,
    primer_cnfg=primer_cnfg,
    out_fusion=out_fusion
)
