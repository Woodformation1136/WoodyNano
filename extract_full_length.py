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
parser.add_argument('--min_length', metavar='', default=150, type=int,
        help='Minimal output read length. (Default: 1)')
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
    min_length,
    primer_cnfg,
    out_fusion
):

    # report start time and variables
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print_to_stdout(f'Work start...\nCurrent time: {current_time}\n')
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

    stats = {
        'InputReads': 0,
        'LengthFailed': 0,
        'QscoreFailed': 0,
        'PassedQC': 0,
        'FullLengthReads': 0,
        'BodyHits': 0,
        'FusionReads': 0
    }

    output = open(output_dir, 'w')
    output = open(output_dir, 'a')

    if out_fusion:
        fusion = open(out_fusion, 'w')
        fusion = open(out_fusion, 'a')

    with open(input_dir, 'r') as f:
        read = seqtools.SeqFastq.read(f)
        while read:
            stats['InputReads'] += 1
            short = False
            low_quality = False

            # check length
            if len(read.seq) < len_cutoff:
                stats['LengthFailed'] += 1
                short = True

            # check qscore
            if read.mean_q() < q_cutoff:
                stats['QscoreFailed'] += 1
                low_quality = True

            # if short or low_quality
            if short or low_quality:
                # get next read
                read = seqtools.SeqFastq.read(f)
                continue

            stats['PassedQC'] += 1

            # align primers, check bodyprimers
            main.primer_alignment(
                read=read,
                fp=fp,
                rp=rp,
                score=error_cutoff,
                ap_length=ap_length
            )
            main.body_hit(read)

            # body_hit == True
            if read.body_hit:
                stats['BodyHits'] += 1

                # if output fusion reads
                if out_fusion:
                    # split read into read1, read2
                    split_read1, split_read2 = main.split_fusion(
                        read=read,
                        primer_cnfg=primer_cnfg
                    )

                    # for new_read in read1, read2
                    for new_read in split_read1, split_read2:

                        # if new_read exist, new_read must be full-length
                        if new_read:
                            # trimm new_read, orient new_readmain.generate_cutpoint(split_read)
                            main.cut_primer(new_read)
                            main.read_orientation(new_read, direction='+')
                            main.cut_ployA(new_read)

                            # output, make sure trimming were done
                            if new_read.cut_ploy:
                                if len(new_read.seq) >= min_length:
                                    stats['FusionReads'] += 1
                                    new_read.write(fusion)

            # body_hit == False
            else:
                # check if read is full-length
                main.is_full_length(read, primer_cnfg)

                # if read is full-length
                if read.strand:
                    # trim read, orient read
                    main.generate_cutpoint(read)
                    main.cut_primer(read)
                    main.read_orientation(read, direction='+')
                    main.cut_ployA(read)

                    # output, make sure trimming were done
                    if read.cut_ploy:
                        if len(read.seq) >= min_length:
                            stats['FullLengthReads'] += 1
                            read.write(output)

            # get the next read
            read = seqtools.SeqFastq.read(f)

    output.close()
    if out_fusion:
        fusion.close()
    
    print_to_stdout(f"InputReads: {stats['InputReads']}")
    print_to_stdout(f"LengthFailed: {stats['LengthFailed']}")
    print_to_stdout(f"QscoreFailed: {stats['QscoreFailed']}")
    print_to_stdout(f"PassedQC: {stats['PassedQC']}")
    print_to_stdout(f"FullLengthReads: {stats['FullLengthReads']}")
    print_to_stdout(f"BodyHits: {stats['BodyHits']}")
    print_to_stdout(f"FusionReads: {stats['FusionReads']}")
    print_to_stdout("")


    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print_to_stdout(f'Finished...\nCurrent time: {current_time}')

    return

# get arguments
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
    min_length=args.min_length,
    primer_cnfg=primer_cnfg,
    out_fusion=out_fusion
)
