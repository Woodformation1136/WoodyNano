# WoodyNano

WoodyNano is a tool for extracting clean reads (trimmed, oriented, full-length reads) from Nanopore SQK-PCS109 sequencing results.

## Requirements

* Edlib 

## Default usage

```
python extract_full_length.py [input_fastq] [output_fastq]
```


## Split fusion reads

```
python extract_full_length.py --fusion_read split_reads.fastq [input_fastq] [output_fastq]
```


## Full usage

```
usage: extract_full_length.py [-h] [-l len_cutoff] [-q q_cutoff]
                              [-e error_rate_cutoff] [-p primer_fasta]
                              [--ap_length ] [--primer_cnfg] [--fusion_read]
                              input_fastq output_fastq

Arguments availible to WoodyNano

positional arguments:
  input_fastq           Input fastq.
  output_fastq          Output fastq.

optional arguments:
  -h, --help            show this help message and exit
  -l len_cutoff         Reads length cutoff. (Default: sum(ap_length))
  -q q_cutoff           Reads qscore cutoff. (Default: 7)
  -e error_rate_cutoff  Primer error rate cutoff. (Default: 0.3)
  -p primer_fasta       Primer fasta. (Not required)
  --ap_length           AP length, separated by space (Default: 130 60)
  --primer_cnfg         Primer configuration. (Not required)
  --fusion_read         Write splitted fusion reads to this file. (Not required)
```

