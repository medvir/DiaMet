<!-- README.md is generated from README.Rmd. Please edit that file -->

# DiaMet

This repository aims to analyze those metagenomic Illumina reads that
were unable to be classified by
[VirMet](https://github.com/medvir/VirMet)
(`undetermined_reads.fastq.gz`) by aligning them on protein level using
[DIAMOND](https://github.com/bbuchfink/diamond)/BLASTx.

`diamet.py` analyzes all reads from `undetermined_reads.fastq.gz`

- as single reads, and

- as contigs created by *de novo* assembly using
  [megahit](https://github.com/voutcn/megahit).

### How to run

1.  Enter timavo.

`ssh timavo`

2.  Move into the directory of the sample whose undetermined reads you
    want to analyze.

`cd /analysis/VirMet/<run>/<sample>/`

3.  Run the python script.

`python <path to script>/diamet.py`

### Input

To run `diamet.py`, you need:

- the `diamond` unix executable file which can be found
  [here](https://github.com/bbuchfink/diamond);

- [megahit](https://github.com/voutcn/megahit) installed on the server;

- a protein database (defined in the code; we are using *swissprot*);

- `undetermined_reads.fastq.gz`, which should be in the current working
  directory.

### Output

`diamet.py` will output the following files:

- `undetermined_reads_diamet.pdf` which plots taxonomic classification
  distribution of all hits;

- `undetermined_reads_diamet.tsv` which lists all hits and their Query
  Seq - id (qseqid), Query sequence length (qlen), Alignment length
  (length), Unique Subject Scientific Name (sscinames), and Unique
  Subject Super Kingdom (sskingdoms);

- `undetermined_reads_diamet_viral.csv` which lists only the viral hits
  and their counts;

- `undetermined_contigs_diamet.tsv` which lists all hits of the contigs
  and their Query Seq - id (qseqid), Query sequence length
  (qlen),Alignment length (length), Unique Subject Scientific Name
  (sscinames), and Unique Subject Super Kingdom (sskingdoms).
