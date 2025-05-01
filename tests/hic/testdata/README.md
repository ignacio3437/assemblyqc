# Test data for HiC sub workflow

- `HYh1h2_Chr01_6000000_6015000.fasta.gz:` 15000 bases from chr1 of each haplotype starting at 6 MB
- `HYh1h2_Chr01_6000000_6015000.R{1,2}.fastq.gz:` Raw HiC data was downloaded from SRA, passed through `fastp` with args `--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50`, and mapped to the original full assembly. The BAM file was then sliced to the chr1 region, and reads were extracted
- `HYh1h2_Chr01_6000000_6015000.bam:` Obtained by mapping `fastq` to `fasta` using bwa-mem

## Data sources

- `fasta:` <kiwifruitgenome.atcgn.com>
- `fastq:` SRR21031552
