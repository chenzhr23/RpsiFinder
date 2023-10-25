\name{RpsiFinder}
\alias{RpsiFinder}
\title{R API for psiFinder}
\usage{
RpsiFinder <- function(genome_fa,
                       genome_fai,
                       treat_bam,
                       input_bam,
                       gene_bed12,
                       output_dir,
                       output_name,
                       options = "-p 1.5 -t 5 -r 0.05 -M 1 -f 1 -m 0 -s -n -w 20")
}
\arguments{
  \item{genome_fa}{fasta file of reference genome, usually downloaded form databases like UCSC, ENSEMBL, etc.}
  \item{genome_fai}{fai file of reference genome, usually obtained by samtools, etc.}
  \item{treat_bam}{bam file of treatment sample, usually derived from CMC-treatment protocol, etc.}
  \item{inpt_bam}{bam file of input sample, usually derived from CMC-control protocol, etc.}
  \item{gene_bed12}{bed12 file of annotation, usually downloaded form databases like UCSC, ENSEMBL, etc., or organize the file by bioinformatics tools}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
  \item{options}{options for psiFinder:
    -v/--verbose                   : verbose information
    -V/--version                   : rtsSeeker version
    -h/--help                      : help informations
    -U/--psi                       : search psi sites
    -s/--strand                    : strand-specific[default=false]
    -S/--soft                      : skip read with soft-clip[default=false]
    -n/--norm                      : normalized reads to the locus number
    -p/--pval                      : p-value[default=1.5]
    -c/--collapser                 : keep duplication, deault is false
    -M/--model <int>               : model[0 for genome, 1 for gene, 2 for both, default=1]
    --treat <string>               : treatment file<BAM format>
    --input <string>               : input file<BAM format>
    --gene <string>                : gene file <BED12 format>
    -o/--outfile <string>          : output file
    -t/--min-tag <double>          : minimum tag number for each psi, default>=5.0 read
    -r/--rpm <double>              : minimum rpm value for each psi, default>=0.05
    -f/--fold <int>                : minimum fold-change[default>=1.0]
    -m/--mfold <int>               : minimum fold-change for stop-tag/mean-tag
    -w/--window <int>              : window size around the rts position[default=50]
    -l/--min-len <int>             : minimum length of reads, default=15
    -L/--max-len <int>             : maximum length of reads, default=1000}
}
\description{Run psiFinder for pseudouridylation identification.}
\examples{
RpsiFinder(genome_fa="../test_data/genome/hg38.fa",
            genome_fai="../test_data/genome/hg38.fa.fai",
            treat_bam="../test_data/A2.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gzAligned.sortedByCoord.out.bam",
            input_bam = "../test_data/A1.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gzAligned.sortedByCoord.out.bam",
            gene_bed12 = "../test_data/annotation/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12",
            output_dir="/home/chenzhr23/test_out/RpsiFinder_out",
            output_name = "RpsiFinder_out")
}
\value{tab-delimited text files include stop read information for each CMC-treatment and CMC-control group, e.g. treatStopRatio, ctrlStopRatio}
