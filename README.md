# metagenomics
metagenomics/16S data analysis

```

  '16s.R' - OTU analysis based on mothur

  this script apply for generate OTU contigs data from
  paired-end short reads files, which could be assembled via mothur
  program and run alignment with greengenes database or silva
  database.

SYNOPSIS
Rscript "16s.R" --left <*.fq> --right <*.fq> [--outputdir <directory, default="${dirname(&left)}/16s_result/">] [--num_threads <integer, default=72>] [--workflow-debug <boolean, default=FALSE>] [--mothur <string, default="/opt/metagenomics/mothur/mothur">] [--blast+ <string, default="/opt/metagenomics/ncbi-blast-2.12.0+/bin">] [--evalue <double, default=1E-10>] [--skip-mothur <boolean, default=FALSE>] [--greengenes <string, default="/opt/metagenomics/greengenes/taxonomy/gg_13_8_99.fasta+gg_13_8_99.gg.tax">] [--template <*.fasta, default="/opt/metagenomics/greengenes/refalign/gg_13_8_99.refalign">]

CommandLine Argument Values:

 --left:           the left paired-end short reads file.
 --right:          the right paired-end short reads file.
 --outputdir:      the directory folder path for output data result files.
 --num_threads:    an optional argument for setting up the processors that will be used
                   in mothur program, by default is using 72 processor threads.

 --workflow-debug: run debug of the workflow?
 --mothur:         file path to the mothur program executable file.
 --blast+:         the folder path of the NCBI blast+ suite for run OTU contig alignment
                   with the greengenes or silva database.

 --evalue:         the e-value cutoff of the blastn search pipeline. default e-value
                   is 1e-10.

 --skip-mothur:    skip of run mothur workflow for generate OTU sequence dataset? this
                   options is usually apply for run workflow debug.

 --greengenes:     the reference OTU sequnece database for run taxonomy annotation of
                   the OTU contigs which is generated from the mothur software. this
                   database file can be download from the mothur release page: https://mothur.org/wiki/greengenes-formatted_databases/

 --template:       the reference template file for run mothur alignment of the generated
                   contig fasta sequence file.


Authors:
[1]	xieguigang <xie.guigang@gcmodeller.org> 

Dependency List:

 Loading: 
[1]	Metagenomics 
```