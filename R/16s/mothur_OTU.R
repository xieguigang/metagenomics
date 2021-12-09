
#' Build OTU based on mothur
#' 
#' @param left the left paired-end sequence file.
#' @param right the right paired-end sequence file.
#' @param outputdir the result output dir for save temp file
#' @param num_threads the number of the processor threads that 
#'    will be used in the mothur software.
#' 
const mothur_OTU as function(left, right, 
                             outputdir   = "./", 
                             num_threads = 32) {

    const work16s = list(
        left        = left, 
        right       = right, 
        outputdir   = outputdir, 
        num_threads = num_threads
    );

    make.contigs(left, right, outputdir, num_threads = num_threads);
    work16s$contig.fasta = write_contig("16s.trim.contigs.fasta");

    # RunAutoScreen
    # summary.seqs + screen.seqs 
    # 
    work16s$contig.fasta |> screen.seqs(num_threads = num_threads);
    work16s$contigs = "contig.good.fasta";
    work16s$groups  = "16s.contigs.good.groups";

    unique.seqs(work16s$contigs, logfile = "[4]unique.seqs.txt");

    work16s$names = "contig.good.names"; 
    work16s$names |> count.seqs(groups= work16s$groups);

    work16s$count_table = "contig.good.count_table";
    summary.seqs(
        seqfile     = "contig.good.unique.fasta", 
        count_table = work16s$count_table,
        num_threads = num_threads, 
        logfile     = "[6]summary.seqs.txt"
    );
    
    # contig.good.unique.summary
    write_contig("contig.good.unique.fasta");

    work16s$contigs = "contig.fasta";
    work16s$contigs |> align.seqs(        
        silva      = work16s$silva,
        flip       = "T",
        processors = work16s$num_threads, 
        logfile    = "[7]align.seqs.txt"
    );

    # contig.align
    # contig.align.report
    # contig.flip.accnos
    work16s$align = "contig.align";
    work16s$align |> filter.seqs(       
        num_threads = work16s$num_threads, 
        logfile     = "[8]filter.seqs.txt"
    );

    # contig.filter
    # contig.filter.fasta
    write_contig("contig.filter.fasta");
    unique.seqs(
        contigs = work16s$contigs, 
        logfile = "[9]unique.seqs.txt"
    );

    # contig.names
    # contig.unique.fasta
    work16s$align = "contig.unique.fasta";
    work16s$align |> dist.seqs(
        calc        = "onegap",
        countends   = "F",
        cutoff      = 0.03,
        output      = "lt",
        num_threads = work16s$num_threads, 
        logfile     = "[10]dist.seqs.txt"
    );

    # contig.unique.phylip.dist
    work16s$dist = "contig.unique.phylip.dist";
    work16s$dist |> cluster(
        method      = "furthest",
        cutoff      = 0.03,
        num_threads = work16s$num_threads, 
        logfile     = "[11]cluster.txt"
    );

    # contig.unique.phylip.fn.sabund
    # contig.unique.phylip.fn.rabund
    # contig.unique.phylip.fn.list
    work16s$list = "contig.unique.phylip.fn.list";
    work16s$list |> bin.seqs(
        contigs = work16s$contigs,
        name    = "contig.names",
        logfile = "[12]bin.seqs.txt"
    );

    # contig.unique.phylip.fn.unique.fasta
    # contig.unique.phylip.fn.0.01.fasta
    # contig.unique.phylip.fn.0.02.fasta
    # contig.unique.phylip.fn.0.03.fasta
    work16s$dist |> get.oturep(
        contig.unique.fasta = "contig.unique.fasta",
        list    = work16s$list,
        label   = 0.03, 
        logfile = "[13]get.oturep.txt"
    );

    print("Mothur job done!");
}