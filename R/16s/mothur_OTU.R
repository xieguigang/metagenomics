
#' Build OTU based on mothur
#' 
#' @param left the left paired-end sequence file.
#' @param right the right paired-end sequence file.
#' @param outputdir the result output dir for save temp file
#' @param num_threads the number of the processor threads that 
#'    will be used in the mothur software.
#' 
const mothur_OTU as function(left, right, refalign,
                             outputdir   = "./", 
                             num_threads = 32) {

    work16s = list(
        left        = left, 
        right       = right, 
        outputdir   = outputdir, 
        num_threads = num_threads,
        silva       = refalign
    );

    using workdir as workdir(outputdir) {
        work16s = work16s |> mothur_workflow();
    }

    print("Mothur job done!");

    return(work16s);
}

#' commandline workflow for run mothur
#' 
#' @param work16s the 16s workflow workspace object, contains necessary
#'    parameter files.
#' 
#' @details you must change of current work directory to the output 
#'    directory, due to the reason of all of the data file path in 
#'    this workflow function is relative path to the result data directory
#'    ``outputdir``.
#' 
const mothur_workflow as function(work16s) {
    print("View of the mothur workflow parameters:");
    str(work16s);

    print("current workspace for run mothur:");
    print(getwd());

    cat("\n\n");
    cat("-----------------------------=======================================================================-----------------------------------\n");
    cat("Schloss, P.D., et al.,\n");
    cat("Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities.\n");
    cat("Appl Environ Microbiol, 2009. 75(23):7537-41.\n");
    cat("\n");

    make.contigs(work16s$left, work16s$right, num_threads = work16s$num_threads);
    work16s$contig.fasta = write_contig("16s.trim.contigs.fasta");

    # RunAutoScreen
    # summary.seqs + screen.seqs 
    # 
    work16s$contig.fasta |> screen.seqs(num_threads = work16s$num_threads);
    work16s$contigs = "contig.good.fasta";
    work16s$groups  = "16s.contigs.good.groups";

    unique.seqs(work16s$contigs, logfile = "[4]unique.seqs.txt");

    work16s$names = "contig.good.names"; 
    work16s$names |> count.seqs(groups= work16s$groups);

    work16s$count_table = "contig.good.count_table";
    summary.seqs(
        seqfile     = "contig.good.unique.fasta", 
        count_table = work16s$count_table,
        num_threads = work16s$num_threads, 
        logfile     = "[6]summary.seqs.txt"
    );

    # contig.good.unique.summary
    write_contig("contig.good.unique.fasta");

    work16s$contigs = "contig.fasta";
    work16s$contigs |> align.seqs(        
        reference   = work16s$silva,
        flip        = "T",
        num_threads = work16s$num_threads, 
        logfile     = "[7]align.seqs.txt"
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
        contigs      = work16s$contigs,
        contig.names = "contig.names",
        logfile      = "[12]bin.seqs.txt"
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

    work16s;
}