
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
    cat("-------------===============================================--------------\n");
    cat("Schloss, P.D., et al.,\n");
    cat("Introducing mothur: Open-source, platform-independent, community-supported\n");
    cat("software for describing and comparing microbial communities.\n");
    cat("Appl Environ Microbiol, 2009. 75(23):7537-41.\n");
    cat("\n");

    # make.contigs
    # 16s_result/16s.files
    # 16s_result/16s.scrap.contigs.fasta
    # 16s_result/16s.contigs.groups
    # 16s_result/16s.contigs.report
    # 16s_result/16s.trim.contigs.fasta
    # 16s_result/[1]make.contigs.txt
    make.contigs(work16s$left, work16s$right, num_threads = work16s$num_threads);
    # result file renames to contig.fasta
    work16s$contig.fasta = write_contig("16s.trim.contigs.fasta");

    # RunAutoScreen
    # summary.seqs + screen.seqs
    #
    # 16s_result/contig.bad.accnos
    # 16s_result/contig.good.fasta
    # 16s_result/16s.contigs.good.groups
    # 16s_result/[3]screen.seqs.txt
    work16s$contig.fasta |> screen.seqs(num_threads = work16s$num_threads);
    work16s$contigs = "contig.good.fasta";
    work16s$groups  = "16s.contigs.good.groups";

    # 16s_result/[4]unique.seqs.txt
    # 16s_result/contig.good.names
    # 16s_result/contig.good.unique.fasta
    unique.seqs(work16s$contigs, logfile = "[4]unique.seqs.txt");

    # 16s_result/[5]count.seqs.txt
    # 16s_result/contig.good.count_table
    work16s$names = "contig.good.names";
    work16s$names |> count.seqs(groups= work16s$groups);

    # 16s_result/[6]summary.seqs.txt
    # 16s_result/contig.good.unique.summary
    work16s$count_table = "contig.good.count_table";
    summary.seqs(
        seqfile     = "contig.good.unique.fasta",
        count_table = work16s$count_table,
        num_threads = work16s$num_threads,
        logfile     = "[6]summary.seqs.txt"
    );

    # contig.good.unique.summary
    #  -> contig.fasta
    write_contig("contig.good.unique.fasta");

    # 16s_result/contig.align
    # 16s_result/[7]align.seqs.txt
    # 16s_result/contig.align.report
    # 16s_result/contig.flip.accnos
    work16s$contigs = "contig.fasta";
    work16s$contigs |> align.seqs(
        reference   = work16s$silva,
        flip        = "T",
        num_threads = work16s$num_threads,
        logfile     = "[7]align.seqs.txt"
    );

    # 16s_result/[8]filter.seqs.txt
    # 16s_result/contig.filter
    # 16s_result/contig.filter.fasta
    work16s$align = "contig.align";
    work16s$align |> filter.seqs(
        num_threads = work16s$num_threads,
        logfile     = "[8]filter.seqs.txt"
    );

    # contig.filter
    # contig.filter.fasta
    #  -> contig.fasta
    write_contig("contig.filter.fasta", work16s$contigs);
    # 16s_result/contig.unique.fasta
    # 16s_result/[9]unique.seqs.txt
    # 16s_result/contig.names
    unique.seqs(
        contigs = work16s$contigs,
        logfile = "[9]unique.seqs.txt"
    );

    # contig.names
    # contig.unique.fasta
    work16s$align = "contig.unique.fasta";
    # 16s_result/[10]dist.seqs.txt
    # 16s_result/contig.unique.phylip.dist
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
    # 16s_result/contig.unique.phylip.fn.list
    # 16s_result/contig.unique.phylip.fn.rabund
    # 16s_result/contig.unique.phylip.fn.sabund
    # 16s_result/[11]cluster.txt
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
    # 16s_result/contig.count_table
    # 16s_result/contig.unique.phylip.fn.0.01.pick.list
    # 16s_result/contig.unique.phylip.fn.unique.list
    # 16s_result/contig.unique.phylip.fn.0.02.pick.list
    # 16s_result/contig.unique.phylip.fn.0.03.pick.list
    # 16s_result/[12]bin.seqs.txt
    work16s$list |> bin.seqs(
        contigs      = work16s$contigs,
        contig.names = "contig.names",
        logfile      = "[12]bin.seqs.txt"
    );

    # sftp://xieguigang@192.168.0.207/home/xieguigang/16s/test/16s_result/contig.unique.phylip.fn.unique.unique.fasta
    # sftp://xieguigang@192.168.0.207/home/xieguigang/16s/test/16s_result/contig.unique.phylip.fn.unique.rep.fasta
    # sftp://xieguigang@192.168.0.207/home/xieguigang/16s/test/16s_result/contig.unique.phylip.fn.unique.rep.names
    # sftp://xieguigang@192.168.0.207/home/xieguigang/16s/test/16s_result/[13]get.oturep.txt
    work16s$dist |> get.oturep(
        contig.unique.fasta = "contig.unique.fasta",
        list    = work16s$list,
        label   = 0.03,
        logfile = "[13]get.oturep.txt"
    );

    work16s;
}
