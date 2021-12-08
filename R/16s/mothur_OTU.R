
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
    screen.seqs(work16s$contig.fasta, num_threads = num_threads);

    work16s$contigs = "contig.good.fasta";
    work16s$groups  = "16s.contigs.good.groups";

    unique.seqs(work16s$contigs, logfile = "[4]unique.seqs.txt");

    work16s$names = "contig.good.names"; 
    count.seqs(work16s$names, groups= work16s$groups);

    work16s$count_table = "contig.good.count_table";
    summary.seqs(
        seqfile="contig.good.unique.fasta", 
        count_table=work16s$count_table,
        num_threads = num_threads, 
        logfile = "[6]summary.seqs.txt"
    );
    # contig.good.unique.summary

    write_contig("contig.good.unique.fasta");

    work16s$contigs = "contig.fasta";

    
    align.seqs(fasta=work16s$contigs,reference=work16s$silva,flip="T",processors=work16s$num_threads, logfile = "[7]align.seqs.txt");
    # contig.align
    # contig.align.report
    # contig.flip.accnos

    work16s$align = "contig.align";
    filter.seqs(fasta=work16s$align,processors=work16s$num_threads, logfile = "[8]filter.seqs.txt");
    # contig.filter
    # contig.filter.fasta

    write_contig("contig.filter.fasta");
    unique.seqs(fasta=work16s$contigs, logfile = "[9]unique.seqs.txt");
    # contig.names
    # contig.unique.fasta
    
    work16s$align = "contig.unique.fasta";
    dist.seqs(fasta=work16s$align,calc="onegap",countends="F",cutoff=0.03,output="lt",processors=work16s$num_threads, logfile = "[10]dist.seqs.txt");
    # contig.unique.phylip.dist

    work16s$dist = "contig.unique.phylip.dist";
    cluster(phylip=work16s$dist,method="furthest",cutoff=0.03,processors=work16s$num_threads, logfile = "[11]cluster.txt");

    # contig.unique.phylip.fn.sabund
    # contig.unique.phylip.fn.rabund
    # contig.unique.phylip.fn.list

    work16s$list = "contig.unique.phylip.fn.list";
    bin.seqs(list=work16s$list,fasta=work16s$contigs,name="contig.names", logfile = "[12]bin.seqs.txt");

    # contig.unique.phylip.fn.unique.fasta
    # contig.unique.phylip.fn.0.01.fasta
    # contig.unique.phylip.fn.0.02.fasta
    # contig.unique.phylip.fn.0.03.fasta
    get.oturep(phylip=work16s$dist,fasta="contig.unique.fasta",list=work16s$list,label=0.03, logfile = "[13]get.oturep.txt");

    print("Mothur job done!");
}

