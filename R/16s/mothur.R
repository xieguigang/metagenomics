#' The ``make.contigs`` command reads a forward fastq file and a reverse
#' fastq file and outputs new fasta and report files.
#'
#' @details The make.contigs command parameters are file, ffastq, rfastq,
#' ffasta, rfasta, fqfile, rqfile, findex, rindex, oligos, format, tdiffs,
#' bdiffs, pdiffs, align, match, mismatch, gapopen, gapextend, insert,
#' deltaq, maxee, allfiles and processors.
#'
#' Here we create a combo file for given input data to mothur: the 3 column
#' format is used for datasets where the sequences have already had the
#' barcodes and primers removed and been split into separate files.
#'
#' The first column is the group, the second is the forward fastq and
#' the third column contains the reverse fastq.
#'
#' @param algorithm [align] The align parameter allows you to specify the 
#'     alignment method to use. Your options are: gotoh and needleman. 
#'     The default is needleman.
#' @param score [match & mismatch & gapopen & gapextend] These parameters 
#'    are used while aligning the sequence read to determine the overlap. 
#'    The match parameter allows you to specify the bonus for having the 
#'    same base. The default is 1.0. The mistmatch parameter allows you to 
#'    specify the penalty for having different bases. The default is -1.0. 
#'    The gapopen parameter allows you to specify the penalty for opening a 
#'    gap in an alignment. The default is -2.0. The gapextend parameter 
#'    allows you to specify the penalty for extending a gap in an alignment. 
#'    The default is -1.0.
#' 
#' @param insert The insert parameter allows you to set a quality scores 
#'    threshold. When we are merging the overlapping regions, in the case 
#'    where we are trying to decide whether to keep a base or remove it 
#'    because the base is compared to a gap in the other fragment, if the 
#'    base has a quality score below the threshold we eliminate it. 
#'    Default=20.
#' 
const make.contigs = function(left, right, 
                              algorithm   = ["needleman", "gotoh"], 
                              score       = [1.0, -1.0, -2.0, -1.0],
                              insert      = 20,  
                              num_threads = 8) {

    # write raw data sequence file inputs
    # at here
    cat(`16s\t${left}\t${right}`, file = "./16s.files");

    # run mothur make.contigs command
    runMothur(
        command = "make.contigs",
        argv    = list(
            file        = "16s.files", 
            checkorient = "T", 
            align       = algorithm[1],
            match       = score[1],
            mismatch    = score[2],
            gapopen     = score[3],
            gapextend   = score[4],
            insert      = insert,
            processors  = num_threads
        ),
        log     = "[1]make.contigs.txt"
    );

    # result content files should be appears in the
    # output directory:

    # -rw-r--r-- 1 root root  3256536 Dec 10 17:24 16s.contigs.groups
    # -rw-r--r-- 1 root root  4341698 Dec 10 17:24 16s.contigs.report
    # -rw-r--r-- 1 root root       49 Dec 10 17:28 16s.files
    # -rw-r--r-- 1 root root        0 Dec 10 17:23 16s.scrap.contigs.fasta
    # -rw-r--r-- 1 root root        0 Dec 10 17:23 16s.scrap.contigs.qual
    # -rw-r--r-- 1 root root 31614742 Dec 10 17:24 16s.trim.contigs.fasta
    # -rw-r--r-- 1 root root 86838453 Dec 10 17:24 16s.trim.contigs.qual
}

#' Contig file renames
#' 
#' @param seqfile the source sequence file
#' @param contig.fasta the target renamed sequence file.
#' 
#' @return the parameter value comes from the ``contig.fasta`` parameter.
#'
const write_contig = function(seqfile, contig.fasta = "contig.fasta") {
    if (file.exists(contig.fasta)) {
        print(`${contig.fasta} is already exists, going to remove it...`);
        unlink(contig.fasta);

        if (file.exists(contig.fasta)) {
            print("unlink of target file is not working!");
        } else {
            print("yes!");
        }
    }

    print(`${seqfile} => ${contig.fasta}...`);

    seqfile
    |> readLines()
    |> writeLines(con = contig.fasta)
    ;
    cat("done~!\n");

    contig.fasta;
}

#' The ``summary.seqs`` command will summarize the quality of sequences in
#' an unaligned or aligned fasta-formatted sequence file.
#'
#' @param num_threads The processors option enables you to accelerate the
#' summary process by using multiple processors. Default processors=Autodetect
#' number of available processors and use all available.
#'
const summary.seqs = function(seqfile, count_table,
                              num_threads = 8,
                              logfile     = "log.txt") {

    runMothur(
        command = "summary.seqs",
        argv    = list(fasta=seqfile,processors=num_threads),
        log     = logfile
    );
}

#' RunAutoScreen
#'
#' @description The screen.seqs command enables you to keep sequences that
#' fulfill certain user defined criteria. Furthermore, it enables you to
#' cull those sequences not meeting the criteria from a names, group,
#' contigsreport, alignreport and summary file.
#'
const screen.seqs = function(contigs, num_threads = 8) {
    using workdir as workdir(dirname(contigs)) {
        const file as string = "[2]summary.seqs.txt";
        const summary = list(min = 0, max = 0);
        const groups  = "16s.contigs.groups";

        # sftp://xieguigang:biodeep123@192.168.0.207/home/xieguigang/16s/test/16s_result/[2]summary.seqs.txt
        # sftp://xieguigang:biodeep123@192.168.0.207/home/xieguigang/16s/test/16s_result/contig.summary
        summary.seqs(
            seqfile     = contigs,
            count_table = NULL,
            num_threads = num_threads,
            logfile     = file
        );

        print("sequence summary result:");
        cat("\n");
        cat(readText(file));
        cat("\n\n");

        # filter out all empty string line:
        # argument is of length zero
        for(line in readLines(file) |> which(l -> l != "")) {
            const cols as string = unlist(strsplit(line, "\t", fixed = FALSE));
            const header as string = cols[1];

            if (header == "2.5%-tile:") {
                summary$min = cols[4];
            } else if (header == "97.5%-tile:") {
                summary$max = cols[4];
            }
        }

        print("Parsed result:");
        str(summary);

        runMothur(
            command = "screen.seqs",
            argv    = list(
                fasta     = contigs,
                group     = groups,
                maxambig  = 0,
                minlength = summary$min,
                maxlength = summary$max
            ),
            log     = "[3]screen.seqs.txt"
        );
        # contig.good.fasta
        # contig.bad.accnos
        # 16s.contigs.good.groups
    }
}

#' The unique.seqs command returns only the unique sequences found in
#' a fasta-formatted sequence file and a file that indicates those
#' sequences that are identical to the reference sequence. Often times
#' a collection of sequences will have a significant number of identical
#' sequences. It sucks up considerable processing time to have to align,
#' calculate distances, and cluster each of these sequences individually.
#'
const unique.seqs = function(contigs, logfile = "[4]unique.seqs.txt") {
    using workdir as workdir(dirname(contigs)) {
        runMothur(
            command = "unique.seqs",
            argv    = list(fasta = contigs),
            log     = logfile
        );
        # contig.good.names
        # contig.good.unique.fasta
    }
}

#' The count.seqs / make.table command counts the number of sequences
#' represented by the representative sequence in a name file. If a group
#' file is given, it will also provide the group count breakdown.
#'
const count.seqs = function(names, groups) {
    runMothur(
        command = "count.seqs",
        argv    = list(name = names, group = groups),
        log     = "[5]count.seqs.txt"
    );
    # contig.good.count_table
}

#' The align.seqs command aligns a user-supplied fasta-formatted candidate
#' sequence file to a user-supplied fasta-formatted template alignment.
#' The general approach is to i) find the closest template for each candidate
#' using kmer searching, blastn, or suffix tree searching; ii) to make a
#' pairwise alignment between the candidate and de-gapped template sequences
#' using the Needleman-Wunsch, Gotoh, or blastn algorithms; and iii) to re-insert
#' gaps to the candidate and template pairwise alignments using the NAST
#' algorithm so that the candidate sequence alignment is compatible with the
#' original template alignment. We provide several alignment databases for
#' 16S and 18S rRNA gene sequences that are compatible with the greengenes or
#' silva alignments; however, custom alignments for any DNA sequence can be
#' used as a template and users are encouraged to share their alignment for
#' others to use. In general the alignment is very fast - we are able to align
#' over 186,000 full-length sequences to the SILVA alignment in less than 3 hrs
#' with a quality as good as the SINA aligner. Furthermore, this rate can be
#' accelerated using multiple processors. While the aligner doesn’t explicitly
#' take into account the secondary structure of the 16S rRNA gene, if the
#' template database is based on the secondary structure, then the resulting
#' alignment will at least be implicitly based on the secondary structure.
#'
const align.seqs = function(contigs, reference,
                            flip        = "F",
                            num_threads = 8,
                            logfile     = "[7]align.seqs.txt") {
    runMothur(
        command = "align.seqs",
        argv    = list(
            candidate  = contigs,
            template   = reference,
            flip       = flip,
            processors = num_threads
        ),
        log = logfile
    );
}

#' filter.seqs removes columns from alignments based on a criteria defined by the
#' user. For example, alignments generated against reference alignments (e.g. from
#' RDP, SILVA, or greengenes) often have columns where every character is either a
#' ‘.’ or a ‘-‘. These columns are not included in calculating distances because
#' they have no information in them. By removing these columns, the calculation of
#' a large number of distances is accelerated. Also, people also like to mask their
#' sequences to remove variable regions using a soft or hard mask (e.g. Lane’s mask).
#' This type of masking is only encouraged for deep-level phylogenetic analysis, not
#' fine level analysis such as that needed with calculating OTUs.
#'
const filter.seqs = function(align, num_threads, logfile = "[8]filter.seqs.txt") {
    runMothur(
        command = "filter.seqs",
        argv    = list(
            fasta      = align,
            processors = num_threads
        ),
        log   = logfile
    );
}

#' create sequence distance matrix
#' 
#' @details The dist.seqs command will calculate uncorrected pairwise 
#' distances between aligned DNA sequences. This approach is better 
#' than the commonly used DNADIST because the distances are not stored 
#' in RAM, rather they are printed directly to a file.
#' 
#' Furthermore, it is possible to ignore “large” distances that one 
#' might not be interested in. The command will generate a column-formatted 
#' distance matrix that is compatible with the column option in the 
#' various cluster commands. The command is also able to generate a 
#' phylip-formatted distance matrix.
#'
const dist.seqs as function(align_fasta,
                            calc        = "onegap",
                            countends   = "F",
                            cutoff      = 0.03,
                            output      = "lt",
                            num_threads = 32,
                            logfile     = "[10]dist.seqs.txt" ) {
    runMothur(
        command = "dist.seqs",
        argv    = list(
            fasta      = align_fasta,
            calc       = calc,
            countends  = countends,
            cutoff     = cutoff,
            output     = output,
            processors = num_threads
        ),
        log  = logfile
    );
}

#' assign sequences to OTUs
#'
#' @description Once a distance matrix gets read into mothur, the cluster command can be
#' used to assign sequences to OTUs. Presently, mothur implements three clustering
#' methods:
#'
#' + OptiClust (opti): OTUs are assembled using metrics to determine the quality of
#'                     clustering (the default setting).
#' + Nearest neighbor (nearest): Each of the sequences within an OTU are at most X% distant
#'                               from the most similar sequence in the OTU.
#' + Furthest neighbor (furthest): All of the sequences within an OTU are at most X%
#'                                 distant from all of the other sequences within the OTU.
#' + Average neighbor (average): This method is a middle ground between the other two
#'                               algorithms.
#' + AGC (agc): Abundance-based greedy clustering.
#' + DGC (dgc): Distance-based greedy clustering.
#' + Unique (unique): Creates a list file from a name or count file where every unique
#'                    sequence is assigned to it’s own OTU (i.e. Amplicon Sequence
#'                    Variant)
#'
#' If there is an algorithm that you would like to see implemented, please consider either
#' contributing to the mothur project or contacting the developers and we’ll see what we
#' can do. The opticlust algorithm is the default option.
#'
const cluster as function(dist,
                          method      = "furthest",
                          cutoff      = 0.03,
                          num_threads = 32,
                          logfile     = "[11]cluster.txt") {
    runMothur(
        command = "cluster",
        argv = list(phylip = dist,
                    method = method,
                    cutoff = cutoff,
                    processors = num_threads
        ),
        log  = logfile
    );
}

#' 
#' 
#' @details ``bin.seqs`` prints out a fasta-formatted file where 
#' sequences are ordered according to the OTU that they belong to.
#' Such an output may be helpful for generating primers specific 
#' to an OTU or for classification of sequences.
#'
const bin.seqs = function(list, contigs, contig.names, logfile = "[12]bin.seqs.txt") {
    runMothur(
        command = "bin.seqs",
        argv    = list(list = list, fasta = contigs, name = contig.names),
        log     = logfile
    );
}

#' Get the representative sequence
#' 
#' @details While the ``bin.seqs`` command reports the OTU number for 
#' all sequences, the ``get.oturep`` command generates a fasta-formatted 
#' sequence file containing only a representative sequence for each OTU. 
#' A .rep.fasta and .rep.name file or .rep.count_table file is generated
#' for each OTU definition.
#'
const get.oturep as function(dist, contig.unique.fasta, list,
                             label   = 0.03,
                             logfile = "[13]get.oturep.txt") {
    runMothur(
        command = "get.oturep",
        argv    = list(
            phylip = dist,
            fasta  = contig.unique.fasta,
            list   = list,
            label  = label
        ),
        log     = logfile
    );
}

#' Do OTU annotation based on the reference database
#' 
#' @details The classify.seqs command allows the user to use 
#' several different methods to assign their sequences 
#' to the taxonomy outline of their choice. Current 
#' methods include the Wang approach, using a k-nearest 
#' neighbor consensus and zap. Taxonomy outlines and 
#' reference sequences can be obtained from the 
#' taxonomy outline page. The command requires that 
#' you provide a fasta-formatted input and database 
#' sequence file and a taxonomy file for the reference 
#' sequences.
#' 
const classify.seqs = function(
    fasta, 
    reference, 
    taxonomy, 
    method="knn", 
    numwanted=3,
    processors=2,
    logfile = "[14]classify.seqs.txt"  
) {
    runMothur(
        command = "classify.seqs",
        argv    = list(
            reference = reference,
            fasta  = fasta,
            taxonomy   = taxonomy,
            method  = method,
            numwanted = numwanted,
            processors = processors
        ),
        log     = logfile
    );
}

#' taxonomy summary of the annotation
#' 
#' @details The summary.tax command reads a taxonomy file and 
#' an optional name and or group file, and summarizes 
#' the taxonomy information.
#' 
const summary.tax = function(taxonomy, group, logfile = "[15]summary.tax.txt") {
    runMothur(
        command = "summary.tax",
        argv    = list(
            taxonomy = taxonomy,
            group = group
        ),
        log     = logfile
    );
}

#' count for groups
#' 
#' @details The split.groups command reads a list, fasta, flow or 
#' fastq file and count and generates a list, fasta, flow 
#' or fastq files for each group in the count file.
#' 
const split.groups = function(fasta, group, logfile = "[16]split.groups.txt") {
    runMothur(
        command = "split.groups",
        argv    = list(
            fasta = fasta,
            group = group
        ),
        log     = logfile
    );
}