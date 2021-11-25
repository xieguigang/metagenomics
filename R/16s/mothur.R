#' the mothur cli wrapper
#' 
#' @param command the mothur command string.
#' @param log the logfile file path.
#'
const runMothur as function(command, argv, log) {
    const mothur as string      = getOption("mothur");
    const is_debug as boolean   = getOption("workflow.debug");
    const commandArgs as string = `${mothur} "#${command}(${mothurArgvs(command, argv)})"`;

    print(commandArgs);

    if (!is_debug) {
        commandArgs
        |> system(intern = TRUE)
        |> writeLines(con = log)
        ;
    } else {
        0;
    }
}

#' Create commandline arguments
#' 
#' @param argv a list object that contains commandline parameters for
#'    run a specific mothur command. commandline argument which it has 
#'    ``NULL`` value will be ignored in this function. 
#'
const mothurArgvs as function(argv) {
    for(name in names(argv)) {
        if (is.null(argv[[key]])) {
            warning(`missing value of '${name}' in '${command}'?`);
        }
    }

    argv 
    |> names()
    |> which(key -> !is.null(argv[[key]]))
    |> sapply(key -> `${key}=${argv[[key]]}`)
    |> paste(collapse = ", ")
    ;
}

#' @details The 3 column format is used for datasets where the sequences 
#' have already had the barcodes and primers removed and been split into 
#' separate files. The first column is the group, the second is the forward 
#' fastq and the third column contains the reverse fastq.
#' 
const make.contigs as function(left, right, outputdir, num_threads = 8) {
    using workdir as workdir(outputdir) {
        # write raw data sequence file inputs
        # at here
        cat(`16s\t${left}\t${right}`, file = `${outputdir}/16s.files`);
        
        # run mothur make.contigs command
        runMothur(
            command = "make.contigs",
            argv    = list(file="16s.files", processors=num_threads), 
            log     = "[1]make.contigs.txt"
        );
    }

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

const write_contig as function(seqfile) {
    const contig.fasta = "contig.fasta";

    if (file.exists(contig.fasta)) {
        unlink(contig.fasta);
    }

    print(`${seqfile} => ${contig.fasta}...`);
    file.rename(seqfile, contig.fasta);
    cat("done~!\n");

    contig.fasta;
} 

const summary.seqs as function(seqfile, count_table, 
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
const screen.seqs as function(contigs, num_threads = 8) {
    using workdir as workdir(dirname(file)) {
        const file as string = "[2]summary.seqs.txt";
        const summary = list(min = 0, max = 0);
        const groups  = "16s.contigs.groups";

        summary.seqs(
            seqfile     = contigs, 
            count_table = NULL, 
            num_threads = num_threads,
            logfile     = file
        );

        for(line in readLines(file)) {
            const cols   = unlist(strsplit(line, "\t", fixed = FALSE));
            const header = cols[1];

            if (header == "2.5%-tile:") {
                summary$min = cols[4];
            } else if (header == "97.5%-tile:") {
                summary$max = cols[4];
            }
        }

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

const unique.seqs as function(contigs, logfile = "[4]unique.seqs.txt") {
    using workdir as workdir(dirname(contigs)) {
        runMothur(
            command = "unique.seqs",
            argv    = list(fasta=contigs), 
            log     = logfile
        );
        # contig.good.names
        # contig.good.unique.fasta
    }
}

const count.seqs as function(names, groups) {
    runMothur(
        command = "count.seqs",
        argv    = list(name=names, group=groups), 
        log     = "[5]count.seqs.txt"
    );
    # contig.good.count_table
}

const align.seqs as function(contigs, silva, num_threads = 8) {
    runMothur(
        command = "align.seqs",
        argv    = list(
            fasta      = contigs,
            reference  = silva,
            flip       = "T",
            processors = num_threads
        ),
        log     = "[7]align.seqs.txt"
    );
}

const filter.seqs as function() {
    runMothur(
        command = "filter.seqs",
        argv    = list(
            fasta= align,
            processors= num_threads
        ),
        log   = "[8]filter.seqs.txt"
    );
}

const dist.seqs as function() {
    runMothur(
        command = "dist.seqs",
        argv    = list(
            fasta=$align,
            calc=onegap,
            countends=F,
            cutoff=0.03,
            output=lt,
            processors=$num_threads
        ),
        log  = "[10]dist.seqs.txt"
    );
}

const cluster as function() {
    runMothur(
        command = "dist.seqs",
        argv    = list(phylip=$dist,method=furthest,cutoff=0.03,processors=$num_threads),
        log     = "[11]cluster.txt"
    );
}

const bin.seqs as function() {
    runMothur(
        command = "bin.seqs",
        argv    = list(list=$list,fasta=$contigs,name=contig.names),
        log     = "[12]bin.seqs.txt"
    );
}

const get.oturep as function() {
    runMothur(
        command = "get.oturep",
        argv    = list(phylip=$dist,fasta=contig.unique.fasta,list=$list,label=0.03),
        log     = "[13]get.oturep.txt"
    );
}