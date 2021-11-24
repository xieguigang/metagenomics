#' the mothur cli wrapper
#' 
#' @param command the mothur command string.
#' @param log the logfile file path.
#'
const runMothur as function(command, argv, log) {
    const mothur as string      = getOption("mothur");
    const is_debug as boolean   = getOption("workflow.debug");
    const commandArgs as string = `${mothur} "#${command}(${mothurArgvs(argv)})"`;

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

const mothurArgvs as function(argv) {
    argv 
    |> names()
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
} 

const readSummary as function(contigs, num_threads = 8) {
    using workdir as workdir(dirname(file)) {
        const file as string = "[2]summary.seqs.txt";
        const summary = list(min = 0, max = 0);
        const groups  = "16s.contigs.groups";

        # RunAutoScreen
        runMothur(
            command = "summary.seqs",
            argv    = list(fasta=contigs,processors=num_threads), 
            log     = file
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

const unique.seqs as function(contigs) {
    using workdir as workdir(dirname(contigs)) {
        runMothur(
            command = "unique.seqs",
            argv    = list(fasta=$contigs), 
            log     = "[4]unique.seqs.txt"
        );
        # contig.good.names
        # contig.good.unique.fasta
    }
}