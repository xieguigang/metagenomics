#' the mothur cli wrapper
#' 
#' @param command the mothur command string.
#' @param log the logfile file path.
#'
const runMothur as function(command, log) {
    const mothur as string      = getOption("mothur");
    const is_debug as boolean   = getOption("workflow.debug");
    const commandArgs as string = `${mothur} "#${command}" > ${log}`;

    print(commandArgs);

    if (!is_debug) {
        system(commandArgs);
    } else {
        0;
    }
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
            command = `make.contigs(file=16s.files, processors=${num_threads})`, 
            log     = `[1]make.contigs.txt`
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