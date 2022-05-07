const workflow2 = function(refalign,
                           outputdir     = "./", 
                           num_threads   = 8, 
                           disable_cache = FALSE) {

    # generate file: ./16s.files
    const getSample_id = Metagenomics::mothur_files(getwd(), file.path(outputdir, "16s.files"));
    const check_filecache as function(filename) {
        if (disable_cache) {
            FALSE;
        } else {
            file.exists(filename);
        }
    }

    setwd(outputdir);

    if (!check_filecache("16s.trim.contigs.fasta")) {
        # make.contigs
        runMothur(
            command = "make.contigs",
            argv    = list(
                file       = "16s.files", 
                processors = num_threads
            ),
            log     = "[1]make.contigs.txt"
        );

        screen.seqs(
            contigs="16s.trim.contigs.fasta", 
            num_threads=num_threads
        );
    }

    if (!check_filecache("16s.trim.contigs.good.align")) {
        align.seqs(
            contigs="16s.trim.contigs.good.fasta", 
            reference=refalign$template,
            num_threads = num_threads
        );
    }

    if (!check_filecache("16s.trim.contigs.good.good.align")) {
        runMothur(
            command = "screen.seqs",
            argv    = list(
                fasta     = "16s.trim.contigs.good.align",
                alignreport="16s.trim.contigs.good.align.report", 
                minsim=90, 
                minscore=10, 
                group="16s.contigs.good.groups"
            ),
            log     = "[3]screen.seqs.txt"
        );
    }

    if (!check_filecache("16s.trim.contigs.good.good.gg.knn.taxonomy")) {
        classify.seqs(
            fasta="16s.trim.contigs.good.good.align", 
            reference=refalign$greengenes, 
            taxonomy=refalign$taxonomy, 
            method="knn", 
            processors=num_threads,
            numwanted=3
        );
    }

    if (!check_filecache("16s.trim.contigs.good.good.gg.knn.tax.summary")) {
        summary.tax(
            taxonomy="16s.trim.contigs.good.good.gg.knn.taxonomy", 
            group="16s.contigs.good.good.groups"
        );

        split.groups(
            fasta="16s.trim.contigs.good.fasta", 
            group="16s.contigs.good.groups"
        );
    }

    # result file
    # 16s.trim.contigs.good.good.gg.knn.tax.summary
    #  -> 16s_results.summary

    const file = "./16s.trim.contigs.good.good.gg.knn.tax.summary";
    const tree = taxonomy_kit::read.mothurTree(file);
    const OTU = tree 
    |> as.OTU_table() 
    |> as.data.frame()
    ;

    if ((length(getSample_id) == 1) && ("total" in colnames(OTU))) {
        v = OTU[, "total"];
        OTU[, "total"] = NULL;
        OTU[, getSample_id] = v;
    }

    print("parse mothur result file:");
    print(normalizePath(file));
    print("get the final OTU table result:");
    str(OTU);

    write.csv(OTU, file = "./mothur_OTU_table.csv", row.names = FALSE);
}