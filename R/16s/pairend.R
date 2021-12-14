const pairEnds as function(dir) {
    files = dir
    |> list.files(pattern = "*.fq")
    |> groupBy(function(path) {
        basename(path) 
        |> gsub("_R1", "") 
        |> gsub("_R2", "")
    })
    |> lapply(function(pairend) {      
        list(
            R1 = pairend[pairend like $".+_R1.fq"],
            R2 = pairend[pairend like $".+_R2.fq"]
        );
    }, names = p -> p$key);

    R1 = sapply(files, pair -> pair$R1);
    R2 = sapply(files, pair -> pair$R2);

    data.frame(
        sample_id = names(files), 
        R1 = basename(R1), 
        R2 = basename(R2), 
        dir = dir
    );
}

#' Create mothur file list
#' 
#' @param dir a directory path which contains the 
#'    ``R1/R2`` pair-end fastq sequence file.
#' 
const mothur_files as function(dir, filename = "16s.files") {
    print("scan the fastq raw data files...");

    # load pair-end file list
    files = pairEnds(dir);
    # write raw data sequence file inputs
    # at here
    print(files);
    
    [sample_id, R1, R2] = files;

    sample_id = make.names(sample_id);
    (1:length(sample_id))
    |> sapply(function(i) {
        sprintf("%s\t%s%s.fq\t%s%s.fq", sample_id[i], dir, R1[i], dir, R2[i]);
    })
    |> writeLines(con = file.path(dir, filename))
    ;

    invisible(NULL);
}