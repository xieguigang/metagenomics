imports "gast" from "metagenomics_kit";

const align_silva as function(blastn = getOption("ncbi_blast")) {

}

#' Align greengenes database for taxonomy annotation
#' 
#' @param work16s the workspace object which is comes 
#'      from the ``greengenes_OTUTaxonomy`` function.
#' @param blastn the executable file path of the NCBI 
#'      blastn program.
#' @param is_debug run in debug mode?
#' 
const align_greengenes as function(work16s, 
                                   blastn   = getOption("ncbi_blast"), 
                                   is_debug = getOption("workflow.debug")) {

    # 进行blastn序列比对操作来完成物种鉴定
    is_debug = as.logical(is_debug);
    blastn_cli = `
        ${blastn} 
            -query  "${work16s$OTU_rep}" 
            -db     "${work16s$greengenes[[1]]}" 
            -out    "${work16s$blastnOut}" 
            -evalue "${work16s$evalue}"
            -num_threads ${work16s$num_threads}
    `;
    
    print("commandline for run taxonomy annotation:");
    print(blastn_cli);    

    if (!is_debug) {
        blastn_cli |> system(intern = TRUE);    
    } else {
        print("skip of run blastn in debug model!");
    }

    if (file.exists(work16s$blastnOut)) {
        # run blastn parser
        csv = `${dirname(work16s$blastnOut)}/${basename(work16s$blastnOut)}.csv`;

        using buffer as open.stream(csv, type = "Mapping", ioRead = FALSE) {
            work16s$blastnOut 
            |> getBlastnMapping()
            |> stream.flush(stream = buffer)
            ;
        }

        if (file.exists(csv)) {
            csv 
            |> read.csv() 
            |> print(max.print = 10)
            ;
        }

        # run taxonomy annotation and create OTU 
        # relative abundance table result.
        gast = work16s |> taxonomy_annotation();

        write.csv(gast, file = `${outputdir}/gast_OTUs.csv`);
    } else {
        stop("error during of run taxonomy annotation process.");
    }
}

#' Do taxonomy annotation for the OTU contigs
#' 
#' @param OTU_mapping the blastn search result table.
#' 
#' @return the 16s raw data processing pipeline final result.
#' 
const taxonomy_annotation as function(work16s) {
    taxonomy = parse.greengenes_tax(work16s$greengenes$taxonomy);  
    OTU_rep = work16s$OTU_rep |> parse.mothur_OTUs();
    gast = work16s$blastnOut 
        |> read.blast(type = "nucl") 
        |> OTU.taxonomy(
            OTUs           = OTU_rep, 
            taxonomy       = taxonomy, 
            gast.consensus = FALSE
        )
        ;

    gast;
}