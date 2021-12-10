imports "annotation.workflow" from "seqtoolkit";

require(GCModeller);

# run blast for gene annotation

#' Parse blastn result text
#' 
#' @param outTxt the filepath of the blastn alignment output.
#'    function will parse the blastn mapping table result from
#'    this data file via GCModeller library.
#' 
const getBlastnMapping as function(outTxt) {
    outTxt 
    |> read.blast(type = "nucl")
    |> blastn.maphit()
    ;
}