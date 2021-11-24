require(Metagenomics);

#'
#' 

[@info "the left paired-end short reads file."]
[@type "*.fq"]
const left as string = ?"--left" || stop("At least one paired-end short reads file is required!");

[@info "the right paired-end short reads file."]
[@type "*.fq"]
const right as string = ?"--right";

[@info "the directory folder path for output data result files."]
[@type "directory"]
const outputdir as string = ?"--outputdir" || `${dirname(left)}/16s_result/`;

[@info "an optional argument for setting up the processors that will be 
        used in mothur program, by default is using 32 processor 
        threads."]
const num_threads as integer = ?"--num_threads" || 32;

