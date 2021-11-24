require(Metagenomics);

#' title: OTU analysis based on mothur
#' author: xieguigang <xie.guigang@gcmodeller.org> 
#'
#' description: this script apply for generate OTU contigs data from
#'    paired-end short reads files, which could be assembled via mothur
#'    program and run alignment with greengenes database or silva 
#'    database.

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

[@info "run debug of the workflow?"]
const workflow_debug as boolean = ?"--workflow-debug";

[@info "file path to the mothur program executable file."]
const mothur as string = ?"--mothur" || "/opt/metagenomics/mothur/mothur";

[@info "the folder path of the NCBI blast+ suite for run 
        OTU contig alignment with the greengenes or silva 
        database."]
const blast_bin as string = ?"--blast+" || "/opt/metagenomics/ncbi-blast-2.12.0+/bin";

#region "config workspace"
options(workflow.debug = workflow_debug);
options(mothur = mothur);
options(ncbi_blast = `${blast_bin}/blastn`);
#end region