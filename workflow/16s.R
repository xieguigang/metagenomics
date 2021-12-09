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
const right as string = ?"--right" || stop("Missing of the right pair-end short reads file!");

[@info "the directory folder path for output data result files."]
[@type "directory"]
const outputdir as string = ?"--outputdir" || `${dirname(left)}/16s_result/`;

[@info "an optional argument for setting up the processors that will be 
        used in mothur program, by default is using 72 processor 
        threads."]
const num_threads as integer = ?"--num_threads" || 72;

[@info "run debug of the workflow?"]
const workflow_debug as boolean = ?"--workflow-debug";

[@info "file path to the mothur program executable file."]
const mothur as string = ?"--mothur" || "/opt/metagenomics/mothur/mothur";

[@info "the folder path of the NCBI blast+ suite for run 
        OTU contig alignment with the greengenes or silva 
        database."]
const blast_bin as string = ?"--blast+" || "/opt/metagenomics/ncbi-blast-2.12.0+/bin";

[@info "skip of run mothur workflow for generate OTU sequence dataset? 
        this options is usually apply for run workflow debug."]
const skip_mothur as boolean = ?"--skip-mothur";

[@info "the reference OTU sequnece database for run taxonomy 
        annotation of the OTU contigs which is generated from 
        the mothur software. this database file can be download
        from the mothur release page: 
        https://mothur.org/wiki/greengenes-formatted_databases/"]
const greengenes as string = ?"--greengenes" || "/opt/metagenomics/greengenes/taxonomy/gg_13_8_99.fasta+gg_13_8_99.gg.tax";
[@info "the reference template file for run mothur alignment of 
        the generated contig fasta sequence file."]
[@type "*.fasta"]
const template as string = ?"--template" || "/opt/metagenomics/greengenes/refalign/gg_13_8_99.refalign";

#region "config workspace"
options(workflow.debug = workflow_debug);
options(mothur = mothur);
options(ncbi_blast = `${blast_bin}/blastn`);
options(greengenes = greengenes);
options(mothur_template = template);
#end region

Metagenomics::greengenes_OTUTaxonomy(
        left, right, 
        outputdir   = outputdir, 
        num_threads = num_threads,
        skip_mothur = skip_mothur
);