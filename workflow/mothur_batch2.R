require(Metagenomics);
require(GCModeller);

imports "taxonomy_kit" from "metagenomics_kit";

[@info "A source data directory which contains *.fq raw data files."]
[@type "directory"]
const src as string = ?"--src" || stop("a source data directory path must be provided!");

[@info "number of threads that will be used for running parallel task."]
const num_threads as integer = ?"--num_threads" || 32;

[@info "file path to the mothur program executable file."]
const mothur as string = ?"--mothur" || "/opt/metagenomics/mothur/mothur";

[@info "the data directory path for save the result data files."]
const outputdir as string = ?"--outputdir" || file.path(src, "16s_results/");

[@info "disable of the file cache?"]
const disable_cache as boolean = ?"--cache-disable";

[@info "the reference template file for run mothur alignment of
        the generated contig fasta sequence file."]
[@type "*.fasta"]
const template as string = ?"--template" || "/opt/metagenomics/greengenes/refalign/gg_13_8_99.refalign";

[@info "the reference OTU sequnece database for run taxonomy
        annotation of the OTU contigs which is generated from
        the mothur software. this database file can be download
        from the mothur release page:
        https://mothur.org/wiki/greengenes-formatted_databases/"]
const greengenes as string = ?"--greengenes" || "/opt/metagenomics/greengenes/taxonomy/gg_13_8_99.fasta+gg_13_8_99.gg.tax";
const refalign = greengenes_opts(greengenes);

dir.create(outputdir, recursive = TRUE);

refalign$template = template;

const work16s = list(
    outputdir = normalizePath(outputdir),
    refalign = refalign
);

if (!dir.exists(src)) {
    stop(`invalid data source folder path [${src}]!`);
} else {    
    setwd(src);
}

options(mothur = mothur);

print("get greengenes reference database:");
str(refalign);
print("data result files will be saved at:");
print(outputdir);

refalign |> workflow2(
    outputdir     = outputdir, 
    num_threads   = num_threads, 
    disable_cache = disable_cache
);