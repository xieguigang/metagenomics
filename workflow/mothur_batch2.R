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
print(work16s$outputdir);

# generate file: ./16s.files
Metagenomics::mothur_files(getwd(), file.path(outputdir, "16s.files"));

setwd(work16s$outputdir);

const check_filecache as function(filename) {
    if (disable_cache) {
        FALSE;
    } else {
        file.exists(filename);
    }
}

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

print("parse mothur result file:");
print(normalizePath(file));
print("get the final OTU table result:");
str(OTU);

write.csv(OTU, file = "./mothur_OTU_table.csv", row.names = FALSE);