require(Metagenomics);

[@info "A source data directory which contains *.fq raw data files."]
[@type "directory"]
const src as string = ?"--src" || stop("a source data directory path must be provided!");

[@info "number of threads that will be used for running parallel task."]
const num_threads as integer = ?"--num_threads" || 32;

[@info "the reference OTU sequnece database for run taxonomy
        annotation of the OTU contigs which is generated from
        the mothur software. this database file can be download
        from the mothur release page:
        https://mothur.org/wiki/greengenes-formatted_databases/"]
const greengenes as string = ?"--greengenes" || "/opt/metagenomics/greengenes/taxonomy/gg_13_8_99.fasta+gg_13_8_99.gg.tax";
const refalign = greengenes_opts(greengenes);

if (!dir.exists(src)) {
    stop(`invalid data source folder path [${src}]!`);
} else {
    setwd(src);
}

print("get greengenes reference database:");
str(refalign);

# generate file: ./16s.files
Metagenomics::mothur_files(getwd(), "16s.files");

make.contigs(file="./16s.files", processors=num_threads);
screen.seqs(
    fasta="16s.trim.contigs.fasta", 
    minlength=200, 
    maxlength=450, 
    maxambig=0, 
    group="16s.contigs.groups", 
    processors=num_threads
);
align.seqs(
    candidate="16s.trim.contigs.good.fasta", 
    template=refalign$greengenes
);
screen.seqs(
    fasta="16s.trim.contigs.good.align", 
    alignreport="16s.trim.contigs.good.align.report", 
    minsim=90, 
    minscore=10, 
    group="16s.contigs.good.groups"
);
classify.seqs(
    fasta="16s.trim.contigs.good.good.align", 
    template=refalign$greengenes, 
    taxonomy=refalign$taxonomy, 
    method="knn", 
    processors=num_threads,
    numwanted=3
);
summary.tax(
    taxonomy="16s.trim.contigs.good.good.16s_taxonomy_1_3.knn.taxonomy", 
    group="16s.contigs.good.good.groups"
);
split.groups(
    fasta="16s.trim.contigs.good.fasta", 
    group="16s.contigs.good.groups"
);

# result file
# 16s.trim.contigs.good.good.nematode_taxonomy_1_3.knn.tax.summary 
#  -> 16s_results.summary