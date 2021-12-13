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

[@info "the e-value cutoff of the blastn search pipeline. default e-value is 1e-10."]
const evalue as double = ?"--evalue" || 1e-10;

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

#region "make.contigs assembler"
[@info "The align parameter allows you to specify the alignment method to
        use. Your options are: ``gotoh`` and ``needleman``. The default
        is ``needleman``."]
const alignment as string = ?"--alignment" || "needleman";

[@info "These parameters are used while aligning the sequence read to
        determine the overlap. The match parameter allows you to specify
        the bonus for having the same base. The default is 1.0. The
        mistmatch parameter allows you to specify the penalty for having
        different bases. The default is -1.0. The gapopen parameter allows
        you to specify the penalty for opening a gap in an alignment.
        The default is -2.0. The gapextend parameter allows you to specify
        the penalty for extending a gap in an alignment. The default is
        -1.0."]
[@type "match, mismatch, gapopen, gapextend"]
const scores as string = ?"--scores" || "1.0,-1.0,-2.0,-1.0";

[@info "The insert parameter allows you to set a quality scores threshold.
        When we are merging the overlapping regions, in the case where we
        are trying to decide whether to keep a base or remove it because
        the base is compared to a gap in the other fragment, if the base has
        a quality score below the threshold we eliminate it. Default=20."]
const insert as double = ?"--insert" || 20;
#end region

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
    skip_mothur = skip_mothur,
    evalue      = evalue,
    make.contigs = list(
        algorithm = alignment,
        score     = scores
        |> strsplit("\s*,\s*", fixed = FALSE, perl = TRUE)
        |> as.numeric(),
        insert    = insert
    )
);
