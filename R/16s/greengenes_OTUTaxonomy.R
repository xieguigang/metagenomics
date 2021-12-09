
#' Workflow for 16s
#' 
#' @description Workflow for create 16s data analysis based on the mothur
#'    software and NCBI blastn program. Workflow steps in this function 
#'    includes:
#'  
#'    1. create OTU sequence based on the mothur workflow
#'    2. alignment of the OTU sequence with the greengenes database
#'       via NCBI blastn pipeline
#'
const greengenes_OTUTaxonomy as function(left, right, 
                                         outputdir   = "./", 
                                         num_threads = 32) {

    const is_debug as boolean  = getOption("workflow.debug");
    const blastn as string = getOption("ncbi_blast");
    const refalign as string = getOption("mothur_template");
    const [greengenes, taxonomy] = Metagenomics::greengenes_opts();

    # config to absolute file path.
    left      = normalizePath(left);
    right     = normalizePath(right);
    outputdir = normalizePath(outputdir);

    # 使用mothur程序组装测序结果为contig，生成OTU序列文件
    Metagenomics::mothur_OTU(left, right, refalign, outputdir, num_threads); 

    # 在这里进行SILVA的16S数据库的比对操作，进行OTU序列所属的物种鉴定
    # 首先需要将OTU的fasta文件之中由于前面的mothur程序align的空格和连接符都删除掉
    # 否则blastn程序会报错
    work16s$fasta   = `${outputdir}/contig.unique.phylip.fn.0.03.rep.fasta`;
    work16s$OTU_rep = `${outputdir}/OTU.rep.fasta`; 

    work16s$fasta
    |> readLines()
    |> writeLines(con = work16s$OTU_rep)
    ;

    # 进行blastn序列比对操作来完成物种鉴定
    const blastn_cli as string = `
        ${blastn} 
            -query  ${work16s$OTU_rep} 
            -db     ${greengenes} 
            -out    "${outputdir}/OTU_greengene_99.txt" 
            -evalue 1e-50 
            -num_threads ${num_threads}
    `;
    
    print("commandline for run taxonomy annotation:");
    print(blastn_cli);    

    if (!is_debug) {
        blastn_cli |> system(intern = TRUE);    
    }

    print("greengenes OTU Taxonomy align job done!");
}

#' parse greengenes option value
#' 
#' @return a list object that contains two elements:
#'    1. fasta: is the file path of the fasta sequence file
#'    2. taxonomy: is the file path of the taxonomy data for 
#'                 the data annotation of the OTU sequence
#' 
const greengenes_opts as function() {
    const greengenes as string = getOption("greengenes");
    const name as string       = basename(greengenes);
    const parts as string      = strsplit(name, "+", fixed = TRUE);
    const repository as string = dirname(greengenes);

    list(
        greengenes = `${repository}/${parts[1]}`,
        taxonomy   = `${repository}/${parts[2]}`
    );
}