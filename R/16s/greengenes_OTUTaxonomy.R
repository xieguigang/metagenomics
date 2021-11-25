const greengenes_OTUTaxonomy as function(left, right, 
                                         outputdir   = "./", 
                                         num_threads = 32) {

    const is_debug as boolean   = getOption("workflow.debug");

    mothur_OTU(left, right, outputdir, num_threads); 

    # 在这里进行SILVA的16S数据库的比对操作，进行OTU序列所属的物种鉴定
    # 首先需要将OTU的fasta文件之中由于前面的mothur程序align的空格和连接符都删除掉
    # 否则blastn程序会报错
    work16s$fasta   = "./contig.unique.phylip.fn.0.03.rep.fasta";
    work16s$OTU_rep = "./OTU.rep.fasta"; 

    work16s$fasta
    |> readLines()
    |> writeLines(con = $OTU_rep)
    ;

    # 进行blastn序列比对操作来完成物种鉴定
    const blastn_cli as string = `$blastn -query $OTU_rep -db $greengene -out ./OTU_greengene_99.txt -evalue 1e-50 -num_threads $num_threads`;
    
    print(blastn_cli);    

    if (!is_debug) {
        blastn_cli |> system(intern = TRUE);    
    }

    print("greengenes OTU Taxonomy align job done!");
}