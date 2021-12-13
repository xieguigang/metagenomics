
const align_silva as function() {

}

const align_greengenes as function(work16s) {
    # 进行blastn序列比对操作来完成物种鉴定
    blastn_cli = `
        ${blastn} 
            -query  "${work16s$OTU_rep}" 
            -db     "${greengenes[[1]]}" 
            -out    "${work16s$blastnOut}" 
            -evalue "${evalue}"
            -num_threads ${num_threads}
    `;
    
    print("commandline for run taxonomy annotation:");
    print(blastn_cli);    

    if (!is_debug) {
        blastn_cli |> system(intern = TRUE);    
    }

    if (file.exists(work16s$blastnOut)) {
        # run blastn parser
        csv = `${dirname(work16s$blastnOut)}/${basename(work16s$blastnOut)}.csv`;

        using buffer as open.stream(csv, type = "Mapping", ioRead = FALSE) {
            work16s$blastnOut 
            |> getBlastnMapping()
            |> stream.flush(stream = buffer)
            ;
        }

        # run taxonomy annotation and create OTU relative abundance table result.
        
    } else {
        stop("error during of run taxonomy annotation process.");
    }
}