const pairEnds as function(dir) {
    dir
    |> list.files(pattern = "*.fq")
    |> groupBy(function(path) {
        basename(path) 
        |> gsub("_R1", "") 
        |> gsub("_R2", "")
    })
    |> lapply(function(pairend) {
        print(pairend);
        
        list(
            R1 = pairend[pairend like $".+_R1.fq"],
            R2 = pairend[pairend like $".+_R2.fq"]
        );
    }, names = p -> p$key);
}