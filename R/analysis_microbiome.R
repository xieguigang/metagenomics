let analysis_microbiome = function(otu_table, sampleinfo, output_dir) {
    let bootstrap = file.path(@datadir, "microbiome","bootstrap.R");
    let interop = function(otu_table, sampleinfo, output_dir, pkg_dir) {
        bootstrap_r(pkg_dir);
        process_analysis(data = otu_table, sample_info = sampleinfo, output_dir = output_dir);
    }

    native_r(interop, list(
        otu_table = otu_table, 
        sampleinfo = sampleinfo, 
        output_dir = normalizePath(output_dir), 
        pkg_dir = dirname(bootstrap)
    ));
}