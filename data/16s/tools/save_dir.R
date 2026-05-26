save_dir = function(output_dir, name) {
    output_dir = file.path(output_dir, name);
    dir.create(output_dir, showWarnings = FALSE);
    return(output_dir);
}