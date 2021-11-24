
#' Build OTU based on mothur
#' 
const mothur_OTU as function(left, right, 
                             outputdir   = "./", 
                             num_threads = 32) {

    make.contigs(left, right, outputdir, num_threads = num_threads);
    const contig.fasta = write_contig("16s.trim.contigs.fasta");
    # RunAutoScreen
    screen.seqs(contig.fasta, num_threads = num_threads);
}