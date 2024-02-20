// export R# source type define for javascript/typescript language
//
// package_source=Metagenomics

declare namespace Metagenomics {
   module _ {
      /**
      */
      function onLoad(): object;
   }
   module align {
      /**
        * @param flip default value Is ``F``.
        * @param num_threads default value Is ``8``.
        * @param logfile default value Is ``[7]align.seqs.txt``.
      */
      function seqs(contigs: any, reference: any, flip?: any, num_threads?: any, logfile?: any): object;
   }
   /**
     * @param blastn default value Is ``Call "getOption"("ncbi_blast")``.
     * @param is_debug default value Is ``Call "getOption"("workflow.debug")``.
   */
   function align_greengenes(work16s: any, blastn?: any, is_debug?: any): object;
   /**
     * @param blastn default value Is ``Call "getOption"("ncbi_blast")``.
   */
   function align_silva(blastn?: any): object;
   module bin {
      /**
        * @param logfile default value Is ``[12]bin.seqs.txt``.
      */
      function seqs(list: any, contigs: any, contig.names: any, logfile?: any): object;
   }
   module classify {
      /**
        * @param method default value Is ``knn``.
        * @param numwanted default value Is ``3``.
        * @param processors default value Is ``2``.
        * @param logfile default value Is ``[14]classify.seqs.txt``.
      */
      function seqs(fasta: any, reference: any, taxonomy: any, method?: any, numwanted?: any, processors?: any, logfile?: any): object;
   }
   /**
     * @param method default value Is ``furthest``.
     * @param cutoff default value Is ``0.03``.
     * @param num_threads default value Is ``32``.
     * @param logfile default value Is ``[11]cluster.txt``.
   */
   function cluster(dist: any, method?: any, cutoff?: any, num_threads?: any, logfile?: any): object;
   module count {
      /**
      */
      function seqs(names: any, groups: any): object;
   }
   module dist {
      /**
        * @param calc default value Is ``onegap``.
        * @param countends default value Is ``F``.
        * @param cutoff default value Is ``0.03``.
        * @param output default value Is ``lt``.
        * @param num_threads default value Is ``32``.
        * @param logfile default value Is ``[10]dist.seqs.txt``.
      */
      function seqs(align_fasta: any, calc?: any, countends?: any, cutoff?: any, output?: any, num_threads?: any, logfile?: any): object;
   }
   module filter {
      /**
        * @param logfile default value Is ``[8]filter.seqs.txt``.
      */
      function seqs(align: any, num_threads: any, logfile?: any): object;
   }
   module get {
      /**
        * @param label default value Is ``0.03``.
        * @param logfile default value Is ``[13]get.oturep.txt``.
      */
      function oturep(dist: any, contig.unique.fasta: any, list: any, label?: any, logfile?: any): object;
   }
   /**
   */
   function getBlastnMapping(outTxt: any): object;
   /**
     * @param greengenes default value Is ``Call "getOption"("greengenes")``.
   */
   function greengenes_opts(greengenes?: any): object;
   /**
     * @param outputdir default value Is ``./``.
     * @param num_threads default value Is ``32``.
     * @param evalue default value Is ``1E-10``.
     * @param make.contigs default value Is ``Call "list"("algorithm" <- "needleman", "score" <- [1, -1, -2, -1], "insert" <- 20)``.
     * @param skip_mothur default value Is ``false``.
   */
   function greengenes_OTUTaxonomy(left: any, right: any, outputdir?: any, num_threads?: any, evalue?: any, make.contigs?: any, skip_mothur?: any): object;
   module make {
      /**
        * @param algorithm default value Is ``["needleman", "gotoh"]``.
        * @param score default value Is ``[1, -1, -2, -1]``.
        * @param insert default value Is ``20``.
        * @param num_threads default value Is ``8``.
      */
      function contigs(left: any, right: any, algorithm?: any, score?: any, insert?: any, num_threads?: any): object;
   }
   /**
     * @param filename default value Is ``./16s.files``.
   */
   function mothur_files(dir: any, filename?: any): object;
   /**
     * @param outputdir default value Is ``./``.
     * @param assembler default value Is ``Call "list"("algorithm" <- "needleman", "score" <- [1, -1, -2, -1], "insert" <- 20)``.
     * @param num_threads default value Is ``32``.
   */
   function mothur_OTU(left: any, right: any, refalign: any, outputdir?: any, assembler?: any, num_threads?: any): object;
   /**
   */
   function mothur_workflow(work16s: any): object;
   /**
   */
   function mothurArgvs(command: any, argv: any): object;
   /**
   */
   function mothurErrorParser(stdout: any): object;
   /**
   */
   function pairEnds(dir: any): object;
   /**
   */
   function runMothur(command: any, argv: any, log: any): object;
   module screen {
      /**
        * @param num_threads default value Is ``8``.
      */
      function seqs(contigs: any, num_threads?: any): object;
   }
   module split {
      /**
        * @param logfile default value Is ``[16]split.groups.txt``.
      */
      function groups(fasta: any, group: any, logfile?: any): object;
   }
   module summary {
      /**
        * @param num_threads default value Is ``8``.
        * @param logfile default value Is ``log.txt``.
      */
      function seqs(seqfile: any, count_table: any, num_threads?: any, logfile?: any): object;
      /**
        * @param logfile default value Is ``[15]summary.tax.txt``.
      */
      function tax(taxonomy: any, group: any, logfile?: any): object;
   }
   /**
   */
   function taxonomy_annotation(work16s: any): object;
   module unique {
      /**
        * @param logfile default value Is ``[4]unique.seqs.txt``.
      */
      function seqs(contigs: any, logfile?: any): object;
   }
   /**
     * @param outputdir default value Is ``./``.
     * @param num_threads default value Is ``8``.
     * @param disable_cache default value Is ``false``.
   */
   function workflow2(refalign: any, outputdir?: any, num_threads?: any, disable_cache?: any): object;
   /**
     * @param contig.fasta default value Is ``contig.fasta``.
   */
   function write_contig(seqfile: any, contig.fasta?: any): object;
}
