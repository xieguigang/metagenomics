imports "gast" from "metagenomics_kit";

const align_silva as function(blastn = getOption("ncbi_blast")) {

}

#' Align greengenes database for taxonomy annotation
#' 
#' @param work16s the workspace object which is comes 
#'      from the ``greengenes_OTUTaxonomy`` function.
#' @param blastn the executable file path of the NCBI 
#'      blastn program.
#' @param is_debug run in debug mode?
#' 
const align_greengenes as function(work16s, 
                                   blastn   = getOption("ncbi_blast"), 
                                   is_debug = getOption("workflow.debug")) {

    # 进行blastn序列比对操作来完成物种鉴定
    is_debug = as.logical(is_debug);
    blastn_cli = `
        ${blastn} 
            -query  "${work16s$OTU_rep}" 
            -db     "${work16s$greengenes[[1]]}" 
            -out    "${work16s$blastnOut}" 
            -evalue "${work16s$evalue}"
            -num_threads ${work16s$num_threads}
    `;
    
    print("commandline for run taxonomy annotation:");
    print(blastn_cli);    

    if (!is_debug) {
        blastn_cli |> system(intern = TRUE);    
    } else {
        print("skip of run blastn in debug model!");
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

        if (file.exists(csv)) {
            csv 
            |> read.csv() 
            |> print(max.print = 10)
            ;
        }

        # run taxonomy annotation and create OTU 
        # relative abundance table result.
        gast = work16s |> taxonomy_annotation();

        write.csv(gast, file = `${outputdir}/gast_OTUs.csv`);
    } else {
        stop("error during of run taxonomy annotation process.");
    }
}

#' Do taxonomy annotation for the OTU contigs
#' 
#' @param OTU_mapping the blastn search result table.
#' 
#' @return the 16s raw data processing pipeline final result.
#' 
const taxonomy_annotation as function(work16s) {
    print("Do taxonomy annotation for the OTU contigs");
    print("via gast algorithm module!");

    cat("\n\n\n");
    cat("######################################################################################\n");
    cat("# \n");
    cat("# gast: compares trimmed sequences against a reference database For assigning taxonomy\n");
    cat("# \n");
    cat("# Author: Susan Huse, shuse@mbl.edu\n");
    cat("# \n");
    cat("# Date: Fri Feb 25 11:41:10 EST 2011\n");
    cat("# \n");
    cat("# Copyright (C) 2011 Marine Biological Laborotory, Woods Hole, MA\n");
    cat("# \n");
    cat("# This program Is free software; you can redistribute it And/Or\n");
    cat("# modify it under the terms Of the GNU General Public License\n");
    cat("# As published by the Free Software Foundation; either version 2\n");
    cat("# Of the License, Or (at your Option) any later version.\n");
    cat("# \n");
    cat("# This program Is distributed In the hope that it will be useful,\n");
    cat("# but WITHOUT ANY WARRANTY; without even the implied warranty Of\n");
    cat("# MERCHANTABILITY Or FITNESS For A PARTICULAR PURPOSE.  See the\n");
    cat("# GNU General Public License For more details.\n");
    cat("# \n");
    cat("# For a copy Of the GNU General Public License, write To the Free Software\n");
    cat("# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n");
    cat("# Or visit http://www.gnu.org/copyleft/gpl.html\n");
    cat("# \n");
    cat("# Keywords: gast taxonomy refssu\n");
    cat("# \n");
    cat("# Assumptions:\n");
    cat("# \n");
    cat("# Revisions:\n");
    cat("# \n");
    cat("# Programming Notes:\n");
    cat("# \n");
    cat("#####################################################################################\n");
    cat("\n\n\n");

    taxonomy = parse.greengenes_tax(work16s$greengenes$taxonomy);  
    OTU_rep = work16s$OTU_rep |> parse.mothur_OTUs();
    gast = work16s$blastnOut 
        |> read.blast(type = "nucl") 
        |> OTU.taxonomy(
            OTUs           = OTU_rep, 
            taxonomy       = taxonomy, 
            gast.consensus = FALSE
        )
        ;

    gast;
}