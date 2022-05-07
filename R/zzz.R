require(GCModeller);

imports "taxonomy_kit" from "metagenomics_kit";
imports "microbiome" from "metagenomics_kit";

#' loading the required packages
#' 
const .onLoad = function() {
    print("package for processing of the Metagenomics/16s microbiome data");
    print("base on the GCModeller platform");
}
