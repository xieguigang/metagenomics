require(GCModeller);

imports "microbiome" from "metagenomics_kit";

const otu_file as string = ?"--out_table" || stop("A mothur OTU table file must be provided!");
const outfile  as string = ?"--out_KO"    || `${dirname(otu_file)}/metagenome.csv`;
const matrixKO as string = ?"--matrix"    || "/opt/metagenomics/ko_13_5_precalculated.PICRUSt";

using PICRUSt as file(matrixKO) {

	OTUtable = read.csv(otu_file, row.names = "taxonomy");
	OTUtable[, "OTU_num"] = NULL;
	
	metagenome = PICRUSt
	|> read.PICRUSt_matrix()
	|> predict_metagenomes(OTUtable)
	|> as.data.frame()
	;
	
	# run enrichment analysis
	write.csv(metagenome, file = outfile, row.names = TRUE);
}