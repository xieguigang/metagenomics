require(GCModeller);

imports "taxonomy_kit" from "metagenomics_kit";

file = "./16s_results.summary";
tree = taxonomy_kit::read.mothurTree(file);
OTU = tree 
|> as.OTU_table() 
|> as.data.frame()
;

print(OTU, max.print = 10);

write.csv(OTU, file = "./mothur_OTU_table.csv");
