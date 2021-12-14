require(GCModeller);

imports "taxonomy_kit" from "metagenomics_kit";

file = "F:\\16s_test\\16s.trim.contigs.good.good.gg.knn.tax.summary";
tree = taxonomy_kit::read.mothurTree(file);
OTU = tree 
|> as.OTU_table() 
|> as.data.frame()
;

print(OTU, max.print = 10);

write.csv(OTU, file = "F:\\16s_test\\OTU.csv");
