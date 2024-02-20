require(GCModeller);

imports "taxonomy_kit" from "metagenomics_kit";

let tree = Ncbi.taxonomy_tree("E:\\metagenomics\\data");
let species = ranks(tree) |> taxonomy_ranks("species");
let ncbi_taxid = [unlist(species)]::taxid;

str(species);
print(ncbi_taxid);

for(node in tqdm(ncbi_taxid)) {
    # str(node);
    print(tree |> lineage(node));
}

pause();