require(GCModeller);

imports "bioseq.fasta" from "seqtoolkit";

let seed = read_stockholm("C:\Users\Administrator\Downloads\Pfam-A.seed");

write.fasta(seed, file = "Z:/Pfam-A.fas", lineBreak=60);