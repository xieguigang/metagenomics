
#' Build OTU based on mothur
#' 
const mothur_OTU as function(left, right, 
                             outputdir   = "./", 
                             num_threads = 32) {

    const work16s = list(
        left        = left, 
        right       = right, 
        outputdir   = outputdir, 
        num_threads = num_threads
    );

    make.contigs(left, right, outputdir, num_threads = num_threads);
    work16s$contig.fasta = write_contig("16s.trim.contigs.fasta");

    # RunAutoScreen
    # summary.seqs + screen.seqs 
    # 
    screen.seqs(work16s$contig.fasta, num_threads = num_threads);

    work16s$contigs = "contig.good.fasta";
    work16s$groups  = "16s.contigs.good.groups";

    unique.seqs(work16s$contigs, logfile = "[4]unique.seqs.txt");

    work16s$names = "contig.good.names"; 
    count.seqs(work16s$names, groups= work16s$groups);

    work16s$count_table = "contig.good.count_table";
    summary.seqs(
        seqfile="contig.good.unique.fasta", 
        count_table=work16s$count_table,
        num_threads = num_threads, 
        logfile = "[6]summary.seqs.txt"
    );
    # contig.good.unique.summary

    write_contig("contig.good.unique.fasta");

    work16s$contigs = "contig.fasta";

    
    align.seqs(fasta=$contigs,reference=$silva,flip=T,processors=$num_threads)","[7]align.seqs.txt");
    # contig.align
    # contig.align.report
    # contig.flip.accnos

    $align = "contig.align";
    filter.seqs(fasta=$align,processors=$num_threads)", "[8]filter.seqs.txt");
    # contig.filter
    # contig.filter.fasta

    write_contig("contig.filter.fasta");
    unique.seqs(fasta=$contigs)", "[9]unique.seqs.txt");
    # contig.names
    # contig.unique.fasta

    
$align = "contig.unique.fasta";
runMothur("dist.seqs(fasta=$align,calc=onegap,countends=F,cutoff=0.03,output=lt,processors=$num_threads)", "[10]dist.seqs.txt");
# contig.unique.phylip.dist

my $dist = "contig.unique.phylip.dist";
runMothur("cluster(phylip=$dist,method=furthest,cutoff=0.03,processors=$num_threads)", "[11]cluster.txt");

# contig.unique.phylip.fn.sabund
# contig.unique.phylip.fn.rabund
# contig.unique.phylip.fn.list

my $list = "contig.unique.phylip.fn.list";
runMothur("bin.seqs(list=$list,fasta=$contigs,name=contig.names)", "[12]bin.seqs.txt");

# contig.unique.phylip.fn.unique.fasta
# contig.unique.phylip.fn.0.01.fasta
# contig.unique.phylip.fn.0.02.fasta
# contig.unique.phylip.fn.0.03.fasta
runMothur("get.oturep(phylip=$dist,fasta=contig.unique.fasta,list=$list,label=0.03)", "[13]get.oturep.txt");

print "Mothur job done!\n";

# 在这里进行SILVA的16S数据库的比对操作，进行OTU序列所属的物种鉴定
# 首先需要将OTU的fasta文件之中由于前面的mothur程序align的空格和连接符都删除掉
# 否则blastn程序会报错
my $fasta   = "./contig.unique.phylip.fn.0.03.rep.fasta";
my $data    = read_file $fasta, {binmode => ':utf8'};
my $OTU_rep = "./OTU.rep.fasta"; 

$data =~ s/[.-]//g;
write_file $OTU_rep, {binmode => ':utf8'}, $data;

# 进行blastn序列比对操作来完成物种鉴定
$CLI = "$blastn -query $OTU_rep -db $greengene -out ./OTU_greengene_99.txt -evalue 1e-50 -num_threads $num_threads";
print $CLI."\n";
system($CLI);

print("greengenes OTU Taxonomy align job done!");
}