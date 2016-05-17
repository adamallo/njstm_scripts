library(phangorn)
args <- commandArgs(TRUE)
if (length(args)!=1 || !file.info(args[1])$isdir){
        print("Usage Rscript RF.r directory")
        quit()
}

dir=args[1]
setwd(dir)
rep=as.numeric(basename(dir))
s_tree=read.tree("s_tree.trees")
methods=c("lnjst","onjst","unjst")
options(stringsAsFactors=FALSE)
final_data=data.frame(rep=numeric(0),rf=numeric(0),mdata=character(),method=character())
for (m in 1:length(methods))
{
	method=methods[m]
	trees=list.files(path=".",pattern=paste0("*.",method))
	for (i in 1:length(trees))
	{
		file=trees[i]
		condition=gsub(paste0(".",method),"",file)
		tree=read.tree(file)
		dist=RF.dist(s_tree,tree,check.labels = TRUE)/((s_tree$Nnode+1)*2-6) ##s_tree$Nnode= Number of internal nodes. This tree is rooted, so internal nodes+1 = n_leaves. 2*(n-3) = number of internal branches/bipartitions in an unrooted tree * 2.
		final_data=rbind(final_data,data.frame(rep=rep,rf=dist,mdata=condition,method=method))
	}
}
final_data$rep=as.numeric(final_data$rep)
final_data$rf=as.numeric(final_data$rf)
header=FALSE
if (rep==1)
{
	header=TRUE
}
write.table(final_data,file="rf.csv",sep=",",quote=FALSE,row.names=FALSE,col.names=header)
