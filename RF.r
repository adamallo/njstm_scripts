library(phangorn)
args <- commandArgs(TRUE)
if (length(args)!=1 || !file.info(args[1])$isdir){
        print("Usage Rscript RF.r directory")
        quit()
}

setwd(args[1])
dirs=list.files(path=".",pattern="^[0-9][0-9][0-9][0-9][0-9]$",recursive=FALSE)

final_data=data.frame(rep=numeric(0),rf=numeric(0),mdata=character(),method=character())

for (i in 1:length(dirs)){
	dir=dirs[i]
	print(paste0("Working on ",dir))
	rep=as.numeric(basename(dir))
	s_tree=read.tree(paste0(dir,"/s_tree.trees"))
	labels_original=sort(s_tree$tip.label)
	methods=c("lnjst","onjst","unjst","astral","astridmo","astridmu","astriddef")
	options(stringsAsFactors=FALSE)
	for (m in 1:length(methods))
	{
		method=methods[m]
		trees=list.files(path=dir,pattern=paste0("*.",method,"$"))
		for (j in 1:length(trees))
		{
			file=paste(dir,trees[j],sep="/")
			condition=gsub(paste0(".",method),"",basename(file))
			dist=tryCatch({
				tree=read.tree(file)
				labels2=sort(tree$tip.label)
				todel=setdiff(labels_original,labels2) ##The gene tree may be smaller if all individuals of one species have been dropped from all replicates
				if (length(todel)!=0){
					print(paste0("The file ",file," does not have the same set of labels as the species tree. Different labels: ",paste(todel,collapse=" ")))
					temp_stree=drop.tip(s_tree,tip=todel) #Prunes tips and internal branches
					RF.dist(temp_stree,tree,check.labels=TRUE,normalize=TRUE)
				} else {
					RF.dist(s_tree,tree,check.labels=TRUE,normalize=TRUE)
				}
			},warning=function(war){
				print(paste0("The file ",file," generated a NA"))
				return(NA)
			},error=function(err){
				print(paste0("The file ",file," generated a NA"))
				return(NA)
			})
			#dist=RF.dist(s_tree,tree,check.labels = TRUE)/((s_tree$Nnode+1)*2-6) ##s_tree$Nnode= Number of internal nodes. This tree is rooted, so internal nodes+1 = n_leaves. 2*(n-3) = number of internal branches/bipartitions in an unrooted tree * 2.
			final_data=rbind(final_data,data.frame(rep=rep,rf=dist,mdata=condition,method=method))
		}
	}
}
final_data$rep=as.numeric(final_data$rep)
final_data$rf=as.numeric(final_data$rf)
write.table(final_data,file="rf.csv",sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)
