library(cowplot)
setwd("/Users/Diego/Desktop/NJstM/results")
files=list.files(path=".",pattern="_matrix.csv$")

colors=c("0"="white","1"="black")
plots=list()
for (file in files) {
	name=sub(x=file,pattern=".csv",replacement="")
	title_text=sub(x=name,pattern="r0_(.*)_matrix",replacement="Random \\1%")
	title_text=sub(x=title_text,pattern="r1_(.*)_matrix",replacement="By individual \\1%")
	datai=read.csv(file)
	#plots[[file]]=ggplot(data=datai,aes(x=y,y=x,fill=as.factor(press)))+geom_tile()+scale_fill_manual(values=colors,name="",labels=c("Absent","Present"))+scale_x_discrete(name="Taxon")+scale_y_discrete(name="Gene")+labs(title=title_text)

	plots[[file]]=ggplot(data=datai,aes(x=x,y=y,fill=as.factor(press)))+geom_tile()+scale_fill_manual(guide=FALSE,values=colors,name="",labels=c("Absent","Present"))+scale_x_discrete(name="Gene")+scale_y_discrete(name="Taxon")+labs(title=title_text)
	#save_plot(plot,filename=paste0(name,".pdf"),base_height = 8,base_aspect_ratio = 1.2)
}
final=plot_grid(plotlist=(plots[-1])[c(1,5,2,6,3,7,4,8)],nrow = 4)
save_plot(final,filename="grid_matrices.pdf",base_height=12)
save_plot(final,filename="grid_matrices.png",base_height=12)
