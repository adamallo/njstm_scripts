library(RSQLite)
library(cowplot)
library(pogit)
library(xtable)
library(glmnet)
library(ggrepel)
setwd("/Users/Diego/Desktop/NJstM/results")

###MY FUNCTIONS###
compare2_plot=function(x,methods,names,glm_option=c("glm","poly"),normalized=FALSE){
  if(normalized == TRUE)
  {
    rf="normrfdif"
  } else {
    rf="rfdif"
  }
  conditions=levels(x$mdata)
  n_comb=length(conditions)
  wilcox1=list()
  meandif=list()
  prob_new_atleastasgoodas_old=list()
  prob_new_betterthan_old=list()
  prob_old_atleastasgoodas_new=list()
  prob_old_betterthan_new=list()
  ratio_atleastasgoodas_newold=list()
  ratio_betterthan_newold=list()
  for (cond in conditions)
  {  
    tempdata1=x[x$method==methods[1] & x$mdata==cond,c("rep","rf")]
    colnames(tempdata1)=c("rep","rf1")
    tempdata2=x[x$method==methods[2] & x$mdata==cond,c("rep","rf")]
    colnames(tempdata2)=c("rep","rf2")
    tempdata=merge(tempdata1,tempdata2)
    tempdata=tempdata[complete.cases(tempdata),]
    vector_original=tempdata$rf1
    vector_new=tempdata$rf2
    
    p_value=wilcox.test(vector_original,vector_new,paired=TRUE)$p.value
    add=ifelse(p_value > 0.05/n_comb, "", ifelse(p_value > 0.01/n_comb, " *", 
                                                 ifelse(p_value > 0.001/n_comb, " **", " ***")))
    wilcox1[cond]=paste0(signif(p_value,digits=2),add)
    meandif[cond]=mean(vector_original-vector_new,na.rm=TRUE)
    prob_new_atleastasgoodas_old[cond]=mean(vector_new<=vector_original,na.rm=TRUE)
    prob_new_betterthan_old[cond]=mean(vector_new<vector_original,na.rm=TRUE)
    prob_old_atleastasgoodas_new[cond]=mean(vector_original<=vector_new,na.rm=TRUE)
    prob_old_betterthan_new[cond]=mean(vector_original<vector_new,na.rm=TRUE)
    ratio_atleastasgoodas_newold[cond]=mean(vector_new<=vector_original,na.rm=TRUE)/mean(vector_original<=vector_new,na.rm=TRUE)
    ratio_betterthan_newold[cond]=mean(vector_new<vector_original,na.rm=TRUE)/mean(vector_original<vector_new,na.rm=TRUE)
  }
  
  datafig=x[x$method==methods[1],c(1,2,3,5,6)]
  datafig=datafig[with(datafig,order(rep,mdata)),]
  temp=x[x$method==methods[2],]
  datafig$normrfdif=datafig$rf
  datafig$rfdif=datafig$rf-temp[with(temp,order(rep,mdata)),]$rf
  
  ##Duplicating the normal dataset since it is common for the two missing modes
  temp=datafig[datafig$pmissing==0,]
  temp$missingmode="taxa"
  datafig=rbind(datafig,temp)
  
  ##Some things are hardcoded in this section
  textplot1=data.frame(mdata=vector(mode = "character",length = length(prob_new_atleastasgoodas_old)+1),pmissing=vector(mode="double",length = length(prob_new_atleastasgoodas_old)+1),missingmode=vector(mode="character",length = length(prob_new_atleastasgoodas_old)+1),text=vector(mode = "character",length = length(prob_new_atleastasgoodas_old)+1),meany=vector(mode = "double",length=length(prob_new_atleastasgoodas_old)+1),ycoord=vector(mode = "double",length=length(prob_new_atleastasgoodas_old)+1),stringsAsFactors = FALSE)
  outi=1
  typetext=""
  #while (i<=length(prob_new_atleastasgoodas_old))
  for (i in 1:length(prob_new_atleastasgoodas_old))
  {
    st=gsub(x=wilcox1[i],"[^\\*]*(\\*+)[^\\*]*","\\1")
    if(!is.na(as.numeric(st)))
    {
      st="" ##no significant
    }
    
    name=names(prob_new_atleastasgoodas_old[i])[1]
    meany=mean(datafig[datafig$mdata==name,]$rfdif,na.rm=TRUE)
    ycoord=meany
    v=as.numeric(sub("r._(.*)","\\1",x=name))*0.01
    if (is.na(v))
    {
      ycoord=meany##+0.05
      textplot1[outi,]=c(name,0.0,"taxa",paste0("I:",paste(paste(round(as.numeric(prob_new_atleastasgoodas_old[i]),digits=2),round(as.numeric(ratio_betterthan_newold[i]),digits=2),sep=" | "),st)),meany,ycoord) #Duplicating for original, since there is no missing data(or both with r=0)
      v=0.0
      datafig[datafig$mdata==name & datafig$missingmode=="taxa",]$normrfdif=datafig[datafig$mdata==name & datafig$missingmode=="taxa",]$rfdif/mean(datafig[datafig$mdata==name & datafig$missingmode=="taxa",]$normrfdif)
      outi=i+1
    }
    pmissing=v
    if(sub("r(.)_.*","\\1",x=name)=="1")
    {
      v="taxa"
      ycoord=meany##+0.05
      typetext="I:"
    } else {
      v="random"
      ycoord=meany##-0.05 I am using ggrepel now
      typetext="R:"
    }
    missingmode=v
    textplot1[outi,]=c(name,pmissing,missingmode,paste(typetext,paste(paste(round(as.numeric(prob_new_atleastasgoodas_old[i]),digits=2),round(as.numeric(ratio_betterthan_newold[i]),digits=2),sep=" | "),st)),meany,ycoord)
    outi=outi+1
    datafig[datafig$mdata==name & datafig$missingmode==missingmode,]$normrfdif=datafig[datafig$mdata==name & datafig$missingmode==missingmode,]$rfdif/mean(datafig[datafig$mdata==name & datafig$missingmode==missingmode,]$normrfdif)
  }
  textplot1$mdata=as.factor(textplot1$mdata)
  textplot1$missingmode=as.factor(textplot1$missingmode)
  textplot1$meany=as.numeric(textplot1$meany)
  textplot1$ycoord=as.numeric(textplot1$ycoord)
  textplot1$pmissing=as.numeric(textplot1$pmissing)
  if(glm_option=="poly"){
    plot1=ggplot(data=datafig,aes_string(x="pmissing",y=rf,group="as.factor(interaction(pmissing,missingmode))",color="missingmode"))+geom_violin(bw=0.005,aes(fill=missingmode))+stat_summary(aes(shape=missingmode),color="black",size=0.8)+geom_smooth(method="glm",formula=y~poly(x,2),aes(group=NULL))+geom_text_repel(data = textplot1,aes(x=pmissing,y=ycoord,label=text,group=NULL,color=NULL))+scale_y_continuous(name=paste0("RF difference between ",names[1]," and ",names[2]))+scale_x_continuous("Proportion of missing taxa")+scale_color_brewer(name="Missing data mode",type="qual",palette = 2,labels=c("Random","By individual"),guide=FALSE)+scale_fill_brewer(name="Missing data mode",type="qual",palette = 2,labels=c("Random","By individual"))+scale_shape_discrete(name="Missing data mode",labels=c("Random","By individual"))
  }else {
    plot1=ggplot(data=datafig,aes_string(x="pmissing",y=rf,group="as.factor(interaction(pmissing,missingmode))",color="missingmode"))+geom_violin(bw=0.005,aes(fill=missingmode))+stat_summary(aes(shape=missingmode),color="black",size=0.8)+geom_smooth(method="glm",aes(group=NULL))+geom_text_repel(data = textplot1,aes(x=pmissing,y=ycoord,label=text,group=NULL,color=NULL))+scale_y_continuous(name=paste0("RF difference between ",names[1]," and ",names[2]))+scale_x_continuous("Proportion of missing taxa")+scale_color_brewer(name="Missing data mode",type="qual",palette = 2,labels=c("Random","By individual"),guide=FALSE)+scale_fill_brewer(name="Missing data mode",type="qual",palette = 2,labels=c("Random","By individual"))+scale_shape_discrete(name="Missing data mode",labels=c("Random","By individual"))
  }
   
  return(list(plot1,datafig,textplot1))
}

RsqGLMdm=function(obs = NULL, pred = NULL, model = NULL) { ##Not really mine, bugfix of the modEvA function
  model.provided <- ifelse(is.null(model), FALSE, TRUE)
  
  if (model.provided) {
    if (!("glm" %in% class(model))) stop ("'model' must be of class 'glm'.")
    if (!is.null(pred)) message("Argument 'pred' ignored in favour of 'model'.")
    if (!is.null(obs)) message("Argument 'obs' ignored in favour of 'model'.")
    obs <- model$y
    pred <- model$fitted.values
    
  } else { # if model not provided
    if (is.null(obs) | is.null(pred)) stop ("You must provide either 'obs' and 'pred', or a 'model' object of class 'glm'")
    if (length(obs) != length(pred)) stop ("'obs' and 'pred' must be of the same length (and in the same order).")
    ##DM bugfix in this stop
    if (any(!(obs %in% c(0, 1))) | any(pred < 0) | any(pred > 1)) stop ("Sorry, 'obs' and 'pred' options currently only implemented for binomial GLMs (binary response variable with values 0 or 1) with logit link.")
    logit <- log(pred / (1 - pred))
    model <- glm(obs ~ logit, family = "binomial")
  }
  
  null.mod <- glm(obs ~ 1, family = family(model))
  loglike.M <- as.numeric(logLik(model))
  loglike.0 <- as.numeric(logLik(null.mod))
  N <- length(obs)
  
  # based on Nagelkerke 1991:
  CoxSnell <- 1 - exp(-(2 / N) * (loglike.M - loglike.0))
  Nagelkerke <- CoxSnell / (1 - exp((2 * N ^ (-1)) * loglike.0))
  
  # based on Allison 2014:
  McFadden <- 1 - (loglike.M / loglike.0)
  Tjur <- mean(pred[obs == 1]) - mean(pred[obs == 0])
  sqPearson <- cor(obs, pred) ^ 2
  
  return(list(CoxSnell = CoxSnell, Nagelkerke = Nagelkerke, McFadden = McFadden, Tjur = Tjur, sqPearson = sqPearson))
}

mcmctoglm=function(mcmc,min_posterior=0.8,names=NULL) {
  table=summary(mcmc)$modTable
  colmasq=table[,2]>min_posterior
  colmasq[1]=FALSE ##Intercept
  xdata=mcmc$data$X
  if (is.null(names)==FALSE & length(names)==ncol(xdata)-1)
  {
    colnames(xdata)=c("Intercept",names)
  }
  xdata=xdata[,colmasq]
  return(glm(mcmc$data$y ~ xdata,family=binomial()))
}

mcmcdetercoef=function(mcmc,min_posterior=0.8) {
  table=summary(mcmc)$modTable
  colmasq=table[,2]>min_posterior
  colmasq[1]=TRUE ##Intercept
  coeffs=table[colmasq,1]
  x=mcmc$data$X[,colmasq]
  yobs=mcmc$data$y
  ypred=1/(1+exp(-x%*%coeffs))
  return(RsqGLMdm(obs = yobs,pred=ypred))
}

mcmclatextable=function(mcmc,outfile,min_posterior=0.8,names) {
  require(xtable)
  table=summary(mcmc)$modTable
  colmasq=table[,2]>min_posterior
  colmasq[1]=TRUE ##Intercept
  if (length(rownames(table)) != length(names)+1 ) {
    stop("The names to construct the output latex table of the Bayesian variable selection do not have the appropriate length")
  }
  rownames(table)=c(rownames(table)[1],names)
  table=table[colmasq,]
  table=cbind(rownames(table),table)
  colnames(table)[1]="Regressors"
  sink(outfile)
  cat("\\documentclass{article}
\\usepackage[landscape]{geometry}
\\usepackage{adjustbox}

\\begin{document}
\\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
\\begin{tabular}{lllll}")
  print(xtable(table,digits=3),include.rownames=FALSE,only.contents=TRUE)
  cat(paste("Coefficient of discrimination= ",round(mcmcdetercoef(mcmc,min_posterior)$Tjur,digits = 3)))
  cat("\\end{tabular}
\\end{adjustbox}
\\end{document}")
  sink()
}

lassotextable=function(cvfit,outfile="outfile.tex",lambda=c("lambda.1se","lambda.min")) {
  table=coef(cvfit,s=lambda)
  table=as.matrix(table[table[,1]!=0,])
  ncv=which(cvfit$lambda==cvfit[lambda])
  missclass=paste(round(x=cvfit$cvm[ncv],digits = 3),"$\\pm$",round(x=cvfit$cvsd[ncv],digits = 3))
  
  sink(outfile)
  cat("\\documentclass{article}
\\usepackage[landscape]{geometry}
\\usepackage{adjustbox}

\\begin{document}
\\begin{adjustbox}{width={\\textwidth},totalheight={\\textheight},keepaspectratio}%
\\begin{tabular}{ll}
Regressor & Coefficient\\\\")
  print(xtable(table,digits=3),include.rownames=TRUE,only.contents=TRUE,include.colnames=FALSE)
  cat(paste("Missclassification error= ",missclass))
  cat("\\end{tabular}
\\end{adjustbox}
\\end{document}")
  sink()
}


multiplelasso=function(y,x,sx,type.measure="class",nfolds=100) {
  require(glmnet)
  require(doMC)
  registerDoMC(cores=4)
  big=glmnet(x,y, family="binomial")
  small=glmnet(sx,y, family="binomial")
  cvsmall=cv.glmnet(sx, y, family="binomial", type.measure = type.measure, nfolds = nfolds, parallel = TRUE)
  cvbig=cv.glmnet(x, y, family="binomial", type.measure = type.measure, nfolds = nfolds, parallel = TRUE)
  return(list(biglasso=big,smalllasso=small,cvbig=cvbig,cvsmall=cvsmall))
}


###HARDCODED VARIABLES IN THIS FUNCTION!!!
preparedata=function(data) {
  data=data[,-6]
  data=data[complete.cases(data),]
  data=data[data$rfdif != 0,]
  x=data[,c(4,5,7:22)]
  x$random=as.numeric(x$missingmode=="random")
  x=x[,-2]
  y=data$rfdif
  yworst=as.numeric(data$rfdif<0)
  x$pyule=(1-exp(-x$SB_rate*x$F_length_gen))^x$Leaves ##Probability of the yule process. This is not correct, since F_length_gen is the tMRCA instead of tOrigin. They should be strongly correlated, though
  x$emut=x$Mu*x$Length_cu*x$Ne/(x$Leaves*2) ##Expected number of mutations per internal branch. It should be (leaves-1)*2, but we have an extra leave, the outgroup 
  x=as.matrix(x)
  yworst=as.vector(yworst)
  scaledx=scale(x)
  smallx=x[,c(1,2,5,6,8,9,10,11,14,15,18,19,20)]
  scaledsmallx=scale(smallx)
  colnames(scaledsmallx)=c("MissingData","Species","HeightCu","LengthCu","Individuals","Loci","SpeciesRateHeterogeneity","GeneRateHeterogeneity","PropExtraLineagesM","PropExtraLineagesSD","Random","ProbYule","MutperBranch")
  colnames(scaledx)=c("MissingData","Species","SpeciationRate","TreeHeightGen","HeightCu","LengthCu","Outgroup","Individuals","Loci","SpeciesRateHeterogeneity","GeneRateHeterogeneity","Ne","Mu","PropExtraLineagesM","PropExtraLineagesSD","ExtraLineagesM","ExtraLineagesSD","Random","ProbYule","MutperBranch")
  return(list(yworst=yworst,scaledx=scaledx,scaledsmallx=scaledsmallx))
}



##Data loading
njst_data=read.csv("rf.csv")
njst_times=read.csv("time.stats",sep = " ")
njst_times=njst_times[complete.cases(njst_times),]
sqlite    <- dbDriver("SQLite")
db=dbConnect(sqlite,"sim_broadsims2.db")

###Looking for missing data and fixing it

for (method in unique(njst_data$method)) {
  for (data in unique(njst_data$mdata)) {
    if (nrow(njst_data[njst_data$method==method & njst_data$mdata==data,]) != 10000) {
      toadd=c(1:10000)[!c(1:10000) %in% njst_data[njst_data$method==method & njst_data$mdata==data,]$rep]
      print(paste("Missing data in ",method,data,":",toadd))
      for (reptoadd in toadd) {
        njst_data=rbind(njst_data,c(reptoadd,NA,data,method))
      }
    }
  }
}

##Data reorganization
print("This will generate NAs, it is OK")
njst_data$pmissing=rep(0,nrow(njst_data))
njst_data$missingmode=factor(x=rep(0,nrow(njst_data)),labels=c("random","taxa"),levels=c(0,1))
for (i in levels(njst_data$mdata))
{
    v=as.numeric(sub("r._(.*)","\\1",x=i))*0.01
    if (is.na(v))
    {
      v=0
    }
    njst_data[njst_data$mdata==i,]$pmissing=v
    if(sub("r(.)_.*","\\1",x=i)=="1")
    {
      njst_data[njst_data$mdata==i,]$missingmode="taxa"
    }
    else
    {
      njst_data[njst_data$mdata==i,]$missingmode="random"
    }
}
print("Generating NAs may stop being OK")
njst_data$rep=as.numeric(njst_data$rep)
njst_data$rf=as.numeric(njst_data$rf)

##Importing data from DB
Species_Trees=dbGetQuery(db,paste0("select ",paste(dbListFields(db, "Species_Trees"),collapse=",")," from Species_Trees"))
##Locus_Trees=dbGetQuery(db,paste0("select ",paste(dbListFields(db, "Locus_Trees"),collapse=",")," from Locus_Trees")) #Nothing interesting here
Gene_Trees=dbGetQuery(db,paste0("select ",paste(dbListFields(db, "Gene_Trees"),collapse=",")," from Gene_Trees"))
Gene_Trees$maxextral=choose(Gene_Trees$n_leaves-1,2)
Gene_Trees$propextral=Gene_Trees$Extra_l/Gene_Trees$maxextral
propextral=aggregate(propextral ~ SID, FUN=mean,data=Gene_Trees)
names(propextral)=c("SID","propextra_lm")
propextra_lsd=aggregate(propextral ~ SID, FUN=sd,data=Gene_Trees)
names(propextra_lsd)=c("SID","propextra_lsd")
extra_lm=aggregate(Extra_l ~ SID, FUN=mean,data=Gene_Trees)
names(extra_lm)=c("SID","extra_lm")
extra_lsd=aggregate(Extra_l ~ SID, FUN=sd,data=Gene_Trees)
names(extra_lsd)=c("SID","extra_lsd")
addata=merge(merge(propextral,propextra_lsd,by="SID"),merge(extra_lm,extra_lsd,by="SID"),by="SID")
datadb=merge(Species_Trees,addata)

##Hardcoded##
datadb=datadb[,c(1,2,3,5,6,7,8,9,10,15,16,18,19,24,25,26,27)] ##I eliminate here variables that do not change in this simulation study and those that affect at the sequence level (we are using true trees, so this does not modify anything)
names(datadb)=c("rep",names(datadb)[-1])
###

##Figure 1: NJst vs NJstmo
##########################
##########################

list1=compare2_plot(x=njst_data,methods=c("lnjst","onjst"),names=c("NJst","NJstmo"),glm_option = "poly",normalized = FALSE)
plot1=list1[[1]]
njst_dataf1=list1[[2]]
table1=list1[[3]]
save_plot(plot1+labs(title="NJstmo vs NJstm"),filename = "plot1.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot1+labs(title="NJstmo vs NJstm"),filename = "plot1.png",base_height = 6,base_aspect_ratio = 1.6)

##LASSO
##If we use all data, draws support the 0 variable and make the estimation depend a lot on that (which conditons generate more or less draws, instead of if one is better or worse than the other)
##We could use multinomial regression, but the multiple coefficients per predictor are difficult to plot and understand.
##I decided to remove draws.

totaldatam1=merge(njst_dataf1,datadb)
list1r=preparedata(totaldatam1)
y1=list1r$yworst
x1=list1r$scaledx
sx1=list1r$scaledsmallx

#Time consuming, but much less than the MCMC
lasso1=multiplelasso(y1,x1,sx1)
############################################

##Taking a look to the results
plot(lasso1$biglasso,xvar="dev",label=TRUE)
plot(lasso1$smalllasso,xvar="dev",label=TRUE)
coef(lasso1$smalllaso, s = "lambda.min")
coef(lasso1$smalllaso, s = "lambda.1se")

lassotextable(lasso1$cvsmall,outfile = "LASSO1.tex",lambda = "lambda.min")

coef(lasso1$cvsmall,s="lambda.1se")

predy=predict(lasso1$cvsmall,newx=sx1,s= "lambda.min",type="class")
mean(y1[predy==y1])*length(y1[predy==y1]) ##Number of correctly predicted replicates in which njstmo was worse
mean(y1==1)*length(y1) ##Number of replicates in which njstmo was worse
mean(!y1[predy==y1])*length(y1[predy==y1])
mean(y1==0)*length(y1)
#Bayesian variable selection

###TIME CONSUMING####
mcmc=logitBvs(y1,rep(1,length(y1)),sx1)
######################

plot(mcmc)
summary(mcmc,IAT = TRUE)

mcmclatextable(mcmc = mcmc,outfile = "mcmc1.tex",names=colnames(sx1))


##Figure 2: NJstmo vs NJstmou
##########################
##########################
list2=compare2_plot(x=njst_data,methods=c("onjst","unjst"),names=c("NJstmo","NJstmu"),glm_option = "glm",normalized = FALSE)
plot2=list2[[1]]
njst_dataf2=list2[[2]]
table2=list2[[3]]
save_plot(plot2+labs(title="NJstmu vs NJstmo"),filename = "plot2.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot2+labs(title="NJstmu vs NJstmo"),filename = "plot2.png",base_height = 6,base_aspect_ratio = 1.6)

#LASSO

totaldatam2=merge(njst_dataf2,datadb)
list2r=preparedata(totaldatam2)
y2=list2r$yworst
x2=list2r$scaledx
sx2=list2r$scaledsmallx

lasso2=multiplelasso(y2,x2,sx2)
predy2=predict(lasso2$cvsmall,newx=sx2,s= "lambda.min",type="class")
mean(y2[predy2==y2])*length(y2[predy2==y2]) ##Number of correctly predicted replicates in which njstmu was worse
mean(y2==1)*length(y2) ##Number of replicates in which njstmu was worse
mean(!y2[predy2==y2])*length(y2[predy2==y2])
mean(y2==0)*length(y2)
##Taking a look to the results
plot(lasso2$biglasso,xvar="dev",label=TRUE)
plot(lasso2$smalllasso,xvar="dev",label=TRUE)
coef(lasso2$smalllaso, s = "lambda.min")
coef(lasso2$smalllaso, s = "lambda.1se")

lassotextable(lasso2$cvsmall,outfile = "LASSO2.tex",lambda = "lambda.min")

#Bayesian variable selection

###TIME CONSUMING####
mcmc2=logitBvs(y2,rep(1,length(y2)),sx2)
######################

##Taking a look to the results
plot(mcmc2)
summary(mcmc2,IAT = TRUE)

mcmclatextable(mcmc = mcmc2,outfile = "mcmc2.tex",names=colnames(sx2))
n_chunks=100

datails=totaldatam2[totaldatam2$rfdif!=0 & totaldatam2$mdata!="original_g_trees",]
#datails=totaldatam2
test=seq(from=min(datails$propextra_lm),to=max(datails$propextra_lm),length.out = n_chunks)
tocenter=(test[2]-test[1])/2

dat=data.frame(missingmode=vector(mode = "character",length = (n_chunks-1)*length(unique(datails$mdata))),pmissing=vector(mode="double",length = (n_chunks-1)*length(unique(datails$mdata))),x=vector(mode="double",length = (n_chunks-1)*length(unique(datails$mdata))),mean=vector(mode="double",length = (n_chunks-1)*length(unique(datails$mdata))),stringsAsFactors = FALSE)

i=1
for (missing in unique(datails$pmissing))
{
  for (type in unique(datails$missingmode)) {
    for (ival in 2:(length(test))) {
      val=test[ival]
      dat[i,]=c(type,missing,val-tocenter,mean(datails[datails$pmissing==missing &datails$missingmode==type & datails$propextra_lm>test[ival-1] & datails$propextra_lm<val,]$rfdif))
      i=i+1  
    }
  }
}
dat$missingmode=as.factor(dat$missingmode)
dat$pmissing=as.numeric(dat$pmissing)
dat$x=as.numeric(dat$x)
dat$mean=as.numeric(dat$mean)
plotnjstmomu=ggplot(data=dat,aes(x=x,y=mean,color=pmissing,linetype=missingmode,group=interaction(missingmode,pmissing)))+geom_point()+geom_smooth(method = "loess",se = FALSE)+scale_color_distiller(palette = "Spectral",name="% missing data",breaks=c(0.10,0.25,0.50,0.75))+scale_y_continuous(name="Mean RF difference between NJstmo and NJstmu")+scale_x_continuous(name="Relative number of extra lineages")+scale_linetype(name="Missing data type",labels=c("Random","By-individual"))
save_plot(plotnjstmomu+labs(title="NJstmu vs NJstmo"),filename = "plot2b.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plotnjstmomu+labs(title="NJstmu vs NJstmo"),filename = "plot2b.png",base_height = 6,base_aspect_ratio = 1.6)

##Number of extra lineages for the mean tree
choose(90,2)*0.07

##ASTRID

##Figure 3
list3=compare2_plot(x=njst_data,methods=c("astriddef","astridmo"),names=c("ASTRID","ASTRIDmo"),glm_option = "glm",normalized = FALSE)
plot3=list3[[1]]
astrid_data=list3[[2]]
table3=list3[[3]]
save_plot(plot3+labs(title="ASTRIDmo vs ASTRID"),filename = "plot3.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot3+labs(title="ASTRIDmo vs ASTRID"),filename = "plot3.png",base_height = 6,base_aspect_ratio = 1.6)

#LASSO

totaldatam3=merge(astrid_data,datadb)

list3r=preparedata(totaldatam3)
y3=list3r$yworst
x3=list3r$scaledx
sx3=list3r$scaledsmallx

lasso3=multiplelasso(y3,x3,sx3)

##Taking a look to the results
plot(lasso3$biglasso,xvar="dev",label=TRUE)
plot(lasso3$smalllasso,xvar="dev",label=TRUE)
coef(lasso3$smalllaso, s = "lambda.min")
coef(lasso3$smalllaso, s = "lambda.1se")

lassotextable(lasso3$cvsmall,outfile = "LASSO3.tex",lambda = "lambda.min")

#Bayesian variable selection

###TIME CONSUMING####
mcmc3=logitBvs(y3,rep(1,length(y3)),sx3)
######################

##Taking a look to the results
plot(mcmc3)
summary(mcmc3,IAT = TRUE)

mcmclatextable(mcmc = mcmc3,outfile = "mcmc3.tex",names=colnames(sx3))

##ASTRID

##Figure 4
list4=compare2_plot(x=njst_data,methods=c("astridmo","astridmu"),names=c("ASTRIDmo","ASTRIDmu"),glm_option = "glm",normalized = FALSE)
plot4=list4[[1]]
astrid_data2=list4[[2]]
table4=list4[[3]]
save_plot(plot4+labs(title="ASTRIDmu vs ASTRIDmo"),filename = "plot4.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot4+labs(title="ASTRIDmu vs ASTRIDmo"),filename = "plot4.png",base_height = 6,base_aspect_ratio = 1.6)


#LASSO

totaldatam4=merge(astrid_data2,datadb)

list4r=preparedata(totaldatam4)
y4=list4r$yworst
x4=list4r$scaledx
sx4=list4r$scaledsmallx

lasso4=multiplelasso(y4,x4,sx4)

##Taking a look to the results
plot(lasso4$biglasso,xvar="dev",label=TRUE)
plot(lasso4$smalllasso,xvar="dev",label=TRUE)
coef(lasso4$smalllaso, s = "lambda.min")
coef(lasso4$smalllaso, s = "lambda.1se")

lassotextable(lasso4$cvsmall,outfile = "LASSO4.tex",lambda = "lambda.min")

#Bayesian variable selection

###TIME CONSUMING####
mcmc4=logitBvs(y4,rep(1,length(y4)),sx4)
######################

##Taking a look to the results
plot(mcmc4)
summary(mcmc4,IAT = TRUE)

mcmclatextable(mcmc = mcmc4,outfile = "mcmc4.tex",names=colnames(sx4))


##Figure 5
list5=compare2_plot(x=njst_data,methods=c("astral","astridmu"),names=c("ASTRAL2multiind","ASTRIDmu"),glm_option = "glm",normalized = FALSE)
plot5=list5[[1]]+labs(title="ASTRIDmu vs ASTRAL2")
astralastrid_data=list5[[2]]
table5=list5[[3]]
save_plot(plot5,filename = "plot5.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot5,filename = "plot5.png",base_height = 6,base_aspect_ratio = 1.6)

isc=read.csv("isc.csv",header=FALSE)
colnames(isc)=c("rep","nmissinspecies","mdata")
final_njst_data=merge(njst_data,isc)

list5b=compare2_plot(x=final_njst_data[final_njst_data$nmissinspecies==0,],methods=c("astral","astridmu"),names=c("ASTRAL2multiind","ASTRIDmu"),glm_option = "glm",normalized = FALSE)
plot5b=list5b[[1]]+labs(title="ASTRIDmu vs ASTRAL2 w/o missing sp combs")
astralastrid_datab=list5b[[2]]
table5b=list5b[[3]]
save_plot(plot5b,filename = "plot5b.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot5b,filename = "plot5b.png",base_height = 6,base_aspect_ratio = 1.6)

testdatam5c=merge(merge(njst_data,datadb),isc)

list5c=compare2_plot(x=testdatam5c[testdatam5c$Ind_per_sp>=2,],methods=c("astral","astridmu"),names=c("ASTRAL2multiind","ASTRIDmu"),glm_option = "glm",normalized = FALSE)
plot5c=list5c[[1]]+labs(title="ASTRIDmu vs ASTRAL2 with multiple individuals per species")
astralastrid_datac=list5c[[2]]
table5c=list5c[[3]]
save_plot(plot5c,filename = "plot5c.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot5c,filename = "plot5c.png",base_height = 6,base_aspect_ratio = 1.6)

list5d=compare2_plot(x=testdatam5c[testdatam5c$Ind_per_sp==3 & testdatam5c$nmissinspecies==0,],methods=c("astral","astridmu"),names=c("ASTRAL2multiind","ASTRIDmu"),glm_option = "glm",normalized = FALSE)
plot5d=list5d[[1]]+labs(title="ASTRIDmu vs ASTRAL2 with multiple ind per sp w/o missing sp combs")
astralastrid_datad=list5d[[2]]
table5d=list5d[[3]]
save_plot(plot5d,filename = "plot5d.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot5d,filename = "plot5d.png",base_height = 6,base_aspect_ratio = 1.6)
#LASSO

totaldatam5=merge(astralastrid_data,datadb)

list5r=preparedata(totaldatam5)
y5=list5r$yworst
x5=list5r$scaledx
sx5=list5r$scaledsmallx

lasso5=multiplelasso(y5,x5,sx5)

##Taking a look to the results
plot(lasso5$biglasso,xvar="dev",label=TRUE)
plot(lasso5$smalllasso,xvar="dev",label=TRUE)
coef(lasso5$smalllaso, s = "lambda.min")
coef(lasso5$smalllaso, s = "lambda.1se")

lassotextable(lasso5$cvsmall,outfile = "LASSO5.tex",lambda = "lambda.min")

#Bayesian variable selection

###TIME CONSUMING####
mcmc5=logitBvs(y5,rep(1,length(y5)),sx5)
######################

##Taking a look to the results
plot(mcmc5)
summary(mcmc5,IAT = TRUE)

mcmclatextable(mcmc = mcmc5,outfile = "mcmc5.tex",names=colnames(sx5))


##Without ISC

totaldatam5b=merge(astralastrid_datab,datadb)

list5br=preparedata(totaldatam5b)
y5b=list5br$yworst
x5b=list5br$scaledx
sx5b=list5br$scaledsmallx

lasso5b=multiplelasso(y5b,x5b,sx5b)

##Taking a look to the results
plot(lasso5b$biglasso,xvar="dev",label=TRUE)
plot(lasso5b$smalllasso,xvar="dev",label=TRUE)
coef(lasso5b$smalllaso, s = "lambda.min")
coef(lasso5b$smalllaso, s = "lambda.1se")

lassotextable(lasso5b$cvsmall,outfile = "LASSO5b.tex",lambda = "lambda.min")

#Bayesian variable selection

###TIME CONSUMING####
mcmc5b=logitBvs(y5b,rep(1,length(y5b)),sx5b)
######################

##Taking a look to the results
plot(mcmc5b)
summary(mcmc5b,IAT = TRUE)

mcmclatextable(mcmc = mcmc5b,outfile = "mcmc5b.tex",names=colnames(sx5b))

plot5grid=plot_grid(plot5,plot5c,labels = c("A","B"),hjust = 0,nrow = 2)

save_plot(plot5grid,filename = "plot5grid.pdf",base_height = 11,base_aspect_ratio = 1)
save_plot(plot5grid,filename = "plot5grid.png",base_height = 11,base_aspect_ratio = 1)


nlevels=length(unique(totaldatam5$Ind_per_sp))
datindails=totaldatam5[complete.cases(totaldatam5),]
i=1
datind=data.frame(missingmode=vector(mode = "character",length = nlevels*3*length(unique(datindails$mdata))),pmissing=vector(mode="double",length = nlevels*3*length(unique(datindails$mdata))),nind=vector(mode="double",length = nlevels*3*length(unique(datindails$mdata))),type=vector(mode="character",length = nlevels*3*length(unique(datindails$mdata))),y=vector(mode="double",length = nlevels*3*length(unique(datindails$mdata))),stringsAsFactors = FALSE)
for (missing in unique(datindails$pmissing))
{
  for (type in unique(datindails$missingmode)) {
    for (met in c("ASTRAL","Draw","ASTRIDmu")){
      for (nind in 1:nlevels) {
        #datind[i,]=c(type,missing,nind,mean(datindails[datindails$pmissing==missing &datindails$missingmode==type & datindails$Ind_per_sp==nind,]$rfdif<0,na.rm=TRUE))
        val="Draw"
        if (met == "ASTRAL") {
          datind[i,]=c(type,missing,nind,met,mean(datindails[datindails$pmissing==missing &datindails$missingmode==type & datindails$Ind_per_sp==nind,]$rfdif<0,na.rm=TRUE))
        }else if (met=="ASTRIDmu") {
          datind[i,]=c(type,missing,nind,met,mean(datindails[datindails$pmissing==missing &datindails$missingmode==type & datindails$Ind_per_sp==nind,]$rfdif>0,na.rm=TRUE))
        }else {
          datind[i,]=c(type,missing,nind,met,mean(datindails[datindails$pmissing==missing &datindails$missingmode==type & datindails$Ind_per_sp==nind,]$rfdif==0,na.rm=TRUE))
        }
        
        i=i+1  
      }
    }
  }
}
datind$missingmode=as.factor(datind$missingmode)
datind$pmissing=as.numeric(datind$pmissing)
datind$nind=as.numeric(datind$nind)
#datind$mean=as.numeric(datind$mean)
datind$y=as.numeric(datind$y)
datind$type=as.factor(datind$type)

plotastralastridnindastral=ggplot(data=datind[datind$type=="ASTRAL",],aes(x=nind,y=y,color=pmissing,linetype=missingmode,group=interaction(missingmode,pmissing)))+geom_point()+geom_smooth(method="glm",formula=y~poly(x,3),se=FALSE)+scale_color_distiller(palette = "Spectral",name="% missing data",breaks=c(0.10,0.25,0.50,0.75),guide=FALSE)+scale_y_continuous(name="Probability of outcome (win)",limits = c(0,1))+scale_x_continuous(name="Individuals per species")+scale_linetype(name="Missing data type",labels=c("Random","By-individual"),guide=FALSE)+labs(title="ASTRAL2")
plotastralastridninddraw=ggplot(data=datind[datind$type=="Draw",],aes(x=nind,y=y,color=pmissing,linetype=missingmode,group=interaction(missingmode,pmissing)))+geom_point()+geom_smooth(method="glm",formula=y~poly(x,3),se=FALSE)+scale_color_distiller(palette = "Spectral",name="% missing data",breaks=c(0.10,0.25,0.50,0.75),guide=FALSE)+scale_y_continuous(name="Probability of outcome (win)",limits = c(0,1))+scale_x_continuous(name="Individuals per species")+scale_linetype(name="Missing data type",labels=c("Random","By-individual"),guide=FALSE)+labs(title="Draw")
plotastralastridnindastrid=ggplot(data=datind[datind$type=="ASTRIDmu",],aes(x=nind,y=y,color=pmissing,linetype=missingmode,group=interaction(missingmode,pmissing)))+geom_point()+geom_smooth(method="glm",formula=y~poly(x,3),se=FALSE)+scale_color_distiller(palette = "Spectral",name="% missing data",breaks=c(0.10,0.25,0.50,0.75))+scale_y_continuous(name="Probability of outcome (win)",limits = c(0,1))+scale_x_continuous(name="Individuals per species")+scale_linetype(name="Missing data type",labels=c("Random","By-individual"))+labs(title="ASTRIDmu")
plotastralastridindfinal=plot_grid(plotastralastridnindastral,plotastralastridninddraw,plotastralastridnindastrid,nrow=1,rel_widths = c(1,1,1.4),labels = c("A","B","C"))
save_plot(plotastralastridindfinal,filename="plotastralastridind.pdf",base_height=5,base_aspect_ratio=3)
save_plot(plotastralastridindfinal,filename="plotastralastridind.png",base_height=5,base_aspect_ratio=3)
##ggplot(data=datind,aes(x=nind,fill=type,y=y))+geom_area(position="stack")+facet_wrap(~ missingmode +pmissing,nrow=2,ncol=5)

#plotastralastridnind=ggplot(data=datind,aes(x=nind,y=mean,color=pmissing,linetype=missingmode,group=interaction(missingmode,pmissing)))+geom_point()+geom_smooth(method = "loess",se = FALSE)+scale_color_distiller(palette = "Spectral",name="% missing data",breaks=c(0.10,0.25,0.50,0.75))+scale_y_continuous(name="Mean RF difference between ASTRAL2 and ASTRIDmu")+scale_x_continuous(name="Relative number of extra lineages")+scale_linetype(name="Missing datinda type",labels=c("Random","By-individual"))

##Times
timesfinal=merge(njst_times,datadb)
timeloci=ggplot(data=timesfinal[timesfinal$method=="astral" | timesfinal$method=="astridmu",],aes(x=N_loci,y=time,shape=method,color=method))+stat_summary(fun.y = mean,fun.ymin = function(x) mean(x) - sd(x),fun.ymax = function(x) mean(x) + sd(x))+geom_smooth(method = "glm",aes(color=NULL),color="black")+scale_y_log10("CPU time (s, logscale)")+scale_color_brewer(name="Method",type = "qual",palette=2,labels=c("ASTRAL2","ASTRIDmu"))+scale_x_continuous(name="Number of loci")+scale_shape_discrete(name="Method",labels=c("ASTRAL2","ASTRIDmu"))
timeleaves=ggplot(data=timesfinal[timesfinal$method=="astral" | timesfinal$method=="astridmu",],aes(x=Leaves*Ind_per_sp,y=time,shape=method,color=method))+stat_summary(fun.y = mean,fun.ymin = function(x) mean(x) - sd(x),fun.ymax = function(x) mean(x) + sd(x))+geom_smooth(method = "glm",aes(color=NULL),color="black")+scale_y_log10("CPU time (s, logscale)")+scale_color_brewer(name="Method",type = "qual",palette=2,labels=c("ASTRAL2","ASTRIDmu"))+scale_x_continuous(name="Number of gene tree leaves")+scale_shape_discrete(name="Method",labels=c("ASTRAL2","ASTRIDmu"))
timegrid=plot_grid(timeloci,timeleaves,labels = c("A","B"),hjust = 0,ncol = 2)
save_plot(timegrid,filename="times.png",base_height=4,base_aspect_ratio=3.2)

###EXTRA NJST and ASTRID original vs ASTRAL2
list6=compare2_plot(x=njst_data,methods=c("astral","lnjst"),names=c("ASTRAL2multiind","NJst"),glm_option = "poly",normalized = FALSE)
plot6=list6[[1]]
astralastrid_data=list6[[2]]
table6=list6[[3]]
save_plot(plot6+labs(title="NJst vs ASTRAL2"),filename = "plot6.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot6+labs(title="NJst vs ASTRAL2"),filename = "plot6.png",base_height = 6,base_aspect_ratio = 1.6)

list7=compare2_plot(x=njst_data,methods=c("astral","astriddef"),names=c("ASTRAL2multiind","ASTRID"),glm_option = "poly",normalized = FALSE)
plot7=list7[[1]]
astralastrid_data=list7[[2]]
table7=list7[[3]]
save_plot(plot7+labs(title="ASTRID vs ASTRAL2"),filename = "plot7.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot7+labs(title="ASTRID vs ASTRAL2"),filename = "plot7.png",base_height = 6,base_aspect_ratio = 1.6)

list8=compare2_plot(x=njst_data,methods=c("astral","unjst"),names=c("ASTRAL2multiind","NJstmu"),glm_option = "poly",normalized = FALSE)
plot8=list8[[1]]
astralastrid_data=list8[[2]]
table8=list8[[3]]
save_plot(plot8+labs(title="NJstmu vs ASTRAL2"),filename = "plot8.pdf",base_height = 6,base_aspect_ratio = 1.6)
save_plot(plot8+labs(title="NJstmu vs ASTRAL2"),filename = "plot8.png",base_height = 6,base_aspect_ratio = 1.6)

