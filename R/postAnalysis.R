library(ggplot2)
library(cowplot)
library(geiger)
library(mvMORPH)
library(phytools)
library(ape)
library(OneR)
library(viridis)
library(boot)
library(nlme)



CleanData <- function(phy, data) {
  cleaned <- treedata(phy,data,warnings=T)
  return(cleaned)
}

compareERC<-function(model.results,erc.file)
{
  erc <- read.table(erc.file,sep="\t",header=T,stringsAsFactors = F)
  tmp.1 <- merge(model.results,erc,by=c("Prot.1","Prot.2"))
  tmp.2 <- merge(model.results,erc,by.x=c("Prot.1","Prot.2"),by.y=c("Prot.2","Prot.1"))
  tmp <- rbind(tmp.1,tmp.2)
  return(tmp)
}

getUniqueMembers <- function(model.results,size.cutoff)
{
  genes <- unique(c(model.results$Prot.1,model.results$Prot.2))
  total <- 0
  for (g in genes)
  {
    tmp <- model.results[which((model.results$Prot.1 == g) | (model.results$Prot.2 == g)),]
    possible.pairs <- rownames(tmp)
    if (length(possible.pairs) > 1)
    { 
      to.remove <- sample(possible.pairs,size = length(possible.pairs) - 1,replace = F)
      model.results <- model.results[!rownames(model.results) %in% to.remove,]
    }
  }
  size <- nrow(model.results)
  if (size > size.cutoff)
  {
    to.remove <- sample(rownames(model.results),size = size - size.cutoff,replace = F)
    model.results <- model.results[!rownames(model.results) %in% to.remove,]
  }
  return(model.results)
}


getPerformanceMetric <- function(df,method="phylo",metric="pvalue")
{
  df.test <- df[which(df$Linked == "binding"),]
  df.null <- df[which(df$Linked == "None"),]
  if (nrow(df.null) > nrow(df.test))
  {
    df.null <- df.null[sample(1:nrow(df.null),size = nrow(df.test)),]
  } else{
    df.test <- df.test[sample(1:nrow(df.test),size = nrow(df.null)),]
  }
  df2 <- rbind(df.test,df.null)
  df2[,"LRT.pval.BH"] <- p.adjust(df2$LRT.pval,method="BH")
  df2[,"Pearson.pval.BH"] <- p.adjust(df2$Corr.significant,method="BH")
  df.test <- df2[which(df2$Linked == "binding"),]
  df.null <- df2[which(df2$Linked == "None"),]
  if (method == "phylo")
  {
    if (metric == "pvalue")
    {
      tp <- length(which(df.test$LRT.pval.BH < 0.05))
      fn <- nrow(df.test) - tp
      fp <- length(which(df.null$LRT.pval.BH < 0.05))
      tn <- nrow(df.null) - fp
      tpr <- tp/(tp+fn)
      tnr <- tn/(tn+fp)
      fdr <- fp/(fp+tp)
      acc <- (tp+tn)/(tp+tn+fp+fn)
      perf <- list(TPR=tpr,FPR=1-tnr,FDR=fdr,Accuracy=acc)
      return(perf)
    }
  }else if (method == "pearson")
  {
    if (metric == "pvalue")
    {
      tp <- length(which(df.test$Pearson.pval.BH < 0.05))
      fn <- nrow(df.test) - tp
      fp <- length(which(df.null$Pearson.pval.BH < 0.05))
      tn <- nrow(df.null) - fp
      tpr <- tp/(tp+fn)
      tnr <- tn/(tn+fp)
      fdr <- fp/(fp+tp)
      acc <- (tp+tn)/(tp+tn+fp+fn)
      perf <- list(TPR=tpr,FPR=1-tnr,FDR=fdr,Accuracy=acc)
      return(perf)
    }
  }
}


checkPIC <- function(df,expr,tree)
{
  results <- data.frame(Prot.1=df[,"Prot.1"],Prot.2=df[,"Prot.2"],Check.1 = numeric(length=nrow(df)),Check.2=numeric(length=nrow(df)),stringsAsFactors = F)
  for (i in 1:nrow(df))
  {
    prot.1 <- df[i,"Prot.1"]
    prot.2 <- df[i,"Prot.2"]
    results[i,"Prot.1"] <- prot.1
    results[i,"Prot.2"] <- prot.2
    x <- expr[prot.1,]
    y <- expr[prot.2,]
    missing.species.1 <- names(x)[is.na(x)]
    missing.species.2 <- names(y)[is.na(y)]
    missing <- union(missing.species.1,missing.species.2)
    if (length(missing) != 0) next
    tmp <- tree
    x <- x[tmp$tip.label]
    y <- y[tmp$tip.label]
    #print(missing)
    # print(x)
    # print(tmp$tip.label)
    #print(tmp$tip.label)

    x.pic <- pic(x,tmp,var.contrasts = T)
    y.pic <- pic(y,tmp,var.contrasts = T)
    check.x <- cor.test(abs(x.pic[,"contrasts"]),sqrt(x.pic[,"variance"]))$p.value
    check.y <- cor.test(abs(y.pic[,"contrasts"]),sqrt(y.pic[,"variance"]))$p.value
    results[i,"Check.1"] <- check.x
    results[i,"Check.2"] <- check.y
  }
  df <- merge(df,results,by=c("Prot.1","Prot.2"))
  return(df)
}


getPIC <- function(df,expr,tree)
{
  pic.corr <- c()
  pic.sig <- c()
  for (i in 1:nrow(df))
  {
    if (df[i,"Num.missing.1"] ==0 && df[i,"Num.missing.2"] == 0)
    {
      #print("yep")
      prot.1 <- df[i,"Prot.1"]
      prot.2 <- df[i,"Prot.2"]
      x <- expr[prot.1,]
      y <- expr[prot.2,]
      x <- x[tree$tip.label]
      y <- y[tree$tip.label]
      
      x.pic <- pic(x,tree)
      y.pic <- pic(y,tree)
      cor.pic <- cor.test(x.pic,y.pic)
      pic.corr <- c(pic.corr,cor.pic$estimate)
      pic.sig <- c(pic.sig,cor.pic$p.value)
    } else
    {
      pic.corr <- c(pic.corr,NA)
      pic.sig <- c(pic.sig,NA)
    }
  }
  df$Pearson.pic <- pic.corr
  df$Pic.sig <-pic.sig
  return(df)
}

removeBadFits <- function(model.results)
{
  good <- which(model.results[,"Conv.Ind"] == 0 & model.results[,"Conv.Cor"]==0 & model.results[,"Reliable.corr"]==0 & model.results[,"Reliable.ind"]==0 & model.results[,"LRT.score"]>=0)
  model.results <- model.results[good,]
  return(model.results)
}

removeResultsWMissingData<-function(model.results,cutoff=0)
{
  model.results <- model.results[which(model.results[,"Num.missing.1"] <= cutoff & model.results[,"Num.missing.2"] <= cutoff),] 
  return(model.results)
}

pullConfidentInteraction<-function(model.results,cutoff=500)
{
  model.results <- model.results[which((model.results$Score >= cutoff & model.results$Linked!="None")|model.results$Linked=="None"),]
  return(model.results)
}


##https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
plotBack2Back <- function(model.results,main)
{
  n1 <- nrow(model.results[which(model.results$Linked == "binding"),])
  n2 <- nrow(model.results[which(model.results$Linked == "None"),])
  #model.results <- model.results[,c("Linked","Evol.cor","Pearson")]
  model.results <- model.results[which(model.results$Linked == "binding" | model.results$Linked == "None"),]
 
  n1.char <- gsub('^([0-9]+)([0-9]{3})$','\\1,\\2',as.character(n1))
  n2.char <- gsub('^([0-9]+)([0-9]{3})$','\\1,\\2',as.character(n2))
  
  evol<-(ggplot(model.results) 
         + geom_histogram(binwidth=0.1,alpha=0.4,aes(x=Pearson,y=0.1*..density..,fill=Linked),size=1.5,position="identity",color="black")
         + geom_histogram(binwidth=0.1,alpha=0.4,aes(x=Evol.cor,y=-0.1*..density..,fill=Linked),size=1.5,position="identity",color="black")
         + xlab(expression(rho[C]))
         + ylab("Relative Density")
         
         + coord_flip()
         + scale_x_continuous(sec.axis = dup_axis(name = expression(rho[U])))
         + scale_y_continuous(breaks=c(-0.2,-0.1,0,0.1,0.2),labels=c(0.2,0.10,0.00,0.1,0.2),limits=c(-0.2,0.2)) ## hardcoding this value can cause problems
         
         + ggtitle(main)
         + scale_fill_viridis(discrete=T,breaks=c("binding", "None"),labels=c(paste0("Binding (N=",n1.char,")"),paste0("Control (N=",n2.char,")")))
         + theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.text = element_text(size=10,face="bold"),
                 axis.title.y=element_text(size=14,face="bold"),
                 axis.title.x = element_text(size=14,face="bold"),
                 axis.text = element_text(size=12,face="bold"),
                 legend.title=element_blank(),
                 legend.position = c(0.2,0.1)))
  
  ymax <- ggplot_build(evol)$layout$panel_scales_x[[1]]$range$range[2]
  xmin <- ggplot_build(evol)$layout$panel_scales_y[[1]]$range$range[1]
  ymin <- ggplot_build(evol)$layout$panel_scales_x[[1]]$range$range[1]
  xmax <- ggplot_build(evol)$layout$panel_scales_y[[1]]$range$range[2]
  
  mean.cat.1 <- as.character(round(mean(model.results[which(model.results$Linked == "binding"),"Pearson"]),2))
  mean.cat.2 <- as.character(round(mean(model.results[which(model.results$Linked == "None"),"Pearson"]),2))
  
  text.pos.x <- xmax * 0.75
  text.pos.y <- ymin*0.4
  
  evol <- evol + annotate("text",x=text.pos.y,y=text.pos.x,label=sprintf("bold(bar(rho)[U])"),parse=T,size=8,fontface=2)
  
  evol <- (evol + annotate("text",x=text.pos.y - ymax*0.1,y=text.pos.x,label=sprintf("bold(Binding == %s)", mean.cat.1),parse=T,size=5,fontface=2)
           + annotate("text",x=text.pos.y - ymax*0.2,y=text.pos.x,label=sprintf("bold(Control == %s)", mean.cat.2),parse=T,size=5,fontface=2))
  
  mean.cat.1 <- as.character(round(mean(model.results[which(model.results$Linked == "binding"),"Evol.cor"]),2))
  mean.cat.2 <- as.character(round(mean(model.results[which(model.results$Linked == "None"),"Evol.cor"]),2))
  
  text.pos.x <- xmin * 0.95
  text.pos.y <- ymin*0.4
  
  
  evol <- evol + annotate("text",x=text.pos.y,y=text.pos.x,label=sprintf("bold(bar(rho)[C])"),parse=T,size=8,fontface=2)
  evol <- (evol + annotate("text",x=text.pos.y - ymax*0.1,y=text.pos.x,label=sprintf("bold(Binding == %s)", mean.cat.1),parse=T,size=5,fontface=2)
           + annotate("text",x=text.pos.y - ymax*0.2,y=text.pos.x,label=sprintf("bold(Control == %s)", mean.cat.2),parse=T,size=5,fontface=2))
  return(evol)
}


djIndex <- function(row,go_terms)
{
  go.terms.1 <- go_terms[which(go_terms$Protein == row["Prot.1"]),"GO"]
  go.terms.2 <- go_terms[which(go_terms$Protein == row["Prot.2"]),"GO"]
  dj.index <- length(intersect(go.terms.1,go.terms.2))/length(union(go.terms.1,go.terms.2))
  return(dj.index)
}

getDJIndex <- function(model.results,go_terms)
{
  model.results["DJ.Index"] <- apply(model.results,1,FUN = djIndex,go_terms=go_terms)
  return(model.results)
}

getResults<-function(path,model="BM",tree.file="../Data/fungi_tree.tre",expr.file="../Data/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv",categories=c("binding","None"),group="Linked",suffix=c("_binding_dated.tsv","_null_dated.tsv"),filter.missing=NULL,filter.score=NULL,remove.bad.go=F,check.model.violations=T)
{
  if (length(suffix) == 2)
  {

    test.case <- list.files(path=path,pattern=suffix[1],full.names = T)
    null.case <- list.files(path=path,pattern=suffix[2],full.names = T)
    print(test.case)
    print(null.case)
    df.1<- read.table(test.case,sep="\t",header=T,stringsAsFactors = F)
    df.2<- read.table(null.case,sep="\t",header=T,stringsAsFactors = F)
    df <- rbind(df.1,df.2)
  } else{
    test.case <- list.files(path=path,pattern=suffix[1],full.names = T) 
    df <-  read.table(test.case,sep="\t",header=T,stringsAsFactors = F)
  }
  print(length(which(df$Linked=="binding")))
  df <- removeBadFits(df)
  print(length(which(df$Linked=="binding")))
  x<-which(duplicated(df[,c("Prot.1","Prot.2")],fromLast = T))
  if (length(x) > 0)
  {
    df <- df[-x,]
  }
  if (!is.null(filter.score))
  {  
    df <- pullConfidentInteraction(df,filter.score)
  }
  if (!is.null(filter.missing))
  {
    df <- removeResultsWMissingData(df,filter.missing)
  }
  print(length(which(df$Linked=="binding")))
  tree <- read.tree(tree.file)
  tree <- drop.tip(tree,c("S.pombe"))
  expr <- read.table(expr.file,sep="\t",header=T,stringsAsFactors=F,row.names=1)
  df <- df[which((df$Linked == categories[1]) | (df$Linked == categories[2])),]
  if (check.model.violations)
  {
    df <- checkModelViolations(df,expr,tree)
  }
  print(length(which(df$Linked=="binding")))
  if(remove.bad.go == T)
  {
    go.terms <- read.table("../Data/go_annotation_prot_ids.txt",sep="\t",header=T)
    df <- getDJIndex(df,go.terms)
    df <- df[which((df$Linked == categories[1]) | (df$Linked == categories[2] & df$DJ.Index == 0)),]
  }
  return(df)
}

checkModelViolations <- function(df,expr,tree,model="BM")
{
  cleaned <- CleanData(tree,t(expr))
  tree <- cleaned$phy
  exp.df <- t(cleaned$data)
  df <- checkPIC(df,exp.df,tree)
  filt <- df[which(df$Check.1 > 0.05 & df$Check.2 > 0.05),]
  return(filt)
}


getGoodGenes <- function()
{
  genes <- unlist(read.table("ortholog_final_matrix_0927.tsv",sep="\t",header=T,stringsAsFactors=F)[,"scerevisiae"])
  genes <- genes[which(genes != "")]
  expr <- read.table("TPM_control/Fixed_gene_expression/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv",sep="\t",header=T,stringsAsFactors=F,row.names=1)
  tree <- read.tree("Tree/new_dated_fungi.tre")
  tree <- drop.tip(tree,"S.pombe")
  cleaned <- CleanData(tree,t(expr))
  expr <- t(cleaned$data)
  tree <- cleaned$phy
  #print(expr)
  print(genes)
  good.genes <- character(length(genes))
  for (i in 1:length(genes))
  {
    prot.1 <- genes[i]

    #print(prot.1)
    tryCatch({
    x <- expr[prot.1,]
    missing <- names(x)[is.na(x)]
    if (length(missing) != 0)
    {
      #tmp <- drop.tip(tree,missing)
      next
    } else{
      tmp <- tree
    }
    x <- x[tmp$tip.label]
    x.pic <- pic(x,tmp,var.contrasts = T)
    check.x <- cor.test(abs(x.pic[,"contrasts"]),sqrt(x.pic[,"variance"]))$p.value
    if (check.x > 0.05)
    {
      good.genes[i] <- prot.1
    }},error=function(){next})
  }
  good.genes <- good.genes[which(good.genes != "")]
  return(good.genes)
}


calculateWeights <- function(model.results,group="binding")
{
  treat <- model.results[which(model.results$Linked == group),]
  prots <- union(treat$Prot.1,treat$Prot.2)
  num.interactions <- data.frame(Protein=prots,Num.Interactions=numeric(length(prots)),Evol.cor.se=numeric(length(prots)),Theta.var=numeric(length(prots)),row.names = "Protein")
  for (i in prots)
  {
    tmp <- treat[which(treat$Prot.1 == i | treat$Prot.2 == i),c("Prot.1","Prot.2","Theta.ind.1","Theta.ind.2")]
    num.interactions[i,"Num.Interactions"] <- nrow(tmp)
  }
  treat["Log.Mean.Degrees"] <- numeric(length=nrow(treat))
  for (i in 1:nrow(treat))
  {
   num.1 <- num.interactions[treat[i,"Prot.1"],"Num.Interactions"]
   num.2 <- num.interactions[treat[i,"Prot.2"],"Num.Interactions"]
   treat[i,"Log.Mean.Degrees"] <- log((num.1 + num.2)/2)
   treat[i,"Weight"] <- 0.5*(1/num.1 + 1/num.2)
  }
  return(treat)
}


numberInteractions <- function(model.results,group="binding",dataset="newest_no_paralogs_interactions.txt")
{
  treat <- model.results[which(model.results$Linked == group),]
  inter <- read.table(dataset,header = T,sep="\t",stringsAsFactors = F)
  inter <- inter[which(inter$Type == group),]
  prots <- union(treat$Prot.1,treat$Prot.2)
  num.interactions <- data.frame(Protein=prots,Num.Interactions=numeric(length(prots)),Evol.cor.se=numeric(length(prots)),Theta.var=numeric(length(prots)),row.names = "Protein")
  for (i in prots)
  {
    tmp <- inter[which(inter$Protein_1 == i | inter$Protein_2 == i),c("Protein_1","Protein_2")]
    num.interactions[i,"Num.Interactions"] <- nrow(tmp)
  }
  treat["Log.Mean.Degrees.Full"] <- numeric(length=nrow(treat))
  for (i in 1:nrow(treat))
  {
    num.1 <- num.interactions[treat[i,"Prot.1"],"Num.Interactions"]
    num.2 <- num.interactions[treat[i,"Prot.2"],"Num.Interactions"]
    treat[i,"Log.Mean.Degrees.Full"] <- log((num.1 + num.2)/2)
  }
  return(treat)
}

wCorr.boot <- function(data,indices,w.obs)
{
  x <- data[indices,]
  w.corr <- wCorr::weightedCorr(x[,1],x[,2],method="Spearman",weights=w.obs[indices])
  return(w.corr)
}
#p-value calculation taken from https://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r

calculate.Pval.wcorr <- function(treatment,group.1,group.2,R=1000)
{
  data <- treatment[,c(group.1,group.2)]
  weights <- treatment$Weight
  boot.results <- boot(data,statistic = wCorr.boot,R=R,w.obs=weights)
  p.value <- with(boot.results,pnorm(abs((2*t0-mean(t)-0)/sd(t)),lower.tail=F)*2)
  print(boot.ci(boot.results,type="norm"))
  return(p.value)
}


plotMetricWithEvolCor<-function(model.results,metric="Mean.Anc.State",apply.arctanh=F,apply.weights=F,main,num.boot=1000,xlabel="Mean Ancestral Estimate",ylabel="Evolutionary Correlation")
{
  model.results["Mean.Anc.State"] <- (model.results$Theta.ind.1 + model.results$Theta.ind.2)/2
  model.results <- calculateWeights(model.results)
  model.results <- model.results[which(model.results$Linked == "binding"),]
  if (apply.arctanh == T)
  {
    model.results["Evol.cor"] <- atanh(model.results$Evol.cor)
  }
  if (apply.weights)
  {
    weights <- model.results$Weight
    
  } else{
    weights <- rep(1,nrow(model.results))
  }
  xmin <- min(model.results[,metric])
  xmax <- max(model.results[,metric])
  ymin <- min(model.results[,"Evol.cor"])
  ymax <- max(model.results[,"Evol.cor"])
  
  #x <- cor.test(model.results$Mean.Anc.State,model.results$Evol.cor,method="spearman")
  x <- wCorr::weightedCorr(model.results[,metric],model.results$Evol.cor,method="Spearman",weights=weights)
  p <- (ggplot(model.results,aes_string(x=metric,y="Evol.cor")) 
        +geom_point(alpha=0.3)
        + xlab(xlabel)
        + ylab(ylabel)
        + xlim(c(xmin-0.1,xmax+0.1))
        + ylim(c(ymin-0.1,ymax+0.1))
        + ggtitle(main)
        + theme(plot.title = element_text(hjust = 0.5,size=10,face="bold"),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.title.y=element_text(size=10,face="bold"),
                axis.title.x = element_text(size=10,face="bold"),
                axis.text = element_text(size=10,face="bold")))
  xmin <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1]-0.15
  xmax <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2]+0.15
  ymin <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1]-0.15
  ymax <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]+0.15
  height <- ymax - ymin
  width <- xmax - xmin
  text.pos.y <- ymax - (height * 0.04)
  text.pos.x <- xmin + (width * 0.1)
  p <- p + annotate("text",x=text.pos.x,y=text.pos.y,label=sprintf("bold(rho[S]) == %0.2f", (round(x,3))),parse=T,size=5,fontface=2)
  p.value <- calculate.Pval.wcorr(model.results,metric,"Evol.cor",R=num.boot)
  print(p.value)
  return(p)
}

getNumInteractions <- function(model.results,group="binding")
{
  if("Linked" %in% colnames(model.results))
  {
    df <- model.results[which(model.results$Linked == group),]
  } else{
    df <- model.results[which(model.results$Type == group),] 
  }
  if ("Prot.1" %in% colnames(df))
  {
    prots <- union(df$Prot.1,df$Prot.2)
  } else{
    prots <- union(df$Protein_1,df$Protein_2)
  }
  num.interactions <- numeric(length=length(prots))
  for (i in 1:length(prots))
  {
    tmp = df[which(df[,1] == prots[i] | df[,2] == prots[i]),]
    num <- nrow(tmp)
    num.interactions[i] <- num
  }
  return(num.interactions)
}


compareRealSim <- function(real,sim,category="binding",metric="Evol.cor",main,xlabel,ylabel)
{
  real <- real[which(real$Linked==category),]
  sim <- sim[which(sim$Linked==category),]
  real.sim <- merge(real,sim,by=c("Prot.1","Prot.2","Linked"),suffixes=c("_Real","_Sim"))
  metric.1 <- paste0(metric,"_Real")
  metric.2 <- paste0(metric,"_Sim")
  x <- cor.test(real.sim[,metric.1],real.sim[,metric.2],method="spearman")
  #x <- wCorr::weightedCorr(model.results[,metric],model.results$Evol.cor,method="Spearman",weights=weights)
  if (category != "None")
  {
    p <- (ggplot(real.sim,aes_string(x=metric.1,y=metric.2)) 
          +geom_point()
          + xlab(xlabel)
          + ylab(ylabel)
          + labs(title=sprintf(main))
          + theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.title.y=element_text(size=14,face="bold"),
                  axis.title.x = element_text(size=14,face="bold"),
                  axis.text = element_text(size=12,face="bold")))
    xmin <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1]
    xmax <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2]
    ymin <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1]
    ymax <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
    height <- ymax - ymin
    width <- xmax - xmin
    text.pos.y <- ymax - (height * 0.05)
    text.pos.x <- xmin + (width * 0.1)
    p <- p + annotate("text",x=text.pos.x,y=text.pos.y,label=sprintf("bold(rho[S]) == %0.2f", x$estimate),parse=T,size=5,fontface=2)
  } else{
    p <-(ggplot(real.sim,aes_string(x=metric.2)) 
           + geom_histogram(binwidth=0.1,alpha=0.4,aes(y=0.1*..density..),position="identity",color="black")
           + xlab(xlabel)
           + ylab(ylabel)
           + ggtitle(main)
           + theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.title.y=element_text(size=14,face="bold"),
                   axis.title.x = element_text(size=14,face="bold"),
                   axis.text = element_text(size=12,face="bold")))
    ymax <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]
    text.pos.y <- ymax * 0.7
    mean.rho <- mean(real.sim[which(real.sim$Linked == category),metric.2])
    sig.test <- t.test(real.sim[which(real.sim$Linked == category),metric.2])
    p.val <- sig.test$p.value
    p <- p + annotate("text",x=-0.75,y=text.pos.y,label=sprintf("bold(bar(rho)[Simulated]) == %0.2f", mean.rho),parse=T,size=5,fontface=2)
    p <- p + annotate("text",x=-0.75,y=text.pos.y-0.01,label=sprintf("p == %0.2f", p.val),parse=T,size=5,fontface=2)  
  }
  return(p)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

createSpeciesHeatMap<-function()
{
  expr <- read.table("TPM_control/Gene_Expression_Files/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv",sep="\t",header=T,row.names=1)
  expr <- expr[,-which(names(expr) %in% c("S.pombe"))]
  tree <- read.tree("Tree/new_dated_fungi.tre")
  tree <- drop.tip(tree,"S.pombe")
  cormat <- round(cor(expr,use="pairwise.complete.obs"),2)
  cormat <- cormat[tree$tip.label,tree$tip.label]
  upper_tri <- get_upper_tri(cormat)
  mid <- median(upper_tri,na.rm=T)
  cor.min <- min(unlist(upper_tri),na.rm=T)
  cor.max <- max(unlist(upper_tri),na.rm=T)
  melted_cormat <- melt(upper_tri,na.rm=T)
  p <- (ggplot(melted_cormat, aes(Var2, Var1, fill = value))
        +geom_tile(color = "white")
        +scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mid, limit = c(cor.min,cor.max), space = "Lab",  name="Pearson\nCorrelation") 
        +theme_minimal()# minimal theme
        +theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)))
  ggsave(filename="heatmap.pdf",device="pdf")
}


geneExpressionToCoevolutionDistanceMetric<-function(model.results,main,xlabel,ylabel,apply.arctanh=F,apply.weights=T,num.boot=1000)
{
  model.results[,"Dist"] <- (model.results$ERC - model.results$Evol.cor)^2
  model.results["Mean.Anc.State"] <- (model.results$Theta.ind.1 + model.results$Theta.ind.2)/2
  model.results <- calculateWeights(model.results)
  model.results <- model.results[which(model.results$Linked == "binding"),]
  if (apply.arctanh == T)
  {
    model.results["Evol.cor"] <- atanh(model.results$Evol.cor)
  }
  if (apply.weights)
  {
    weights <- model.results$Weight
    
  } else{
    weights <- rep(1,nrow(model.results))
  }
  xmin <- min(model.results[,"Mean.Anc.State"])
  xmax <- max(model.results[,"Mean.Anc.State"])
  ymin <- min(model.results[,"Dist"])
  ymax <- max(model.results[,"Dist"])
  
  #x <- cor.test(model.results$Mean.Anc.State,model.results$Evol.cor,method="spearman")
  x <- wCorr::weightedCorr(model.results[,"Mean.Anc.State"],model.results[,"Dist"],method="Spearman",weights=weights)
  p <- (ggplot(model.results,aes_string(x="Mean.Anc.State",y="Dist")) 
        +geom_point(alpha=0.3)
        + xlab(xlabel)
        + ylab(ylabel)
        + xlim(c(xmin-0.1,xmax+0.1))
        + ylim(c(ymin-0.1,ymax+0.1))
        + ggtitle(main)
        + theme(plot.title = element_text(hjust = 0.5,size=10,face="bold"),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.title.y=element_text(size=12,face="bold"),
                axis.title.x = element_text(size=12,face="bold"),
                axis.text = element_text(size=12,face="bold")))
  xmin <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1]-0.15
  xmax <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2]+0.15
  ymin <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1]-0.15
  ymax <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]+0.15
  height <- ymax - ymin
  width <- xmax - xmin
  text.pos.y <- ymax - (height * 0.1)
  text.pos.x <- xmax - (width * 0.25)
  p <- p + annotate("text",x=text.pos.x,y=text.pos.y,label=sprintf("bold(rho[S]) == %0.2f", (round(x,3))),parse=T,size=5,fontface=2)
  p.value <- calculate.Pval.wcorr(model.results,"Mean.Anc.State","Dist",R=num.boot)
  print(p.value)
  return(p)
}

fraser2004Approach <- function(sim)
{
  sim.treat <- sim[which(sim$Linked == "binding"),]
  #sim.treat <- sim.treat[-which((sim.treat$Pearson < 0 & sim.treat$Corr.significant < 0.05) | (sim.treat$Evol.cor < 0 & sim.treat$LRT.pval < 0.05)),]
  sim.null <- sim[which(sim$Linked == "None"),]
  sim.null <- sim.null[sample(1:nrow(sim.null),size=nrow(sim.treat),replace=F),]
  train.test.divide <- sample(1:nrow(sim.treat),size=0.8*nrow(sim.treat),replace = F)
  sim.treat.train <- sim.treat[train.test.divide,]
  sim.treat.test <- sim.treat[-train.test.divide,]
  
  train.test.divide <- sample(1:nrow(sim.null),size=0.8*nrow(sim.null),replace = F)
  sim.null.train <- sim.null[train.test.divide,]
  sim.null.test <- sim.null[-train.test.divide,]
  
  breaks <- seq(-1,1,0.1)
  grouped.null <- cut(x = sim.null.train$Evol.cor,breaks=breaks)
  grouped.treat <- cut(x=sim.treat.train$Evol.cor,breaks=breaks)

  evol.breaks <- summary(grouped.null)/nrow(sim.null.train)<summary(grouped.treat)/nrow(sim.treat.train)
  evol.cutoff <- NA
  for (i in 1:length(evol.breaks))
  {
    if (evol.breaks[i] == T && is.na(evol.cutoff))
    {
      evol.cutoff <- breaks[i] + 0.1
    } else if (evol.breaks[i] == F && !is.na(evol.cutoff)){
      evol.cutoff <- NA
    }
  }
  evol.power <- length(which(sim.treat.test$Evol.cor > evol.cutoff))/nrow(sim.treat.test)
  evol.fpr <- length(which(sim.null.test$Evol.cor > evol.cutoff))/nrow(sim.null.test)
  evol.fdr <- length(which(sim.null.test$Evol.cor > evol.cutoff))/(length(which(sim.treat.test$Evol.cor > evol.cutoff))+length(which(sim.null.test$Evol.cor > evol.cutoff)))
  evol.accuracy <- (evol.power*nrow(sim.treat.test) + (1-evol.fpr) * nrow(sim.null.test))/(nrow(sim.treat.test) + nrow(sim.null.test))
  
  grouped.null <- cut(x = sim.null.train$Pearson,breaks=breaks)
  grouped.treat <- cut(x=sim.treat.train$Pearson,breaks=breaks)
  pearson.breaks <- summary(grouped.null)/nrow(sim.null.train)<summary(grouped.treat)/nrow(sim.treat.train)
  pearson.cutoff <- NA
  for (i in 1:length(pearson.breaks))
  {
    if (pearson.breaks[i] == T && is.na(pearson.cutoff))
    {
      pearson.cutoff <- breaks[i] + 0.1
    } else if (pearson.breaks[i] == F && !is.na(pearson.cutoff)){
      pearson.cutoff <- NA
    }
  }
  pearson.power <- length(which(sim.treat.test$Pearson > pearson.cutoff))/nrow(sim.treat.test)
  pearson.fpr <- length(which(sim.null.test$Pearson > pearson.cutoff))/nrow(sim.null.test)
  pearson.fdr <- length(which(sim.null.test$Pearson > pearson.cutoff))/(length(which(sim.treat.test$Pearson > pearson.cutoff))+length(which(sim.null.test$Pearson > pearson.cutoff)))
 
  pearson.accuracy <- (pearson.power*nrow(sim.treat.test) + (1-pearson.fpr) * nrow(sim.null.test))/(nrow(sim.treat.test) + nrow(sim.null.test))
  return(list(Phylogenetic.cutoff=evol.cutoff,
    Pearson.cutoff=pearson.cutoff,
    Phylogenetic.power=evol.power,
    Pearson.power=pearson.power,
    Phylogenetic.fpr=evol.fpr,
    Pearson.fpr=pearson.fpr,
    Phylogenetic.fdr = evol.fdr,
    Pearson.fdr=pearson.fdr,
    Phylogenetic.accuracy=evol.accuracy,
    Pearson.accuracy=pearson.accuracy))
}

pValueCutoff <- function(x,treat.p.val,null.p.val)
{
  length(which(null.p.val < x))/(length(which(treat.p.val < x))+length(which(null.p.val < x))) - 0.05
}


martin2018Approach <- function(sim,correlation="Pearson")
{
  sim.treat <- sim[which(sim$Linked == "binding"),]
  sim.null <- sim[which(sim$Linked == "None"),]
  sim.null <- sim.null[sample(1:nrow(sim.null),size=nrow(sim.treat),replace=F),]
  null.p.val <- numeric(length=nrow(sim.null))
  treat.p.val <- numeric(length=nrow(sim.treat))
  for (i in 1:nrow(sim.null))
  {
    null.subset <- sim.null[sample(c(1:i-1,i+1:nrow(sim.null)),size=1000,replace=F),]
    null.subset.cor <- null.subset[,correlation]
    total.higher <- which(null.subset.cor > sim.null[i,correlation])
    p.value <- (1 + length(total.higher))/1001
    null.p.val[i] <- p.value
  }
  for (i in 1:nrow(sim.treat))
  {
    null.subset <- sim.null[sample(1:nrow(sim.null),size=1000,replace=F),]
    null.subset.cor <- null.subset[,correlation]
    total.higher <- which(null.subset.cor > sim.treat[i,correlation])
    p.value <- (1 + length(total.higher))/1001
    treat.p.val[i] <- p.value
  }
  fdr.cutoff <- uniroot(pValueCutoff,interval=c(0,1),treat.p.val=treat.p.val,null.p.val=null.p.val,f.lower=-0.05)
  tpr <- length(which(treat.p.val < fdr.cutoff$root))/length(treat.p.val)
  fpr <- length(which(null.p.val < fdr.cutoff$root))/length(null.p.val)
  fdr <- length(which(null.p.val < fdr.cutoff$root))/(length(which(treat.p.val < fdr.cutoff$root))+length(which(null.p.val < fdr.cutoff$root)))
  accuracy <- (tpr*nrow(sim.treat) + (1-fpr) * nrow(sim.null))/(nrow(sim.treat) + nrow(sim.null))
  return(list(TPR=tpr,
              FPR=fpr,
              FDR=fdr,
              Accuracy=accuracy))
}


printSDValues <- function(tpr,fpr,fdr,accuracy)
{
  print(paste("TPR S.D",sd(tpr)))
  print(paste("FPR S.D",sd(fpr)))
  print(paste("FDR S.D",sd(fdr)))
  print(paste("Accuracy S.D",sd(accuracy)))
}

evaluateFraser2004 <- function(sim)
{
  evol.power <- numeric(length=100)
  pearson.power <- numeric(length=100)
  evol.fpr <- numeric(length=100)
  pearson.fpr <- numeric(length=100)
  evol.fdr <- numeric(length=100)
  pearson.fdr <- numeric(length=100)
  evol.accuracy <- numeric(length=100)
  pearson.accuracy <- numeric(length=100)
  for (i in 1:100)
  {
    power <- fraser2004Approach(sim)
    evol.power[i] <- power$Phylogenetic.power
    pearson.power[i] <- power$Pearson.power
    evol.fpr[i] <- power$Phylogenetic.fpr
    pearson.fpr[i] <- power$Pearson.fpr
    evol.fdr[i] <- power$Phylogenetic.fdr
    pearson.fdr[i] <- power$Pearson.fdr
    evol.accuracy[i] <- power$Phylogenetic.accuracy
    pearson.accuracy[i] <- power$Pearson.accuracy

  }
  print("Evolutionary Correlation SD")
  printSDValues(evol.power,evol.fpr,evol.fdr,evol.accuracy)
  print("Pearson Correlation SD")
  printSDValues(pearson.power,pearson.fpr,pearson.fdr,pearson.accuracy)
  comp.power <- t.test(evol.power,pearson.power)
  comp.fpr <- t.test(evol.fpr,pearson.fpr)
  comp.fdr <- t.test(evol.fdr,pearson.fdr)
  comp.accuracy <- t.test(evol.accuracy,pearson.accuracy)
  print(comp.power)
  print(comp.fpr)
  print(comp.fdr)
  print(comp.accuracy)
}


evaluateMartin2018 <- function(sim)
{
  evol.power <- numeric(length=100)
  pearson.power <- numeric(length=100)
  evol.fpr <- numeric(length=100)
  pearson.fpr <- numeric(length=100)
  evol.fdr <- numeric(length=100)
  pearson.fdr <- numeric(length=100)
  evol.accuracy <- numeric(length=100)
  pearson.accuracy <- numeric(length=100)
  for (i in 1:100)
  {
    #print(i)
    power <- martin2018Approach(sim,"Evol.cor")
    evol.power[i] <- power$TPR
    evol.fpr[i] <- power$FPR
    evol.fdr[i] <- power$FDR
    evol.accuracy[i] <- power$Accuracy
  }
  for (i in 1:100)
  {
    #print(i)
    power <- martin2018Approach(sim,"Pearson")
    pearson.power[i] <- power$TPR
    pearson.fpr[i] <- power$FPR
    pearson.fdr[i] <- power$FDR
    pearson.accuracy[i] <- power$Accuracy

  }
  print("Evolutionary Correlation SD")
  printSDValues(evol.power,evol.fpr,evol.fdr,evol.accuracy)
  print("Pearson Correlation SD")
  printSDValues(pearson.power,pearson.fpr,pearson.fdr,pearson.accuracy)
  comp.power <- t.test(evol.power,pearson.power)
  comp.fpr <- t.test(evol.fpr,pearson.fpr)
  comp.fdr <- t.test(evol.fdr,pearson.fdr)
  comp.accuracy <- t.test(evol.accuracy,pearson.accuracy)
  print(comp.power)
  print(comp.fpr)
  print(comp.fdr)
  print(comp.accuracy)
}


evaluatePCM<- function(sim)
{
  evol.power <- numeric(length=100)
  pearson.power <- numeric(length=100)
  evol.fpr <- numeric(length=100)
  pearson.fpr <- numeric(length=100)
  evol.fdr <- numeric(length=100)
  pearson.fdr <- numeric(length=100)
  evol.accuracy <- numeric(length=100)
  pearson.accuracy <- numeric(length=100)
  for (i in 1:100)
  {
    power <- getPerformanceMetric(sim)
    evol.power[i] <- power$TPR
    evol.fpr[i] <- power$FPR
    evol.fdr[i] <- power$FDR
    evol.accuracy[i] <- power$Accuracy
  }
  for (i in 1:100)
  {
    power <- getPerformanceMetric(sim,method="pearson")
    pearson.power[i] <- power$TPR
    pearson.fpr[i] <- power$FPR
    pearson.fdr[i] <- power$FDR
    pearson.accuracy[i] <- power$Accuracy
    
  }
  print("Evolutionary Correlation SD")
  printSDValues(evol.power,evol.fpr,evol.fdr,evol.accuracy)
  print("Pearson Correlation SD")
  printSDValues(pearson.power,pearson.fpr,pearson.fdr,pearson.accuracy)
  comp.power <- t.test(evol.power,pearson.power)
  comp.fpr <- t.test(evol.fpr,pearson.fpr)
  comp.fdr <- t.test(evol.fdr,pearson.fdr)
  comp.accuracy <- t.test(evol.accuracy,pearson.accuracy)
  print(comp.power)
  print(comp.fpr)
  print(comp.fdr)
  print(comp.accuracy)
}

## Code for comparing distributions of evolutionary and peason correlations


df <- getResults("../Results/",
                 suffix=c("binding.tsv","control.tsv"),
                 group=c("Linked"),
                 check.model.violations = T,
                 filter.missing=0,
                 filter.score=NULL,
                 remove.bad.go=T,
                 model="BM"
                 )

exclude <- read.table("../Data/scer_to_exclude.txt",sep="",header=F,stringsAsFactors=F)
df <- df[which(!df[,1] %in% exclude[,1]),]
df <- df[which(!df[,2] %in% exclude[,1]),]
p<-plotBack2Back(df,"Correlation Distribution Comparing Binding Proteins\nto Randomly-Generated Control")
ggsave2("../Images/Figure2.pdf",width=165,height =165 ,units="mm")
#pdf("../Images/pc_pu_back_to_back_plot_nucl_sub_tree.pdf",width=8)
# p
# dev.off()

# # ## Code for generating correlations with evolutionary correlation
treat <- df[which(df$Linked == "binding"),]

p1 <- plotMetricWithEvolCor(treat,metric = "Score",
                           apply.arctanh = F,
                           apply.weights = T,
                           main="STRING Score vs. Phylogenetically-Corrected Correlation\nper Protein Pair",
                           xlabel="STRING Confidence Score",
                           ylabel=expression(rho[C]),
                           num.boot=1000)
ggsave2("../Images/Figure3.pdf",width=150,height =100 ,units="mm")
# p2 <- plotMetricWithEvolCor(treat,metric = "DJ.Index",
#                            apply.arctanh = F,
#                            apply.weights = T,
#                            main="Jaccard Index vs. Phylogenetically-Corrected Correlation\nper Protein Pair",
#                            xlabel="Jaccard Index",
#                            ylabel=expression(rho[C]),
#                            num.boot=1000)
p3 <- plotMetricWithEvolCor(treat,metric = "Mean.Anc.State",
                           apply.arctanh = F,
                           apply.weights = T,
                           main="Mean Ancestral Gene Expression vs. Phylogenetically-Corrected Correlation\nper Protein Pair",
                           xlabel="Mean Ancestral Gene Expression Estimate (Log scale)",
                           ylabel=expression(rho[C]),
                           num.boot=1000)
p4 <- plotMetricWithEvolCor(treat,metric = "Log.Mean.Degrees",
                           apply.arctanh = F,
                           apply.weights = T,
                           main="Mean Number of Interactions vs. Phylogenetically-Corrected Correlation\nper Protein Pair",
                           xlabel="Log(Mean Degree)",
                           ylabel=expression(rho[C]),
                           num.boot=1000)

treat.erc<-compareERC(treat,erc.file="../Data/erc_protein_ids.tsv")
#treat.erc["ERC"] <- atanh(treat.erc[,"ERC"])
p5 <- plotMetricWithEvolCor(treat.erc,metric = "ERC",
                            apply.arctanh = F,
                            apply.weights = T,
                            main="Protein Sequence Coevolution vs. Gene Expression Coevolution\nper Protein Pair",
                            xlabel="Protein Sequence Coevolution",
                            ylabel=expression(rho[C]),
                            num.boot=1000)
treat.cai<-compareERC(treat,erc.file="../Data/cai_protein_ids.tsv")
#treat.erc["ERC"] <- atanh(treat.erc[,"ERC"])
p6 <- plotMetricWithEvolCor(treat.cai,metric = "ERC",
                            apply.arctanh = F,
                            apply.weights = T,
                            main="CAI Coevolution vs. Gene Expression Coevolution\nper Protein Pair",
                            xlabel="CAI Coevolution",
                            ylabel=expression(rho[C]),
                            num.boot=1000)

func <- plot_grid(p4,p3,labels=c("A","B"),align="v",nrow=2,ncol=1)
ggsave2("../Images/Figure4.pdf",width=150,height =190 ,units="mm")

rate <- plot_grid(p5,p6,labels=c("A","B"),align="v",nrow=2,ncol=1)
ggsave2("../Images/Figure5.pdf",width=120,height =190 ,units="mm")

pdf("../Images/metrics_compared_to_evol_cor_filtered_for_BM.pdf",width=7.5,height = 6.2)
p1
p2
p3
p4
p5
p6
dev.off()
#treat.erc[,"ERC"] <- atanh(treat.erc$ERC)
p7 <- geneExpressionToCoevolutionDistanceMetric(treat.erc, apply.arctanh = T,
                                                main="Effect of Gene Expression on Agreement\nbetween CAI and Empirical-based Measures\nof Coevolution",
                                                xlabel="Mean Ancestral Gene Expression Estimate (Log scale)",
                                                ylabel = expression("(CAI coevolution -"~rho[C]*")"^2))

pdf("../Images/mean_anc_ge_compared_to_distance_between_cai_and_empirical_measures_of_coevolution.pdf",width=7.5,height = 6.2)
p7
dev.off()


## Code for plotting the results of simulated data
# sim.results <- read.table("../Results/Simulated_Results/BM_sim_redo.tsv",sep="\t",header=T,stringsAsFactors=F)
# merged.1 <- merge(df,sim.results,by=c("Prot.1","Prot.2"),suffix=c("_Real",""))
# merged.2 <- merge(df,sim.results,by.x=c("Prot.1","Prot.2"),by.y=c("Prot.2","Prot.1"),suffix=c("_Real",""))
# merged <- rbind(merged.1,merged.2)
# print("Evaluating Fraser et. al. 2004")
# evaluateFraser2004(merged)
# print("Evaluating Martin and Fraser 2018")
# evaluateMartin2018(merged)
# print("Evaluating Cope et. al. 2018")
# evaluatePCM(merged)

# sim.bind <-compareRealSim(df,sim.results,metric="Pearson",category = "binding",main="Pearson Correlation Estimation Performance:\nCorrelated Evolution",xlabel="Parameter Estimated from Real Data",ylabel="Parameter Estimated from Simulated Data")
# sim.control <-compareRealSim(df,sim.results,metric="Pearson",category = "None",main="Pearson Correlation Estimation Performance:\nIndependent Evolution",xlabel="Pearson Correlation (From Simulated Data)",ylabel="Relative Density")

# pdf("../Images/real_sim_compare_pearson.pdf")
# sim.bind
# sim.control
# dev.off()

# sim.bind <-compareRealSim(df,sim.results,metric="Evol.cor",category = "binding",main="Evolutionary Correlation Estimation Performance:\nCorrelated Evolution",xlabel="Parameter Estimated from Real Data",ylabel="Parameter Estimated from Simulated Data")
# sim.control <-compareRealSim(df,sim.results,metric="Evol.cor",category = "None",main="Evolutionary Correlation Estimation Performance:\nIndependent Evolution",xlabel="Pearson Correlation (From Simulated Data)",ylabel="Relative Density")

# pdf("../Images/real_sim_compare_evol_cor.pdf")
# sim.bind
# sim.control
# dev.off()

# ## Code for performing resampling done to control for protein membership
# treat <- df[which(df$Linked=="binding"),]
# control <- df[which(df$Linked=="None"),]

# mean.bind.evol <- numeric(length = 200)
# mean.null.evol <- numeric(length=200)
# p.val.bind.evol <- numeric(length=200)
# p.val.null.evol <- numeric(length=200)
# mean.bind.pear <- numeric(length = 200)
# mean.null.pear <- numeric(length=200)
# p.val.bind.pear <- numeric(length=200)
# p.val.null.pear <- numeric(length=200)
# for (i in 1:200)
# {
#   print(i)
#   bind.unique <- getUniqueMembers(treat,200)
#   null.unique <- getUniqueMembers(control,200)
#   mean.bind.evol[i] <- mean(bind.unique$Evol.cor)
#   mean.null.evol[i] <- mean(null.unique$Evol.cor)
#   p.val.bind.evol[i] <- t.test(bind.unique$Evol.cor)$p.value
#   p.val.null.evol[i] <- t.test(null.unique$Evol.cor)$p.value
#   mean.bind.pear[i] <- mean(bind.unique$Pearson)
#   mean.null.pear[i] <- mean(null.unique$Pearson)
#   p.val.bind.pear[i] <- t.test(bind.unique$Pearson)$p.value
#   p.val.null.pear[i] <- t.test(null.unique$Pearson)$p.value
# }
# #
# x <- round(mean(mean.bind.evol),2)
# y <- round(mean(mean.null.evol),2)
# z <- round(mean(mean.bind.pear),2)
# w <- round(mean(mean.null.pear),2)
# resampled <- data.frame(P=c(mean.bind.evol,mean.null.evol,mean.bind.pear,mean.null.pear),
#                  Metric=c(rep("Evolutionary",400),rep("Pearson",400)),
#                  Group=c(rep(paste0("Evolutionary + Binding, ",x),200),rep(paste0("Evolutionary + Control, ",y),200),rep(paste0("Pearson + Binding, ",z),200),rep(paste0("Pearson + Control, ",w),200)))

# p<-(ggplot(resampled,aes(x=P,fill=Group))
#     + geom_histogram(binwidth=0.01,alpha=0.5,color="black",position="identity",size=1.5)
#     + xlab("Mean Correlation")
#     +scale_fill_viridis(discrete = T)
#     +ylim(c(0,35))
#     + labs(fill=expression("Correlation Metric + Group,"~bar(rho)))
#     + theme(legend.position = c(0.2,0.85))
#     + ggtitle("Mean Correlation Distributions:\n Constraining Protein Membership") + theme(plot.title = element_text(hjust=0.5,size=14,face="bold"),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour="black"),legend.text=element_text(size=10,face="bold"),axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_text(size=14,face="bold"),axis.text=element_text(size=12,face="bold"))
# ) 
# pdf("../Images/resampling_procedure.pdf")
# p
# dev.off()

