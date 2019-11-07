library(ape)
library(mvMORPH)
library(geiger)
library(phytools)
library(phangorn)
library(parallel)

CleanData <- function(phy, data) {
  cleaned <- treedata(phy,data,warnings=T)
  return(cleaned)
}

## TO DO: change to argparser for command line interface
args <- commandArgs(trailingOnly = TRUE)

if(length(args) >= 6){
  tree.file <- args[1]
  gene.exp <- args[2]
  std.exp <- args[3]
  interaction.file <- args[4]
  model.type <- args[5]
  out <- args[6]
  to.drop <- c()
  if (length(args)>6)
  {
    to.drop <- c(args[7:length(args)])
  }
} else
{
  tree.file <- "../Data/fungi_tree.tre"
  gene.exp <- "../Data/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv"
  std.exp <- "../Data/std_err_tpm_matrix_updated_Standard.LogNorm.tsv"
  interaction.file <- "../Data/interactions.txt"
  model.type <- "BM"
  out <- "model_results.tsv"
  to.drop <- c()
}


phylogeneticAnalysis <- function(prot.1,prot.2,exp.df,sem.df,link,score,tree)
{
  tryCatch({
    exp.1 <- exp.df[prot.1,]
    exp.2 <- exp.df[prot.2,]
    sem.1 <- t(sem.df[prot.1,])
    sem.2 <- t(sem.df[prot.2,])
    sem <- data.frame(Sem.1=sem.1,Sem.2=sem.2)
    sem[,1] <- sem[,1]^2
    sem[,2] <- sem[,2]^2
    sem <- sem[tree$tip.label,]
    tmp <- data.frame(Prot.1=prot.1,Prot.2=prot.2,Alpha.cor.1=0,Alpha.cor.1.Error=0,Alpha.cor.2=0,Alpha.cor.2.Error=0,Alpha.cor.shared=0,Alpha.cor.shared.Error=0,Sigma.cor.1=0,Sigma.cor.1.Error=0,Sigma.cor.2=0,Sigma.cor.2.Error=0,Sigma.cov=0,Sigma.cov.Error=0,Alpha.ind.1=0,Alpha.ind.1.Error=0,Alpha.ind.2=0,Alpha.ind.2.Error=0,Sigma.ind.1=0,Sigma.ind.1.Error=0,Sigma.ind.2=0,Sigma.ind.2.Error=0,Likelihood.cor=0,Likelihood.ind=0,LRT.score=0,LRT.pval=0,AICc.cor=0,AICc.ind=0,Conv.Cor=0,Conv.Ind=0,Reliable.corr=0,Reliable.ind=0,Pearson=0,Corr.significant=0,Linked=link,Theta.cor.1=0,Theta.cor.2=0,Theta.ind.1=0,Theta.ind.2=0,Evol.cor=0,Num.missing.1=0,Num.missing.2=0,HalfLife.cor.1=0,HalfLife.cor.2=0,HalfLife.ind.1=0,HalfLife.ind.2=0,Score=score)
    df <- data.frame(Exp.1=exp.1,Exp.2=exp.2)
    df <- df[tree$tip.label,]
    num.obs.1 <- sum(is.na(df[,1]))
    num.obs.2 <- sum(is.na(df[,2]))
    sem <- as.matrix(sem[tree$tip.label,])
      if (model.type != "BM")
      {
        null <- mvOU(tree,df,error = sem,model="OU1",optimization="Nelder-Mead",param=list(decomp="diagonal",decompSigma="diagonal",root=F,vcv="randomRoot"),echo=F,diagnostic = F)
        if (model.type == "OU")
        {
          correlated <- mvOU(tree,df,error = sem,model="OU1",optimization="Nelder-Mead",param=list(root=F,vcv="randomRoot"),echo=F,diagnostic = F) ## Assume entire tree is under same selective regime...later work will allow this to vary
          hl.cor <- halflife(correlated)
          ## note these may not correpsond to the half-life of protein 1 and 2...return in ascending order 
          tmp[,"HalfLife.cor.1"] <- hl.cor[1]
          tmp[,"HalfLife.cor.2"] <- hl.cor[2]
          cor.err <- sqrt(diag(solve(correlated$param$opt$hessian)))
          null.err <- sqrt(diag(solve(null$param$opt$hessian)))
          tmp[,"Alpha.cor.1.Error"] <- cor.err[1]
          tmp[,"Alpha.cor.shared.Error"] <- cor.err[2]
          tmp[,"Alpha.cor.2.Error"] <- cor.err[3]
          tmp[,"Sigma.cor.1.Error"] <- cor.err[4]
          tmp[,"Sigma.cov.Error"] <- cor.err[5]
          tmp[,"Sigma.cor.2.Error"] <- cor.err[6]
        } else if (model.type == "OU_alpha_not_shared")
        {
          correlated <- mvOU(tree,df,error=sem,model="OU1",scale.height=F,optimization="Nelder-Mead",param=list(decomp="diagonal",root=F,vcv="randomRoot"),echo=F,diagnostic = F) ## Assume entire tree is under same selective regime...later work will allow this to vary
          tmp[,"HalfLife.cor.1"] <- log(2)/correlated$alpha[1,1]
          tmp[,"HalfLife.cor.2"] <- log(2)/correlated$alpha[2,2]
          cor.err <- sqrt(diag(solve(correlated$param$opt$hessian)))
          null.err <- sqrt(diag(solve(null$param$opt$hessian)))
          tmp[,"Alpha.cor.1.Error"] <- cor.err[1]
          tmp[,"Alpha.cor.2.Error"] <- cor.err[2]
          tmp[,"Sigma.cor.1.Error"] <- cor.err[3]
          tmp[,"Sigma.cov.Error"] <- cor.err[4]
          tmp[,"Sigma.cor.2.Error"] <- cor.err[5]
        }
        tmp[,"Alpha.ind.1.Error"] <- null.err[1]
        tmp[,"Alpha.ind.2.Error"] <- null.err[2]
        tmp[,"Sigma.ind.1.Error"] <- null.err[3]
        tmp[,"Sigma.ind.2.Error"] <- null.err[4]
        tmp[,"Alpha.cor.1"] <- correlated$alpha[1,1]
        tmp[,"Alpha.cor.2"] <- correlated$alpha[2,2]
        tmp[,"Alpha.cor.shared"] <- correlated$alpha[1,2]
        tmp[,"Alpha.ind.1"] <- null$alpha[1,1]
        tmp[,"Alpha.ind.2"] <- null$alpha[2,2]
        stat.cov <- stationary(correlated)
        evol.cor <- cov2cor(stat.cov)
        tmp[,"HalfLife.ind.1"] <- log(2)/null$alpha[1,1]
        tmp[,"HalfLife.ind.2"] <- log(2)/null$alpha[2,2]
      } else if (model.type == "BM"){
        null <- mvBM(tree,df,error = sem,model="BM1",scale.height=F,optimization="Nelder-Mead",param=list(constraint="diagonal"),echo=F,diagnostic = F)
        correlated <- mvBM(tree,df,error = sem,model="BM1",scale.height=F,optimization="Nelder-Mead",echo=F,diagnostic = F) ## Assume entire tree is under same selective regime...later work will allow this to vary
        evol.cor <- cov2cor(correlated$sigma)
        cor.err <- sqrt(diag(solve(correlated$param$opt$hessian)))
        null.err <- sqrt(diag(solve(null$param$opt$hessian)))
        tmp[,"Sigma.cor.1.Error"] <- cor.err[1]
        tmp[,"Sigma.cov.Error"] <- cor.err[2]
        tmp[,"Sigma.cor.2.Error"] <- cor.err[3]
        tmp[,"Sigma.ind.1.Error"] <- null.err[1]
        tmp[,"Sigma.ind.2.Error"] <- null.err[2]
      }
      lr <- LRT(correlated,null,echo = F)
      corr.test <- cor.test(unlist(exp.df[prot.1,]),unlist(exp.df[prot.2,]),na.action="na.omit")
      tmp[,"Likelihood.cor"] <- correlated$LogLik
      tmp[,"Likelihood.ind"] <- null$LogLik
      tmp[,"Evol.cor"] <- evol.cor[1,2]
      tmp[,"Sigma.cor.1"] <- correlated$sigma[1,1]
      tmp[,"Sigma.cor.2"] <- correlated$sigma[2,2]
      tmp[,"Sigma.cov"] <- correlated$sigma[1,2]
      tmp[,"Sigma.ind.1"] <- null$sigma[1,1]
      tmp[,"Sigma.ind.2"] <- null$sigma[2,2]
      tmp[,"LRT.score"] <- lr$ratio
      tmp[,"LRT.pval"] <- lr$pval
      tmp[,"AICc.cor"] <- correlated$AICc
      tmp[,"AICc.ind"] <- null$AICc
      tmp[,"Conv.Cor"] <- correlated$convergence
      tmp[,"Conv.Ind"] <- null$convergence
      tmp[,"Reliable.corr"] <- correlated$hess.values
      tmp[,"Reliable.ind"] <- null$hess.values
      tmp[,"Pearson"] <- corr.test$estimate
      tmp[,"Corr.significant"] <- corr.test$p.value
      tmp[,"Theta.cor.1"] <- correlated$theta[1]
      tmp[,"Theta.cor.2"] <-correlated$theta[2]
      tmp[,"Theta.ind.1"] <- null$theta[1]
      tmp[,"Theta.ind.2"] <- null$theta[2]
      tmp[,"Num.missing.1"] <- num.obs.1
      tmp[,"Num.missing.2"] <- num.obs.2
      return(tmp)
  },warning= function(war){NULL},error=function(err){print("oops");return(NA)})
}
tree <- read.tree(tree.file)
if (length(to.drop) != 0)
{
  tree <- drop.tip(tree,to.drop)
}
exp.df<- t(read.table(gene.exp,sep="\t",header=T,row.names = 1,stringsAsFactors = F))
cleaned <- CleanData(phy = tree,exp.df)
tree <- cleaned$phy
exp.df <- t(cleaned$data)
interactions <- read.table(interaction.file,sep="\t",header=T,stringsAsFactors = F)
sem.df<- read.table(std.exp,sep="\t",header=T,stringsAsFactors = F,row.names=1)
nb_cores <- 40
results.list <- mclapply(1:nrow(interactions),function(i){parallelAnalysis(prot.1=interactions[i,1],prot.2=interactions[i,2],exp.df=exp.df,sem.df=sem.df,link=interactions[i,3],score=interactions[i,4],tree=tree)},mc.cores = getOption("mc.cores", nb_cores))
results.list <- results.list[!is.na(results.list)]
results <- do.call("rbind",results.list)
write.table(results,out,sep="\t",quote=F,row.names = F)