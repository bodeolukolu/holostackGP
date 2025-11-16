#!/usr/bin/env Rscript

holostackGP <- function(
    wdir = "./",
    projname="proj_GP",
    phenofile="Traits.txt",
    genofile="geno_dose.txt",
    metagenomefile="metagenome.txt",
    ploidy=2,
    traits=c("trait-1","trait-2"),
    covariate=NULL,             # c("trait-9","traits-10")
    kernel=c("genomic"),       # "metagenomic or microbiome", genomic","holobiont", "metagenomic+genomic or microbiome+genomic"
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    gwas_pred=FALSE,
    subsample_markers=NULL,
    topHits=100,
    gene_model=c("Full"),           # "Additive", "Dominance", "metagenome or microbiome", "Full"
    R_libpath=NULL
) {
  load_packages <- function(pkgs) {
    for (p in pkgs) {
      if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }

  # List of required packages
  pkgs <- c("mice", "data.table", "dplyr", "AGHmatrix", "vegan", "compositions",
            "ggcorrplot", "nlme", "lsmeans", "agricolae", "doParallel", "foreach",
            "parallel", "rrBLUP", "BGLR", "GWASpoly")

  # Call it once at the start of your function
  load_packages(pkgs)


  options(error = function() {
    traceback(2)
    quit(status = 1)
  })
  options(warn = 1)
  if (!is.null(R_libpath) && nzchar(R_libpath)) {.libPaths(R_libpath)}

  # gene_model <- match.arg(gene_model, c("Full","Additive","Dominance","metagenomic","microbiome"))


  #############################################################################################################################################################################
  # Specify parameters
  ####################
  # Step 1: Specify covariate if required and import data.
  setwd(wdir)
  GP_run_title <- projname
  traits <- traits
  if (covariate == "" || covariate == "None" || covariate == "NULL") {
    covariate <- NULL
  } else {
    covariate <- strsplit(covariate, ",")[[1]]
  }
  myY <- read.table(phenofile, head = TRUE, sep="\t")
  myG <- read.table(genofile, head = FALSE, sep="\t")
  if(kernel == "metagenomic" || kernel == "microbiome") {kernel <- "gBLUP"}
  if(kernel == "genomic") {kernel <- "GBLUP"}
  if(kernel == "holobiont" || kernel == "microbiome+genomic" || kernel == "metagenomic+genomic") {kernel <- "gGBLUP"}
  gp_model <- kernel
  gwas_Gpred <- gwas_pred
  gwas_model <- "MLM"
  Additional_models <- c("rrBLUP","BRR","BayesA","BayesB","BayesC","BayesianLasso","RHKS")     # "rrBLUP","BRR","BayesA","BayesB","BayesC","BayesianLasso" or NULL

  number_reps <- as.numeric(CVrep)
  nfold_CV <- as.numeric(k_fold)
  maf_threshold <- as.numeric(maf)
  perc_missing <- as.numeric(geno_missing_rate)
  if (is.null(subsample_markers) || subsample_markers %in% c("", "None", "NULL")) {
    subsample_markers <- NULL
  } else {
    subsample_markers <- strsplit(subsample_markers, ",")[[1]]
  }
  if (!is.null(subsample_markers)) {subsample_markers <- as.numeric(subsample_markers)}
  gene_model <- match.arg(gene_model, c("Full","Additive","Dominance","metagenomic","microbiome"))
  select_gwasGPmodel <- NULL                                       #c("2-dom-ref","3-dom-alt","3-dom-ref")
  ploidy_levels <- as.numeric(ploidy)

  weight_by <- "scores"   # use effects, scores (i.e. -log10(pvalues) or pvalues.
  ntop_hits <- 100
  mtraits <- NULL
  ncores <- 5
  LOCO <- FALSE

  #metagenome-based parameters
  mincorr <- 0                                 # minimum correlation coefficient (trait vs subsets of metagenome)
  maxcorr <- 0.8                                 # maximum correlation coefficient (trait vs subsets of metagenome)
  pvalue <- 1
  clr_transform <- TRUE
  pheno_zero_inflated <- FALSE
  metag_zero_inflated <- TRUE
  impute_zeroinflation_metagcov <- FALSE
  corr_coeff <- NULL
  metagenome_data <- metagenomefile
  metag_method <- "Aitchison"
  entire_metagenome <- TRUE
  crosstrait <- NULL
  options(warn=0)



  for (trait in c(traits)) {
    for (ploidy in c(ploidy_levels)) {
      for (gene_model in c(gene_models)) {
        model_selection <- TRUE
        if (gene_model == "full" || gene_model == "FULL" ){ gene_model <- "Full"}
        if (gene_model == "all" || gene_model == "ALL" ){ gene_model <- "All"}
        if (gene_model == "additive" || gene_model == "ADDITIVE" ){ gene_model <- "Additive"}
        if (gene_model == "dominance" || gene_model == "DOMINANCE" ){ gene_model <- "Dominance"}
        if(gp_model == "GBLUP"){metagenome_covariate <- FALSE; metagenome_data <- NULL}
        if(gp_model == "gGBLUP"){metagenome_covariate <- TRUE}
        if(gp_model == "gBLUP"){gene_models <- "metagenome"; metag_method <- "Aitchison"; metagenome_covariate <- FALSE}
        if(is.null(covariate)){multitraitGP <- FALSE} else{multitraitGP <- TRUE}
        metag_pca <- TRUE
        dir.create(GP_run_title)
        colnames(myY)[1] <- "Taxa"
        myY[,2:ncol(myY)] <- sapply(myY[,2:ncol(myY)],as.numeric)
        if (!is.null(covariate)) {
          dY <- subset(myY, select=c("Taxa",trait,covariate))
          if (pheno_zero_inflated == FALSE){dY[][dY[] == "0"] <- NA}
        } else {
          colnames(myY)[1] <- "Taxa"
          dY <- subset(myY, select=c("Taxa",trait))
          if (pheno_zero_inflated == FALSE){dY[][dY[] == "0"] <- NA}
        }
        dY <- na.omit(dY)
        dY <- aggregate(dY[,2],by=list(name=dY$Taxa),data=dY,FUN=mean)
        colnames(dY)[1] <- "Taxa"; colnames(dY)[2] <- trait
        dG <- myG; colnames(dG) <- dG[1,]; dG <- dG[-1,]
        if("pvalue" %in% colnames(dG)){dG <- subset(dG, select=-c(pvalue))}

        if (!is.null(metagenome_data)) {
          metagenome_data <- sub("\\..*", "", metagenome_data)
          metag <- read.table(paste(metagenome_data,".txt",sep=""), header=F, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
          metag[1,] <- sub("_mean", "", metag[1,])
          colnames(metag) <- metag[1,]
          metag <- metag[,-c((ncol(metag)-6):(ncol(metag)))]
          metag <- metag[,c(ncol(metag),1:ncol(metag)-1)]
          metag <- as.data.frame(t(metag))
          colnames(metag) <- metag[1,]; metag <- metag[-1,]
          colnames(metag)[1] <- "Taxa"
          metag <- metag[metag$Taxa %in% dY$Taxa,]
          for (i in c(2:ncol(metag))) { metag[,i] <- as.numeric(metag[,i])}
          metag <- metag[,-1]
          metag$percent <- (rowSums(metag[,1:ncol(metag)] > "0")/ncol(metag))*100
          metag <- subset(metag, percent >= 0)
          metag <- subset(metag, select=-c(percent))
          metag <- t(metag)
          metag <- as.data.frame(metag)
          metag$percent <- (rowSums(metag[,1:ncol(metag)] > "0")/ncol(metag))*100
          metag <- subset(metag, percent >= 5)
          metag <- subset(metag, select=-c(percent))
          metag <- as.data.frame(t(metag))
          metag$Plant_ID <- c(rownames(metag))
          metag <- metag[,c(which(colnames(metag)=="Plant_ID"),which(colnames(metag)!="Plant_ID"))]
          colnames(metag)[1] <- "Taxa"
          keep_taxa <- intersect(dY$Taxa,metag$Taxa)
          metag <- metag[metag$Taxa %in% keep_taxa,]
          dY <- dY[dY$Taxa %in% keep_taxa,]
          row.names(metag) <- metag[,1]; metag <- metag[,-1]
        }
        if (!is.null(metagenome_data) && metag_method != "Aitchison" && length(intersect(colnames(metag),trait)) == "0"){corr_coeff <- "full"}
        if (!is.null(metagenome_data) && metag_method != "Aitchison" && corr_coeff == "full") {
          myYsub <- subset(dY, select=c("Taxa",trait))
          row.names(myYsub) <- c(myYsub[,1]); myYsub <- subset(myYsub, select=-c(1))
          if (!exists("metag_corrsub")){
            metag_corrsub <- merge(myYsub, metag, by=0)
            row.names(metag_corrsub) <- c(metag_corrsub[,1]); metag_corrsub <- subset(metag_corrsub, select=-c(1))
          }
          if (impute_zeroinflation_metagcov == TRUE && !exists("metag_corrsub2")){
            hold_trait <- subset(metag_corrsub, select=c(1))
            metag_corrsub2 <- subset(metag_corrsub, select=-c(1))
            metag_corrsub2[,2:ncol(metag_corrsub2)][metag_corrsub2[,2:ncol(metag_corrsub2)] == "0"] <- NA
            metag_corrsub2 <- metag_corrsub2[rowSums(is.na(metag_corrsub2)) != ncol(metag_corrsub2), ]
            PCpred <- pca(metag_corrsub2, nPcs=1, method="nlpca")
            metag_corrsub2 <- as.data.frame(completeObs(PCpred)); PCpred <- NULL
            metag_corrsub2 <- merge(hold_trait, metag_corrsub2, by = 'row.names', all = TRUE)
          } else {metag_corrsub2 <- metag_corrsub}
          tcorr <- cor(metag_corrsub2[sapply(metag_corrsub2, is.numeric)], method="spearman", use="pairwise.complete.obs")
          tpmat <- cor_pmat(metag_corrsub2[sapply(metag_corrsub2, is.numeric)], method = "spearman", exact=FALSE)
          tcorr <- subset(tcorr[,1], abs(tcorr[,1]) <= maxcorr  & abs(tcorr[,1]) >= mincorr )
          tcorr <- as.data.frame(tcorr)
          tpmat <- subset(tpmat[,1], tpmat[,1] <= pvalue)
          tpmat <- as.data.frame(tpmat)
          tcorr <- merge (tcorr,tpmat,by="row.names")
          print(tcorr)
          if (nrow(tcorr) > 0){
            metag_sub <- metag[, (colnames(metag) %in% tcorr[,1])]
            if (min(metag_sub[], na.rm=T) < 0){
              metag_sub[] <- metag_sub[] + abs(min(metag_sub[], na.rm=T))
            }
          }
          if (nrow(tcorr) > 0){
            metag_sub[is.na(metag_sub)] <- 0
            covariates <- prcomp(metag_sub, center = T,scale= F)
            var_explained <- (covariates$sdev^2/sum(covariates$sdev^2))*100
            var_explained[1:5]
            covariates <- as.data.frame(covariates$x)
            if (min(covariates[!is.na(covariates)])-1 < 0){
              for (i in c(1:ncol(covariates))) {
                covariates[,i] <- covariates[,i] + abs(min(covariates[!is.na(covariates)])-1)
              }
            }
            for (i in c(1:ncol(covariates))) {
              covariates[,i] <- covariates[,i]/max(covariates[], na.rm=TRUE)
            }
            if (ncol(covariates) > 2) {
              covariates <- subset(covariates, select=c(1:3))
              covariates$PC4 <- covariates$PC1 + covariates$PC2 + covariates$PC3
              covariates[] <- covariates[]/(max(covariates$PC4)+(max(covariates$PC4)*0.1))
              covariates <- covariates[,-c(4)]
              covariates$Taxa <- row.names(covariates)
              covariates[is.na(covariates)] <- 0
              covariates <- covariates[,c(4,1,2,3)]
            } else {
              if (ncol(covariates) == 2) {
                covariates$PC3 <- covariates$PC1 + covariates$PC2
                covariates[] <- covariates[]/(max(covariates$PC3)+(max(covariates$PC3)*0.1))
                covariates <- covariates[,-c(3)]
                covariates$Taxa <- row.names(covariates)
                covariates <- covariates[,c(3,1,2)]
              }
              if (ncol(covariates) == 1) {
                covariates$PC2 <- covariates$PC1
                covariates[] <- covariates[]/(max(covariates$PC2)+(max(covariates$PC2)*0.1))
                covariates <- subset(covariates, select=-c(2))
                covariates$Taxa <- row.names(covariates)
                covariates <- covariates[,c(2,1)]
              }
            }
          } else {covariates <- NULL}
        }
        if (!is.null(metagenome_data) && metag_method != "Aitchison" && corr_coeff != "full" && file.exists(paste(corr_coeff,sep=""))){
          tcorr <- read.table(paste(corr_coeff,sep=""), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
          tcorr <- subset(tcorr, tcorr[,1]==trait)
          tcorr <- subset(tcorr, abs(tcorr[,3]) <= maxcorr  & abs(tcorr[,3]) >= mincorr & tcorr[,4] <= pvalue)
          rownames(tcorr) <- tcorr[,2]; tcorr <- subset(tcorr, select=c(3,4)); colnames(tcorr)[1:2] <- c("tcorr","tpmat")
          tcorr <- setdiff(rownames(tcorr), trait)

          if (length(tcorr) > 0){
            metag_sub <- metag[,c(intersect(colnames(metag),tcorr))]
            if (min(metag_sub[], na.rm=T) < 0){
              metag_sub[] <- metag_sub[] + abs(min(metag_sub[], na.rm=T))
            }
          }
          if (length(tcorr) > 0){
            if(metag_pca == TRUE){
              metag_sub[is.na(metag_sub)] <- 0
              covariates <- prcomp(metag_sub, center = T,scale= F)
              var_explained <- (covariates$sdev^2/sum(covariates$sdev^2))*100
              var_explained[1:5]
              covariates <- as.data.frame(covariates$x)
              if (min(covariates[!is.na(covariates)])-1 < 0){
                for (i in c(1:ncol(covariates))) {
                  covariates[,i] <- covariates[,i] + abs(min(covariates[!is.na(covariates)])-1)
                }
              }
              for (i in c(1:ncol(covariates))) {
                covariates[,i] <- covariates[,i]/max(covariates[], na.rm=TRUE)
              }
              if (ncol(covariates) > 2) {
                covariates <- subset(covariates, select=c(1:3))
                covariates$PC4 <- covariates$PC1 + covariates$PC2 + covariates$PC3
                covariates[] <- covariates[]/(max(covariates$PC4)+(max(covariates$PC4)*0.1))
                covariates <- covariates[,-c(4)]
                covariates$Taxa <- row.names(covariates)
                covariates[is.na(covariates)] <- 0
                covariates <- covariates[,c(4,1,2,3)]
              } else {
                if (ncol(covariates) == 2) {
                  covariates$PC3 <- covariates$PC1 + covariates$PC2
                  covariates[] <- covariates[]/(max(covariates$PC3)+(max(covariates$PC3)*0.1))
                  covariates <- covariates[,-c(3)]
                  covariates$Taxa <- row.names(covariates)
                  covariates <- covariates[,c(3,1,2)]
                }
                if (ncol(covariates) == 1) {
                  covariates$PC2 <- covariates$PC1
                  covariates[] <- covariates[]/(max(covariates$PC2)+(max(covariates$PC2)*0.1))
                  covariates <- subset(covariates, select=-c(2))
                  covariates$Taxa <- row.names(covariates)
                  covariates <- covariates[,c(2,1)]
                }
              }
            } else {
              metag_sub[is.na(metag_sub)] <- 0
              covariates <- metag_sub
              covariates$Taxa <- row.names(covariates)
              covariates <- covariates[,c(which(colnames(covariates)=="Taxa"),which(colnames(covariates)!="Taxa"))]
            }
          } else {covariates <- NULL}
        }
        if (metagenome_covariate == TRUE && metag_method != "Aitchison") {
          if (!is.null(covariate)) {
            myCV <- merge(covariate,covariates, by="Taxa")
          } else {
            myCV <- covariates
          }
        }
        if (!is.null(covariate)) {myCV <- subset(myY, select=c("Taxa",covariate))} else {myCV <- NULL}
        if (metagenome_covariate == TRUE && metag_method != "Aitchison") {
          if (!is.null(covariate)) {
            myCV <- merge(covariate,covariates, by="Taxa")
          } else {
            myCV <- covariates
          }
        }
        if (gwas_Gpred == TRUE) {
          if (metagenome_covariate == TRUE) {dir.create(paste(trait,"_",gene_model,ploidy,"x","_gwas_",gene_model,"_metag_",sep="")); setwd(paste(trait,"_",gene_model,ploidy,"x","_gwas_",gene_model,"_metag_",sep=""))
          }else{dir.create(paste(trait,"_",gene_model,ploidy,"x","_gwas_",gene_model,"_nometag_",sep="")); setwd(paste(trait,"_",gene_model,ploidy,"x","_gwas_",gene_model,"_nometag_",sep=""))}
        }
        if (gwas_Gpred != TRUE) {
          if (metagenome_covariate == TRUE) {dir.create(paste(trait,"_",gene_model,ploidy,"x","_",gp_model,"_metag_",sep="")); setwd(paste(trait,"_",gene_model,ploidy,"x","_",gp_model,"_metag_",sep=""))
          } else {dir.create(paste(trait,"_",gene_model,ploidy,"x","_",gp_model,"_nometag_",sep=""));
            setwd(paste(trait,"_",gene_model,ploidy,"x","_",gp_model,"_nometag_",sep=""))}
        }
        normalize_kinmat <- function(kinmat){
          #normalize kinship so that Kij \in [0,1]
          tmp=kinmat - min(kinmat)
          tmp=tmp/max(tmp)
          tmp[1:9,1:9]
          #fix eigenvalues to positive
          diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
          tmp[1:9,1:9]
          return(tmp)
        }
        outdir <- gsub(".*/","",getwd())
        final_result_summary <- data.frame(matrix(nrow = 0, ncol = 4))
        if (file.exists(paste("../",GP_run_title,"/",outdir,"gp_summary_stats.txt",sep=""))){
          stop("Output file(s) alreadfy exist → exiting script.")
        }
        if (!file.exists(paste("../",GP_run_title,"/",outdir,"gp_summary_stats.txt",sep=""))){
          gp_summary_stats <- data.frame(matrix(nrow = 0, ncol = 4))
          colnames(gp_summary_stats) <- c("Traits","Median","Mean","StdErr")
          write.table(gp_summary_stats, paste("../",GP_run_title,"/",outdir,"gp_summary_stats.txt",sep=""), quote = FALSE, sep = "\t")
        }


        #Step 2: Run Genomic Prediction
        #######################################################################################
        #Initial
        #######################################################################################
        if (nrow(final_result_summary) == 0) {
          t=number_reps #total replicates
          s=1/nfold_CV #sample of inference, e.g. set it to 1/5 for five fold cross validation
          Y.raw=dY[,c(1,2)]#choose a trait
          Y.raw=Y.raw[!is.na(Y.raw[,2]),] #Remove missing data

          if (metagenome_covariate == TRUE || gwas_Gpred == TRUE) {
            RMIP <- data.frame(matrix(ncol = 4, nrow = 0))
            colnames(RMIP) <- c("SNP","CHROM","POS","Rep")
            write.table(RMIP, paste(GP_run_title,"_RMIP.txt",sep=""), col.names=FALSE, row.names=FALSE, quote = FALSE, sep = "\t")
          }
          if (gwas_Gpred == TRUE) {
            if (metagenome_covariate == TRUE) {
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_gwas_",gene_model,"_metag.txt",sep=""),append=TRUE)
            } else {
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_gwas_",gene_model,".txt",sep=""),append=TRUE)
            }
          }
          if (gwas_Gpred == FALSE) {
            if (metagenome_covariate == TRUE) {
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_",gene_model,"_metag.txt",sep=""),append=TRUE)
            } else {
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_",gene_model,".txt",sep=""),append=TRUE)
            }
          }

          top_percentile <- 0.95; low_percentile <- 0.05
          topP <- subset(Y.raw, Y.raw[,2] >= quantile(Y.raw[,2],top_percentile)); topP$class <- "high_value"
          lowP <- subset(Y.raw, Y.raw[,2] <= quantile(Y.raw[,2],low_percentile)); lowP$class <- "low_value"
          Precision_in <- NULL; Precision_out <- data.frame(matrix(nrow = 0, ncol = 6))
          colnames(Precision_out) <- c("Trait","Rep","Precision_topP","topP_n","Precision_lowP","lowP_n")

          if (gp_model != "GBLUP" && metag_method == "Aitchison") {
            if(!is.null(corr_coeff) && corr_coeff != "full"){
              corr_coeff <- sub("\\..*", "", corr_coeff)
              tcorr <- read.table(paste("../",corr_coeff,".txt",sep=""), header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
              tcorr <- subset(tcorr, tcorr[,1]==trait)
              tcorr <- subset(tcorr, abs(tcorr[,3]) <= maxcorr  & abs(tcorr[,3]) >= mincorr & tcorr[,4] <= pvalue)
              rownames(tcorr) <- tcorr[,2]; tcorr <- subset(tcorr, select=c(3,4)); colnames(tcorr)[1:2] <- c("tcorr","tpmat")
              tcorr <- setdiff(rownames(tcorr), trait)
            }
            if(entire_metagenome == FALSE){
              if(trait %in% colnames(metag)){metag <- metag[,c(intersect(colnames(metag),c(colnames(metag)[1],tcorr)))]}
            }
            if(trait %in% colnames(metag)){if(entire_metagenome == TRUE){metag <- metag[ , -which(names(metag) %in% c(trait))]}}
            metag$Taxa <- row.names(metag)
            metag <- metag[,c(ncol(metag),1:(ncol(metag)-1))]
            metag <- metag[metag$Taxa %in% dY$Taxa,]
            for (i in c(2:ncol(metag))) { metag[,i] <- as.numeric(metag[,i])}
            metag <- metag[,-1]
            metag$percent <- (rowSums(metag[,1:ncol(metag)] > "0")/ncol(metag))*100
            metag <- subset(metag, percent >= 0)
            metag <- subset(metag, select=-c(percent))
            metag <- t(metag)
            metag <- as.data.frame(metag)
            metag$percent <- (rowSums(metag[,1:ncol(metag)] > "0")/ncol(metag))*100
            metag <- subset(metag, percent >= 0)
            metag <- subset(metag, select=-c(percent))
            metag <- as.data.frame(t(metag))
            metag$Plant_ID <- c(rownames(metag))
            metag <- metag[,c(which(colnames(metag)=="Plant_ID"),which(colnames(metag)!="Plant_ID"))]
            colnames(metag)[1] <- "Taxa"
            keep_taxa <- intersect(dY$Taxa,metag$Taxa)
            metag <- metag[metag$Taxa %in% keep_taxa,]
            dY <- dY[dY$Taxa %in% keep_taxa,]
            row.names(metag) <- metag[,1]; metag <- metag[,-1]
            if(clr_transform == TRUE){
              metag_clr <- as.matrix(clr(metag + 1e-6))
            } else {metag_clr <- metag}
            # metag_clrx <- t(metag_clr)
            # write.table(metag_clrx, paste("../clr_",metagenome_data,".txt",sep=""), col.names=TRUE, row.names=T, quote = FALSE, sep = "\t", append=FALSE)
            # dist_metag <- dist(metag_clr)
            metagKI <- tcrossprod(scale(metag_clr)) / ncol(metag_clr)
            metagKI <- normalize_kinmat(as.matrix(metagKI))
            metagKI <- metagKI
            metagKI <- metagKI[Y.raw[,1],Y.raw[,1]]
            metagKI <- cbind(rownames(metagKI), data.frame(metagKI, row.names=NULL));
            colnames(metagKI) <- c(1:ncol(metagKI))
            metagKI <- as.data.frame(metagKI)
            metagKI[,2:ncol(metagKI)] <- sapply(metagKI[,2:ncol(metagKI)], as.numeric)
            # metagKIx <- metagKI
            # write.table(metagKIx, paste("../metagKI_",metagenome_data,".txt",sep=""), col.names=TRUE, row.names=T, quote = FALSE, sep = "\t", append=FALSE)
            #
            # myKI <- as.data.frame(myKI)
            # myKI$Plant_ID <- rownames(myKI)
            # myKI <- myKI[,c(which(colnames(myKI)=="Plant_ID"),which(colnames(myKI)!="Plant_ID"))]
            # rownames(myKI) <- c(1:nrow(myKI)); colnames(myKI) <- c(1:ncol(myKI))
          }
          if (gene_model != "metagenome") {
            ktaxa <- intersect(Y.raw[,1],colnames(dG))
            ktaxa <- setdiff(ktaxa,drop_ind)
            Y.raw <- Y.raw[is.element(Y.raw$Taxa, ktaxa),]
            ktaxa <- c("SNP","CHROM","POS","REF","ALT",ktaxa)
            dG <- as.matrix(dG[,which(colnames(dG) %in% ktaxa)])
            dG <- as.data.frame(dG)
            if (ploidy == 2){
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*2 + (rowSums(dG == "1", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
            }
            if (ploidy == 4){
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*4 + (rowSums(dG == "1", na.rm = TRUE))*3 +
                (rowSums(dG == "2", na.rm = TRUE))*2 + (rowSums(dG == "3", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2 +
                (rowSums(dG == "3", na.rm = TRUE))*3 + (rowSums(dG == "4", na.rm = TRUE))*4
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
            }
            if (ploidy == 6){
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*6 + (rowSums(dG == "1", na.rm = TRUE))*5 +
                (rowSums(dG == "2", na.rm = TRUE))*4 + (rowSums(dG == "3", na.rm = TRUE))*3 +
                (rowSums(dG == "4", na.rm = TRUE))*2 + (rowSums(dG == "5", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2 +
                (rowSums(dG == "3", na.rm = TRUE))*3 + (rowSums(dG == "4", na.rm = TRUE))*4 +
                (rowSums(dG == "5", na.rm = TRUE))*5 + (rowSums(dG == "6", na.rm = TRUE))*6
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
            }
            if (ploidy == 8){
              dG$freq0 <- (rowSums(dG == "0", na.rm = TRUE))*8 + (rowSums(dG == "1", na.rm = TRUE))*7 +
                (rowSums(dG == "2", na.rm = TRUE))*6 + (rowSums(dG == "3", na.rm = TRUE))*5 +
                (rowSums(dG == "4", na.rm = TRUE))*4 + (rowSums(dG == "5", na.rm = TRUE))*3 +
                (rowSums(dG == "6", na.rm = TRUE))*2 + (rowSums(dG == "7", na.rm = TRUE))*1
              dG$freq1 <- (rowSums(dG == "1", na.rm = TRUE))*1 + (rowSums(dG == "2", na.rm = TRUE))*2 +
                (rowSums(dG == "3", na.rm = TRUE))*3 + (rowSums(dG == "4", na.rm = TRUE))*4 +
                (rowSums(dG == "5", na.rm = TRUE))*5 + (rowSums(dG == "6", na.rm = TRUE))*6 +
                (rowSums(dG == "7", na.rm = TRUE))*7 + (rowSums(dG == "8", na.rm = TRUE))*8
              maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
              dG$min <- apply(dG[,(ncol(dG)-1):ncol(dG)], 1, function(x)x[maxn(2)(x)])
              dG$sum <- rowSums(dG[,c("freq0","freq1")], na.rm=TRUE)
              dG$maf <- as.numeric(dG$min)/as.numeric(dG$sum)
            }

            dG<- subset(dG, maf > maf_threshold)
            dG <- subset(dG, select=-c(freq0,freq1,min,maf,sum))
            dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)],as.numeric)
            if (ploidy == 2 && max(dG[,6:ncol(dG)],na.rm=TRUE) == 4){
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "2", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "3", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "4", replacement = "2")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)],as.numeric)
            }
            if (ploidy == 2 && max(dG[,6:ncol(dG)],na.rm=TRUE) == 6){
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "2", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "3", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "4", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "5", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "6", replacement = "2")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)],as.numeric)
            }
            if (ploidy == 2 && max(dG[,6:ncol(dG)],na.rm=TRUE) == 8){
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "2", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "3", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "4", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "5", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "6", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "7", replacement = "1")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)], gsub, pattern = "8", replacement = "2")
              dG[,6:ncol(dG)] <- lapply(dG[,6:ncol(dG)],as.numeric)
            }
            pop_data <- subset(dG, select=-c(1:5))
            pop_data$no_missing <- apply(pop_data, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
            pop_data <- subset(pop_data, no_missing < ncol(pop_data)*perc_missing)
            pop_data <- subset(pop_data, select=-c(no_missing))
            set.seed(123)
            if (is.numeric(subsample_markers)) {pop_data <- sample_n(pop_data, subsample_markers)}
            pop_data <- as.matrix(t(pop_data))
            #Computing the full-autopolyploid matrix based on Slater 2016 (Eq. 8 and 9)
            if (gene_model == "Additive"){
              G_matrix <- Gmatrix(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA,
                                  maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                  pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                  ratio.check = FALSE)
              myKI <- normalize_kinmat(as.matrix(G_matrix))
            }
            if (gene_model == "Dominance"){
              if (ploidy == 2){Gmethod <- "Vitezica"}
              if (ploidy > 2){Gmethod <- "Slater"}
              G_matrix <- Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
                                  maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                  pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                  ratio.check = FALSE)
              myKI <- normalize_kinmat(as.matrix(G_matrix))
            }
            if (gene_model == "Full" || gene_model == "All" ){
              G_matrix <- Gmatrix(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA,
                                  maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                  pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                  ratio.check = FALSE)
              myKI.Add <- normalize_kinmat(as.matrix(G_matrix))
              if (ploidy == 2){Gmethod <- "Vitezica"}
              if (ploidy > 2){Gmethod <- "Slater"}
              G_matrix <- Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
                                  maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                  pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                  ratio.check = FALSE)
              myKI.Dom <- normalize_kinmat(as.matrix(G_matrix))
            }

            if (gene_model == "QTL_only") {myKI <- NULL}
            if (gwas_Gpred == TRUE){
              predG <- dG
              if (ploidy == 2 && max(predG[,6:ncol(predG)],na.rm=TRUE) == 4){
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "2", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "3", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "4", replacement = "2")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)],as.numeric)
              }
              if (ploidy == 2 && max(predG[,6:ncol(predG)],na.rm=TRUE) == 6){
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "2", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "3", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "4", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "5", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "6", replacement = "2")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)],as.numeric)
              }
              if (ploidy == 2 && max(predG[,6:ncol(predG)],na.rm=TRUE) == 8){
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "2", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "3", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "4", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "5", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "6", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "7", replacement = "1")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "8", replacement = "2")
                predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)],as.numeric)
              }
              colnames(predG)[1:3] <- c("rs","chrom","pos")
              predGx <- predG[,1:5]; predG <- predG[,-c(2:5)]
              predGx$alleles <- paste(predGx$REF,"/",predGx$ALT,sep="")
              predGx$chrom <- gsub("Chr0","",predGx$chrom)
              predGx$chrom <- gsub("Chr","",predGx$chrom)
              predGx <- subset(predGx, select=-c(4:5))
              predGx <- predGx[,c(1,4,2,3)]
              predGx$strand <- NA; predGx$assembly <- NA; predGx$center <- NA
              predGx$protLSID <- NA; predGx$assayLSID <- NA; predGx$panel <- NA
              predGx$QCcode <- NA
              predG <- merge(predGx,predG,by="rs")
              predG <- rbind(colnames(predG),predG); colnames(predG) <- c(1:ncol(predG))
            }
          }
          if (!is.null(covariate)){
            Y.cov <- myY[,c("Taxa",covariate)]
            Y.cov <- merge(Y.raw, Y.cov, by="Taxa")
            imp <- mice(Y.cov, m = 5, method = "pmm", seed = 123)
            Y.cov <- complete(imp, 1)   # first imputed dataset
            Y.cov <- Y.cov[,c("Taxa",covariate)]
          }

          n=nrow(Y.raw)
          n.missing=round(n*s)
          storage.GBLUP=matrix(NA,t,1)
          storage.rrBLUP=matrix(NA,t,1)
          storage.RKHS=matrix(NA,t,1)
          storage.BRR=matrix(NA,t,1)
          storage.BayesA=matrix(NA,t,1)
          storage.BayesB=matrix(NA,t,1)
          storage.BayesC=matrix(NA,t,1)
          storage.BayesL=matrix(NA,t,1)
          storage.mStacked=matrix(NA,t,1)
          storagec.GBLUP=matrix(NA,t,1)
          storagec.rrBLUP=matrix(NA,t,1)
          storagec.RKHS=matrix(NA,t,1)
          storagec.BRR=matrix(NA,t,1)
          storagec.BayesA=matrix(NA,t,1)
          storagec.BayesB=matrix(NA,t,1)
          storagec.BayesC=matrix(NA,t,1)
          storagec.BayesL=matrix(NA,t,1)
          storagec.mStacked=matrix(NA,t,1)


          #Loop on replicates
          for(rep in seq(from=1, to=t, by=1)){
            #Set missing data
            sample.missing=sample(1:n,n.missing)
            if(n.missing>0){
              Y0=Y.raw[-sample.missing,]
              test_ids <- setdiff(Y.raw$Taxa, Y0$Taxa)
              Y.masked <- Y.raw; rownames(Y.masked) <- Y.raw$Taxa
              Y.masked[rownames(Y.masked) %in% test_ids, 2] <- NA
              Y.masked <- subset(Y.masked, select=-c(1))
              rownames(Y.raw) <- Y.raw[,1]
              #
              if (!is.null(covariate)){
                Y0.cov=Y.cov[-sample.missing,]
                Y.covmasked <- Y.cov; rownames(Y.covmasked) <- Y.cov$Taxa
                for (covids in 2:ncol(Y.covmasked)){Y.covmasked[rownames(Y.covmasked) %in% test_ids, covids] <- NA}
                Y.covmasked <- subset(Y.covmasked, select=-c(1))
                rownames(Y.cov) <- Y.cov[,1]
              }
            }else{
              Y0=Y.raw
              Y.masked <- Y.raw
              rownames(Y.raw) <- Y.raw[,1]
              if (!is.null(covariate)){
                Y0.cov=Y.cov
                Y.covmasked <- Y.cov
                rownames(Y.cov) <- Y.cov[,1]
              }
            }

            # GWAS marker select
            if (gwas_Gpred == TRUE) {
              dG0 <- dG
              ktaxa <- c("SNP","CHROM","POS","REF","ALT",Y0$Taxa)
              dG0 <- as.data.frame(dG0[,ktaxa])
              dG0 <- dG0[]
              if (ploidy == 2){
                dG0$freq0 <- (rowSums(dG0 == "0", na.rm = TRUE))*2 + (rowSums(dG0 == "1", na.rm = TRUE))*1
                dG0$freq1 <- (rowSums(dG0 == "1", na.rm = TRUE))*1 + (rowSums(dG0 == "2", na.rm = TRUE))*2
                maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
                dG0$min <- apply(dG0[,(ncol(dG0)-1):ncol(dG0)], 1, function(x)x[maxn(2)(x)])
                dG0$sum <- rowSums(dG0[,c("freq0","freq1")], na.rm=TRUE)
                dG0$maf <- as.numeric(dG0$min)/as.numeric(dG0$sum)
              }
              if (ploidy == 4){
                dG0$freq0 <- (rowSums(dG0 == "0", na.rm = TRUE))*4 + (rowSums(dG0 == "1", na.rm = TRUE))*3 +
                  (rowSums(dG0 == "2", na.rm = TRUE))*2 + (rowSums(dG0 == "3", na.rm = TRUE))*1
                dG0$freq1 <- (rowSums(dG0 == "1", na.rm = TRUE))*1 + (rowSums(dG0 == "2", na.rm = TRUE))*2 +
                  (rowSums(dG0 == "3", na.rm = TRUE))*3 + (rowSums(dG0 == "4", na.rm = TRUE))*4
                maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
                dG0$min <- apply(dG0[,(ncol(dG0)-1):ncol(dG0)], 1, function(x)x[maxn(2)(x)])
                dG0$sum <- rowSums(dG0[,c("freq0","freq1")], na.rm=TRUE)
                dG0$maf <- as.numeric(dG0$min)/as.numeric(dG0$sum)
              }
              if (ploidy == 6){
                dG0$freq0 <- (rowSums(dG0 == "0", na.rm = TRUE))*6 + (rowSums(dG0 == "1", na.rm = TRUE))*5 +
                  (rowSums(dG0 == "2", na.rm = TRUE))*4 + (rowSums(dG0 == "3", na.rm = TRUE))*3 +
                  (rowSums(dG0 == "4", na.rm = TRUE))*2 + (rowSums(dG0 == "5", na.rm = TRUE))*1
                dG0$freq1 <- (rowSums(dG0 == "1", na.rm = TRUE))*1 + (rowSums(dG0 == "2", na.rm = TRUE))*2 +
                  (rowSums(dG0 == "3", na.rm = TRUE))*3 + (rowSums(dG0 == "4", na.rm = TRUE))*4 +
                  (rowSums(dG0 == "5", na.rm = TRUE))*5 + (rowSums(dG0 == "6", na.rm = TRUE))*6
                maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
                dG0$min <- apply(dG0[,(ncol(dG0)-1):ncol(dG0)], 1, function(x)x[maxn(2)(x)])
                dG0$sum <- rowSums(dG0[,c("freq0","freq1")], na.rm=TRUE)
                dG0$maf <- as.numeric(dG0$min)/as.numeric(dG0$sum)
              }
              if (ploidy == 8){
                dG0$freq0 <- (rowSums(dG0 == "0", na.rm = TRUE))*8 + (rowSums(dG0 == "1", na.rm = TRUE))*7 +
                  (rowSums(dG0 == "2", na.rm = TRUE))*6 + (rowSums(dG0 == "3", na.rm = TRUE))*5 +
                  (rowSums(dG0 == "4", na.rm = TRUE))*4 + (rowSums(dG0 == "5", na.rm = TRUE))*3 +
                  (rowSums(dG0 == "6", na.rm = TRUE))*2 + (rowSums(dG0 == "7", na.rm = TRUE))*1
                dG0$freq1 <- (rowSums(dG0 == "1", na.rm = TRUE))*1 + (rowSums(dG0 == "2", na.rm = TRUE))*2 +
                  (rowSums(dG0 == "3", na.rm = TRUE))*3 + (rowSums(dG0 == "4", na.rm = TRUE))*4 +
                  (rowSums(dG0 == "5", na.rm = TRUE))*5 + (rowSums(dG0 == "6", na.rm = TRUE))*6 +
                  (rowSums(dG0 == "7", na.rm = TRUE))*7 + (rowSums(dG0 == "8", na.rm = TRUE))*8
                maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
                dG0$min <- apply(dG0[,(ncol(dG0)-1):ncol(dG0)], 1, function(x)x[maxn(2)(x)])
                dG0$sum <- rowSums(dG0[,c("freq0","freq1")], na.rm=TRUE)
                dG0$maf <- as.numeric(dG0$min)/as.numeric(dG0$sum)
              }
              dG0<- subset(dG0, maf > maf_threshold)
              dG0 <- subset(dG0, select=-c(freq0,freq1,min,maf,sum))
              dG0[,6:ncol(dG0)] <- lapply(dG0[,6:ncol(dG0)], as.numeric)
              dG0$no_missing <- apply(dG0, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
              dG0 <- subset(dG0, no_missing < (ncol(dG0)-5)*perc_missing)
              dG0 <- subset(dG0, select=-c(no_missing))
              pop_data <- subset(dG0, select=-c(1:5))
              pop_data <- as.matrix(t(pop_data))
              if (gene_model == "Additive"){
                G_matrix <- Gmatrix(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA,
                                    maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                    pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                    ratio.check = FALSE)
                CVkinship <- normalize_kinmat(as.matrix(G_matrix))
              }
              if (gene_model == "Additive_Dominance"){
                G_matrix <- Gmatrix(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA,
                                    maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                    pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                    ratio.check = FALSE)
                CVkinship.Add <- normalize_kinmat(as.matrix(G_matrix))
                if (ploidy == 2){Gmethod <- "Vitezica"}
                if (ploidy > 2){Gmethod <- "Slater"}
                G_matrix <- Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
                                    maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                    pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                    ratio.check = FALSE)
                CVkinship.Dom <- normalize_kinmat(as.matrix(G_matrix))
                CVkinship<- CVkinship.Add * CVkinship.Dom
              }
              if (gene_model == "Dominance"){
                if (ploidy == 2){Gmethod <- "Vitezica"}
                if (ploidy > 2){Gmethod <- "Slater"}
                G_matrix <- Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
                                    maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                    pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                    ratio.check = FALSE)
                CVkinship<- normalize_kinmat(as.matrix(G_matrix))
              }
              if (gene_model == "Additive" || gene_model == "Additive_Dominance" || gene_model == "Dominance"){
                CVkinship<- as.data.frame(CVkinship)
                CVkinship$Plant_ID <- rownames(CVkinship)
                CVkinship <- CVkinship[,c(which(colnames(CVkinship)=="Plant_ID"),which(colnames(CVkinship)!="Plant_ID"))]
                rownames(CVkinship) <- c(1:nrow(CVkinship)); colnames(CVkinship) <- c(1:ncol(CVkinship))
              }
              if(!is.null(metagenome_data) && !is.null(myCV)){phenosamp <- merge(Y0,myCV, by="Taxa")}
              write.table(Y0, "phenosamp.csv", row.names=F, quote = FALSE, sep = ",")
              write.table(dG0, "geno.csv", row.names=F, col.names=T, quote = FALSE, sep = ",")
              row.names(CVkinship) <- CVkinship[,1]; CVkinship <- CVkinship[,-1]
              colnames(CVkinship) <- row.names(CVkinship)
              data <- read.GWASpoly(ploidy=ploidy, pheno.file="phenosamp.csv", geno.file="geno.csv", format="numeric", n.traits=1, delim=",")
              Kinship <- set.K(data, LOCO = FALSE, K=as.matrix(CVkinship))
              if(LOCO == TRUE){Kinship <- set.K(data, LOCO = TRUE,n.core=ncores)}
              cores <- ncores
              #perform gwas
              if (gene_model == "Additive"){
                if (ncol(Y0) == 2) {
                  if ( ploidy == "2" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = ncores, quiet = FALSE)
                    models <- c("additive")
                  }
                  if ( ploidy == "4" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = ncores, quiet = FALSE)
                    models <- c("additive")
                  }
                  if ( ploidy == "6" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = ncores, quiet = FALSE)
                    models <- c("additive")
                  }
                  if ( ploidy == "8" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"),
                                            traits=c(traitname), n.core = ncores, quiet = FALSE)
                    models <- c("additive")
                  }
                } else {
                  if(ncol(Y0) == 3){
                    params <- set.params(fixed=c("PC1"), fixed.type=rep("numeric",1),n.PC=1,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }
                  if(ncol(Yo) == 4){
                    params <- set.params(fixed=c("PC1","PC2"), fixed.type=rep("numeric",2),n.PC=2,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }
                  if(ncol(Y0) > 4){
                    params <- set.params(fixed=c("PC1","PC2","PC3"), fixed.type=rep("numeric",3),n.PC=3,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }

                  if ( ploidy == "2" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("additive")
                  }
                  if ( ploidy == "4" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("additive")
                  }
                  if ( ploidy == "6" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("additive")
                  }
                  if ( ploidy == "8" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive"),
                                            traits=c(traitname), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("additive")
                  }
                }
              }
              if (gene_model == "Dominance"){
                if (ncol(Y0) == 2) {
                  if ( ploidy == "2" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom"), n.core = ncores, quiet = FALSE)
                    models <- c("1-dom-ref","1-dom-alt")
                  }
                  if ( ploidy == "4" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom"), n.core = ncores, quiet = FALSE)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt")
                  }
                  if ( ploidy == "6" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","3-dom"), n.core = ncores, quiet = FALSE)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt")
                  }
                  if ( ploidy == "8" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom", "3-dom", "4-dom"),
                                            traits=c(traitname), n.core = ncores, quiet = FALSE)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","4-dom-ref","4-dom-alt")
                  }
                } else {
                  if(ncol(Y0) == 3){
                    params <- set.params(fixed=c("PC1"), fixed.type=rep("numeric",1),n.PC=1,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }
                  if(ncol(Yo) == 4){
                    params <- set.params(fixed=c("PC1","PC2"), fixed.type=rep("numeric",2),n.PC=2,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }
                  if(ncol(Y0) > 4){
                    params <- set.params(fixed=c("PC1","PC2","PC3"), fixed.type=rep("numeric",3),n.PC=3,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }

                  if ( ploidy == "2" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt")
                  }
                  if ( ploidy == "4" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt")
                  }
                  if ( ploidy == "6" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","3-dom"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","additive")
                  }
                  if ( ploidy == "8" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom", "3-dom", "4-dom"),
                                            traits=c(traitname), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","4-dom-ref","4-dom-alt")
                  }
                }
              }
              if (gene_model == "Additive_Dominance"){
                if (ncol(Y0) == 2) {
                  if ( ploidy == "2" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive","1-dom"), n.core = ncores, quiet = FALSE)
                    models <- c("additive","1-dom-ref","1-dom-alt")
                  }
                  if ( ploidy == "4" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive","1-dom","2-dom"), n.core = ncores, quiet = FALSE)
                    models <- c("additive","1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt")
                  }
                  if ( ploidy == "6" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive","1-dom","2-dom","3-dom"), n.core = ncores, quiet = FALSE)
                    models <- c("additive","1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt")
                  }
                  if ( ploidy == "8" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("additive","1-dom","2-dom", "3-dom", "4-dom"),
                                            traits=c(traitname), n.core = ncores, quiet = FALSE)
                    models <- c("additive","1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","4-dom-ref","4-dom-alt")
                  }
                } else {
                  if(ncol(Y0) == 3){
                    params <- set.params(fixed=c("PC1"), fixed.type=rep("numeric",1),n.PC=1,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }
                  if(ncol(Yo) == 4){
                    params <- set.params(fixed=c("PC1","PC2"), fixed.type=rep("numeric",2),n.PC=2,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }
                  if(ncol(Y0) > 4){
                    params <- set.params(fixed=c("PC1","PC2","PC3"), fixed.type=rep("numeric",3),n.PC=3,MAF=0.001,geno.freq=0.99,P3D=TRUE)
                  }

                  if ( ploidy == "2" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","additive"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt")
                  }
                  if ( ploidy == "4" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","additive"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","additive")
                  }
                  if ( ploidy == "6" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom","3-dom","additive"), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","additive")
                  }
                  if ( ploidy == "8" ) {
                    GWAS.fitted <- GWASpoly(Kinship, models = c("1-dom","2-dom", "3-dom", "4-dom","additive"),
                                            traits=c(traitname), n.core = ncores, quiet = FALSE, params=params)
                    models <- c("1-dom-ref","1-dom-alt","2-dom-ref","2-dom-alt","3-dom-ref","3-dom-alt","4-dom-ref","4-dom-alt","additive")
                  }
                }
              }

              #select top hits if additive model
              if (gene_model == "Additive"){
                if(weight_by == "effects") {
                  GWAS_effects <- as.data.frame(GWAS.fitted@effects)
                  GWAS_effects$abs <- abs(GWAS_effects[,1])
                  GWAS_effects <- na.omit(GWAS_effects)
                  GWAS_effects <- GWAS_effects[order(GWAS_effects$abs, decreasing = TRUE),]
                  GWAS_effects <- GWAS_effects[1:(ntop_hits*3),]
                  GWAS_effects$SNP <- rownames(GWAS_effects)
                  wdG <- dG[dG$SNP %in% c(GWAS_effects$SNP),]
                  wdGlist <- c(GWAS_effects$SNP)
                }
                if(weight_by == "scores") {
                  GWAS_scores <- as.data.frame(GWAS.fitted@scores)
                  GWAS_scores <- na.omit(GWAS_scores)
                  GWAS_scores[,1] <- GWAS_scores[order(GWAS_scores[,1], decreasing = TRUE),]
                  GWAS_scores$SNP <- rownames(GWAS_scores)
                  GWAS_scores <- GWAS_scores[1:(ntop_hits*3),]
                  wdG <- dG[dG$SNP %in% c(GWAS_scores$SNP),]
                  wdGlist <- c(GWAS_scores$SNP)
                }
                if(weight_by == "pvalues") {
                  GWAS_scores <- as.data.frame(GWAS.fitted@scores)
                  GWAS_pvalue<- GWAS_scores;  GWAS_pvalue[,1] <- 10^-(GWAS_pvalue[,1])
                  GWAS_pvalue <- na.omit(GWAS_pvalue)
                  GWAS_pvalue[,1] <- GWAS_pvalue[order(GWAS_pvalue[,1], decreasing = FALSE),]
                  GWAS_pvalue$SNP <- rownames(GWAS_pvalue)
                  GWAS_pvalue <- GWAS_pvalue[1:(ntop_hits*3),]
                  wdG <- dG[dG$SNP %in% c(GWAS_pvalue$SNP),]
                  wdGlist <- c(GWAS_pvalue$SNP)
                }
              }
              if (gene_model == "Dominance" || gene_model == "Additive_Dominance"){
                if(weight_by == "effects") {
                  GWAS_effects <- as.data.frame(GWAS.fitted@effects)
                  if(!is.null(select_gwasGPmodel)){
                    select_gwasGPmodel <- gsub("-", ".", select_gwasGPmodel)
                    GWAS_effects <- GWAS_effects %>% select(contains(select_gwasGPmodel))
                  }
                  GWAS_effects_sum <- subset(GWAS_effects, select=c(1))
                  GWAS_effects_sum <- na.omit(GWAS_effects_sum)
                  GWAS_effects_sum$abs <- abs(GWAS_effects_sum[,1])
                  GWAS_effects_sum <- GWAS_effects_sum[order(GWAS_effects_sum$abs, decreasing = TRUE),]
                  GWAS_effects_sum$SNP <- rownames(GWAS_effects_sum)
                  GWAS_effects_sum <- GWAS_effects_sum[1:(ntop_hits*3),]
                  colnames(GWAS_effects_sum) <- c("values","abs","SNP")
                  for (selm in c(2:ncol(GWAS_effects))){
                    GWAS_effects_sumH <- subset(GWAS_effects, select=c(selm))
                    GWAS_effects_sumH <- na.omit(GWAS_effects_sumH)
                    GWAS_effects_sumH$abs <- abs(GWAS_effects_sumH[,1])
                    GWAS_effects_sumH <- GWAS_effects_sumH[order(GWAS_effects_sumH$abs, decreasing = TRUE),]
                    GWAS_effects_sumH$SNP <- rownames(GWAS_effects_sumH)
                    GWAS_effects_sumH <- GWAS_effects_sumH[1:(ntop_hits*3),]
                    colnames(GWAS_effects_sumH) <- c("values","abs","SNP")
                    GWAS_effects_sum <- rbind(GWAS_effects_sum, GWAS_effects_sumH)
                    GWAS_effects_sumH <- NULL
                  }
                  GWAS_effects_sum <- GWAS_effects_sum[order(GWAS_effects_sum$values, decreasing = TRUE),]
                  GWAS_effects_sum <- GWAS_effects_sum[!duplicated(GWAS_effects_sum$SNP),]
                  wdG <- dG[dG$SNP %in% c(GWAS_effects_sum$SNP[1:(ntop_hits*3)]),]
                  wdGlist <- c(GWAS_effects_sum$SNP[1:(ntop_hits*3)])
                }
                if(weight_by == "scores") {
                  GWAS_scores <- as.data.frame(GWAS.fitted@scores)
                  if(!is.null(select_gwasGPmodel)){
                    select_gwasGPmodel <- gsub("-", ".", select_gwasGPmodel)
                    GWAS_scores <- GWAS_scores %>% select(contains(select_gwasGPmodel))
                  }
                  GWAS_scores_sum <- subset(GWAS_scores, select=c(1))
                  GWAS_scores_sum <- na.omit(GWAS_scores_sum)
                  GWAS_scores_sum[,1] <- GWAS_scores_sum[order(GWAS_scores_sum[,1], decreasing = TRUE),]
                  GWAS_scores_sum$SNP <- rownames(GWAS_scores_sum)
                  GWAS_scores_sum <- GWAS_scores_sum[1:(ntop_hits*3),]
                  colnames(GWAS_scores_sum) <- c("values","SNP")
                  for (selm in c(2:ncol(GWAS_scores))){
                    GWAS_scores_sumH <- subset(GWAS_scores, select=c(selm))
                    GWAS_scores_sumH <- na.omit(GWAS_scores_sumH)
                    GWAS_scores_sumH[,1] <- GWAS_scores_sumH[order(GWAS_scores_sumH[,1], decreasing = TRUE),]
                    GWAS_scores_sumH$SNP <- rownames(GWAS_scores_sumH)
                    GWAS_scores_sumH <- GWAS_scores_sumH[1:(ntop_hits*3),]
                    colnames(GWAS_scores_sumH) <- c("values","SNP")
                    GWAS_scores_sum <- rbind(GWAS_scores_sum, GWAS_scores_sumH)
                    GWAS_scores_sumH <- NULL
                  }
                  GWAS_scores_sum <- GWAS_scores_sum[order(GWAS_scores_sum$values, decreasing = TRUE),]
                  GWAS_scores_sum <- GWAS_scores_sum[!duplicated(GWAS_scores_sum$SNP),]
                  wdG <- dG[dG$SNP %in% c(GWAS_scores_sum$SNP[1:(ntop_hits*3)]),]
                  wdGlist <- c(GWAS_scores_sum$SNP[1:(ntop_hits*3)])
                }
                if(weight_by == "pvalues") {
                  GWAS_scores <- as.data.frame(GWAS.fitted@scores); GWAS_pvalue <- GWAS_scores
                  if(!is.null(select_gwasGPmodel)){
                    select_gwasGPmodel <- gsub("-", ".", select_gwasGPmodel)
                    GWAS_pvalue <- GWAS_pvalue %>% select(contains(select_gwasGPmodel))
                  }
                  GWAS_pvalue_sum <- subset(GWAS_pvalue, select=c(1))
                  GWAS_pvalue_sum <- na.omit(GWAS_pvalue_sum)
                  GWAS_pvalue[,1] <- 10^-(GWAS_pvalue[,1])
                  GWAS_pvalue_sum[,1] <- GWAS_pvalue_sum[order(GWAS_pvalue_sum[,1], decreasing = FALSE),]
                  GWAS_pvalue_sum$SNP <- rownames(GWAS_pvalue_sum)
                  GWAS_pvalue_sum <- GWAS_pvalue_sum[1:(ntop_hits*3),]
                  colnames(GWAS_pvalue_sum) <- c("values","SNP")
                  for (selm in c(2:ncol(GWAS_pvalue))){
                    GWAS_pvalue_sumH <- subset(GWAS_pvalue, select=c(selm))
                    GWAS_pvalue_sumH <- na.omit(GWAS_pvalue_sumH)
                    GWAS_pvalue[,1] <- 10^-(GWAS_pvalue[,1])
                    GWAS_pvalue_sumH[,1] <- GWAS_pvalue_sumH[order(GWAS_pvalue_sumH[,1], decreasing = FALSE),]
                    GWAS_pvalue_sumH$SNP <- rownames(GWAS_pvalue_sumH)
                    GWAS_pvalue_sumH <- GWAS_pvalue_sumH[1:(ntop_hits*3),]
                    colnames(GWAS_pvalue_sumH) <- c("values","SNP")
                    GWAS_pvalue_sum <- rbind(GWAS_pvalue_sum, GWAS_pvalue_sumH)
                    GWAS_pvalue_sumH <- NULL
                  }
                  GWAS_scores_pvalue <- GWAS_pvalue_sum[order(GWAS_pvalue_sum$values, decreasing = TRUE),]
                  GWAS_scores_pvalue <- GWAS_pvalue_sum[!duplicated(GWAS_pvalue_sum$SNP),]
                  wdG <- dG[dG$SNP %in% c(GWAS_pvalue_sum$SNP[1:(ntop_hits*3)]),]
                  wdGlist <- c(GWAS_pvalue_sum$SNP[1:(ntop_hits*3)])
                }
              }

              if(exists("mselRMIP")){
                if (model_selection == TRUE){
                  GRM <- as.matrix(CVkinship)
                  fixed.marker <- as.data.frame(t(wdG[,-c(2:5)]))
                  replace_na_with_mode_rowwise <- function(df) {
                    apply(df, 1, function(row) {
                      if (any(is.na(row))) {  # Check if the row has NA values
                        non_na_values <- row[!is.na(row)]  # Get non-NA values
                        if (length(non_na_values) > 0) {
                          mode_value <- names(sort(table(non_na_values), decreasing = TRUE))[1]  # Find mode
                          row[is.na(row)] <- mode_value  # Replace NA with mode
                        }
                      }
                      return(row)
                    }) |> t() |> as.data.frame()
                  }
                  fixed.marker <- replace_na_with_mode_rowwise(fixed.marker)
                  sum(is.na(fixed.marker))
                  colnames(fixed.marker) <- fixed.marker[1,]; fixed.marker <- fixed.marker[-1,]
                  fixed.marker <- fixed.marker[c(Y0[,1]),]
                  fixed.marker <- apply(as.matrix(fixed.marker),2,as.numeric)
                  rownames(fixed.marker) <- Y0[,1]
                  rownames(GRM) <- Y0[,1]
                  rownames(Y0) <- Y0[,1]
                  marker_variance <- apply(fixed.marker, 2, var, na.rm=TRUE)
                  fixed.marker <- fixed.marker[, !is.na(marker_variance) & marker_variance > 0]
                  cor_matrix <- cor(fixed.marker, method="spearman")
                  cor_matrix[lower.tri(cor_matrix)] <- NA; diag(cor_matrix) <- NA
                  threshold <- 0.8
                  highly_correlated <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
                  correlated_pairs <- apply(highly_correlated, 1, function(x) colnames(cor_matrix)[x])
                  to_remove <- NULL
                  if (length(correlated_pairs) > 0){to_remove <- unique(correlated_pairs[2, ])}
                  fixed.marker <- fixed.marker[, !colnames(fixed.marker) %in% to_remove]
                  for (msteps in c(10,9,8,7,6,5,4,3,2)){
                    mlmm_result <- tryCatch({mlmm(Y = Y0[,2], X = fixed.marker, K = GRM,
                                                  nbchunks = 2, maxsteps = msteps, thresh = 1.2 * 10^-5)
                    }, error = function(e) {
                      cat("Error in iteration", msteps, ":", conditionMessage(e), "\n")
                      return(NULL)  # Skip iteration on error
                    })
                    if (!is.null(mlmm_result)) {
                      cat("Success at maxsteps =", msteps, "\n")
                      break  # Stop the loop once a successful run occurs
                    }
                  }
                  if (!is.null(mlmm_result)) {
                    fit <- mlmm_result[["opt_thresh"]][["out"]]
                    fit$abs <- abs(fit[,3])
                    fit <- fit[order(fit[,4], decreasing = TRUE),]
                    fit <- fit[1:ntop_hits,]; fit <- fit[,1]
                    wdG <- dG[dG$SNP %in% c(fit),]
                  }
                  if (is.null(mlmm_result)) {wdG <- dG[dG$SNP %in% c(wdGlist[1:ntop_hits]),]}
                }
              }
              if (model_selection == TRUE){
                GRM <- as.matrix(CVkinship)
                fixed.marker <- as.data.frame(t(wdG[,-c(2:5)]))
                replace_na_with_mode_rowwise <- function(df) {
                  apply(df, 1, function(row) {
                    if (any(is.na(row))) {  # Check if the row has NA values
                      non_na_values <- row[!is.na(row)]  # Get non-NA values
                      if (length(non_na_values) > 0) {
                        mode_value <- names(sort(table(non_na_values), decreasing = TRUE))[1]  # Find mode
                        row[is.na(row)] <- mode_value  # Replace NA with mode
                      }
                    }
                    return(row)
                  }) |> t() |> as.data.frame()
                }
                fixed.marker <- replace_na_with_mode_rowwise(fixed.marker)
                sum(is.na(fixed.marker))
                colnames(fixed.marker) <- fixed.marker[1,]; fixed.marker <- fixed.marker[-1,]
                fixed.marker <- fixed.marker[c(Y0[,1]),]
                fixed.marker <- apply(as.matrix(fixed.marker),2,as.numeric)
                rownames(fixed.marker) <- Y0[,1]
                rownames(GRM) <- Y0[,1]
                rownames(Y0) <- Y0[,1]
                marker_variance <- apply(fixed.marker, 2, var, na.rm=TRUE)
                fixed.marker <- fixed.marker[, !is.na(marker_variance) & marker_variance > 0]
                cor_matrix <- cor(fixed.marker, method="spearman")
                cor_matrix[lower.tri(cor_matrix)] <- NA; diag(cor_matrix) <- NA
                threshold <- 0.8
                highly_correlated <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
                correlated_pairs <- apply(highly_correlated, 1, function(x) colnames(cor_matrix)[x])
                to_remove <- NULL
                if (length(correlated_pairs) > 0){to_remove <- unique(correlated_pairs[2, ])}
                fixed.marker <- fixed.marker[, !colnames(fixed.marker) %in% to_remove]

                num_cores <- ncores
                cl <- makeCluster(num_cores)
                registerDoParallel(cl)
                final_fit <- foreach(repm = 1:100, .combine = 'rbind') %dopar% {
                  n=nrow(Y0)
                  n.missing=round(n*0.1)
                  sample.missing=sample(1:n,n.missing)
                  Y0m <- Y0
                  rownames(Y0m) <- 1:nrow(Y0m)
                  train_indices <- setdiff(1:n,  -sample.missing)
                  X_train <- GRM[-sample.missing,  -sample.missing]
                  y_train <- Y0m[c( -sample.missing),]
                  fixed.marker_train <- fixed.marker
                  fixed.marker_train <- fixed.marker_train[-sample.missing,]
                  for (msteps in c(10,9,8,7,6,5,4,3,2)){
                    mlmm_result <- tryCatch({mlmm(Y=y_train[,2], X=fixed.marker_train, K=X_train,
                                                  nbchunks = 2, maxsteps = msteps, thresh = 1.2 * 10^-5)
                    }, error = function(e) {
                      cat("Error in iteration", msteps, ":", conditionMessage(e), "\n")
                      return(NULL)  # Skip iteration on error
                    })
                    if (!is.null(mlmm_result)) {
                      fit0 <-as.data.frame(mlmm_result[["opt_thresh"]][["out"]])
                      fit0$abs <- abs(fit0[,3])
                      fit0 <- fit0[order(fit0[,4], decreasing = TRUE),]
                      fit0 <- fit0[1:ntop_hits,]
                      fit0 <- na.omit(fit0)
                    } else {
                      fit0 <- data.frame(matrix(ncol = 4, nrow = 0))
                      colnames(fit0) <- c("SNP","pval","effect","abs")
                    }
                    return(fit0)
                  }
                }
                stopCluster(cl)
                if (nrow(final_fit) == 0) {
                  wdG <- dG[dG$SNP %in% wdGlist[1:ntop_hits],]
                } else {
                  fit <- as.data.frame(table(final_fit$SNP))
                  fit <- fit[order(fit[,2], decreasing = TRUE),]
                  fit <- fit[1:ntop_hits,]; fit <- fit[,1]
                  wdG <- dG[dG$SNP %in% c(fit),]
                }
              }

              RMIP_iterate <- as.data.frame(wdG[,1:3])
              names(RMIP_iterate) <- c("SNP","CHROM","POS")
              RMIP_iterate$Rep <- rep
              write.table(RMIP_iterate, paste(GP_run_title,"_RMIP.txt",sep=""), col.names=FALSE, row.names=FALSE, quote = FALSE, sep = "\t", append=TRUE)


              #Recompute the full-autopolyploid matrix based on Slater 2016 (Eq. 8 and 9)
              ktaxa <- c("SNP","CHROM","POS","REF","ALT",Y.raw$Taxa)
              wdG <- as.data.frame(wdG[,ktaxa])
              pop_data <- wdG
              pop_data <- subset(pop_data, select=-c(1:5))
              pop_data$no_missing <- apply(pop_data, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
              pop_data <- subset(pop_data, no_missing < ncol(pop_data)*perc_missing)
              pop_data <- subset(pop_data, select=-c(no_missing))
              if (is.numeric(subsample_markers)) {pop_data <- sample_n(pop_data, subsample_markers)}
              pop_data <- as.matrix(t(pop_data))
              if (gene_model == "Additive"){
                G_matrix <- Gmatrix(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA,
                                    maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                    pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                    ratio.check = FALSE)
                myKI <- normalize_kinmat(as.matrix(G_matrix))
              }
              if (gene_model == "Dominance"){
                if (ploidy == 2){Gmethod <- "Vitezica"}
                if (ploidy > 2){Gmethod <- "Slater"}
                G_matrix <- Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
                                    maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                    pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                    ratio.check = FALSE)
                myKI <- normalize_kinmat(as.matrix(G_matrix))
              }
              if (gene_model == "Additive" || gene_model == "Additive_Dominance" || gene_model == "Dominance"){
                myKI <- as.data.frame(myKI)
                myKI$Plant_ID <- rownames(myKI)
                myKI <- myKI[,c(which(colnames(myKI)=="Plant_ID"),which(colnames(myKI)!="Plant_ID"))]
                rownames(myKI) <- c(1:nrow(myKI)); colnames(myKI) <- c(1:ncol(myKI))
              }
              # if (gene_model == "QTL_only") {myKI <- NULL}
              if (gwas_Gpred == TRUE){
                predG <- dG
                if (ploidy == 2 && max(predG[,6:ncol(predG)],na.rm=TRUE) == 4){
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "2", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "3", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "4", replacement = "2")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)],as.numeric)
                }
                if (ploidy == 2 && max(predG[,6:ncol(predG)],na.rm=TRUE) == 6){
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "2", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "3", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "4", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "5", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "6", replacement = "2")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)],as.numeric)
                }
                if (ploidy == 2 && max(predG[,6:ncol(predG)],na.rm=TRUE) == 8){
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "2", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "3", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "4", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "5", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "6", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "7", replacement = "1")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)], gsub, pattern = "8", replacement = "2")
                  predG[,6:ncol(predG)] <- lapply(predG[,6:ncol(predG)],as.numeric)
                }
                colnames(predG)[1:3] <- c("rs","chrom","pos")
                predGx <- predG[,1:5]; predG <- predG[,-c(2:5)]
                predGx$alleles <- paste(predGx$REF,"/",predGx$ALT,sep="")
                predGx$chrom <- gsub("Chr0","",predGx$chrom)
                predGx$chrom <- gsub("Chr","",predGx$chrom)
                predGx <- subset(predGx, select=-c(4:5))
                predGx <- predGx[,c(1,4,2,3)]
                predGx$strand <- NA; predGx$assembly <- NA; predGx$center <- NA
                predGx$protLSID <- NA; predGx$assayLSID <- NA; predGx$panel <- NA
                predGx$QCcode <- NA
                predG <- merge(predGx,predG,by="rs")
                predG <- rbind(colnames(predG),predG); colnames(predG) <- c(1:ncol(predG))
              }
              if (metagenome_covariate == TRUE) {
                if (!is.null(myCV)) {
                  myCV <- merge(Y.raw,myCV, by="Taxa", all.x=TRUE)
                  myCV <- na.omit(myCV)
                  Y.raw <- subset(myCV, select=c(1,2))
                  myCV <- subset(myCV, select=-c(2))
                }
              }
            }

            #Prediction
            nIter <- 10000           # number of iterations: Bayesian methods using BGLR package
            burnIn <- 2000          # number of burnin iterations: Bayesian methods using BGLR package
            verboseBGLR <- FALSE
            if(is.null(myCV)){
              if (gene_model == "Additive" || gene_model == "Dominance"){
                myCV <- as.data.frame(myKI[,1]); myCV$dummy <- "0.95"
                colnames(myCV)[1] <- "Taxa"; myCV$dummy <- as.numeric(myCV$dummy)
              }
              if (gene_model == "metagenome"){
                myCV <- as.data.frame(metagKI[,1]); myCV$dummy <- "0.95"
                colnames(myCV)[1] <- "Taxa"; myCV$dummy <- as.numeric(myCV$dummy)
              }
              if (gene_model == "Full" || gene_model == "All")    {
                myCV <- as.data.frame(myKI.Add[,1]); myCV$dummy <- "0.95"
                colnames(myCV)[1] <- "Taxa"; myCV$dummy <- as.numeric(myCV$dummy)
              }

            }

            # format genomic/metagenomic relationship matrix
            if(gp_model == "GBLUP" || gp_model == "gGBLUP"){
              if (gene_model == "Full" || gene_model == "All"){
                myKI.A <- t(myKI.Add)
                myKI.A <- as.matrix(myKI.A); myKI.A <- apply(myKI.A, 2, as.numeric); storage.mode(myKI.A) <- "numeric"
                rownames(myKI.A) <- colnames(myKI.A)
                common_ids <- intersect(colnames(myKI.A), rownames(Y.masked))
                myKI.A <- myKI.A[common_ids, common_ids]

                myKI.D <- t(myKI.Dom)
                myKI.D <- as.matrix(myKI.D); myKI.D <- apply(myKI.D, 2, as.numeric); storage.mode(myKI.D) <- "numeric"
                rownames(myKI.D) <- colnames(myKI.D)
                common_ids <- intersect(colnames(myKI.D), rownames(Y.masked))
                myKI.D <- myKI.D[common_ids, common_ids]
              } else {
                myKIx <- t(myKI)
                myKIx <- as.matrix(myKIx); myKIx <- apply(myKIx, 2, as.numeric); storage.mode(myKIx) <- "numeric"
                rownames(myKIx) <- colnames(myKIx)
                common_ids <- intersect(colnames(myKIx), rownames(Y.masked))
                myKIx <- myKIx[common_ids, common_ids]
              }
            }
            if(gp_model == "gBLUP" || gp_model == "gGBLUP"){
              metagKIx <- t(metagKI); colnames(metagKIx) <- metagKIx[1,]; metagKIx <- metagKIx[-1,]; rownames(metagKIx) <- colnames(metagKIx)
              metagKIx <- as.matrix(metagKIx); metagKIx <- apply(metagKIx, 2, as.numeric); storage.mode(metagKIx) <- "numeric"
              rownames(metagKIx) <- colnames(metagKIx)
              common_ids <- intersect(colnames(metagKIx), rownames(Y.masked))
              metagKIx <- metagKIx[common_ids, common_ids]
              make.positive.definite <- function(K, tol = 1e-8) {
                eig <- eigen(K, symmetric = TRUE)
                eig$values[eig$values < tol] <- tol
                K.new <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
                return(K.new)
              }
              K.psd <- make.positive.definite(metagKIx)
              rownames(K.psd) <- rownames(metagKIx); colnames(K.psd) <- rownames(metagKIx)
              metagKIx <- K.psd
            }
            Y.masked <- Y.masked[rownames(Y.masked) %in% common_ids, , drop = FALSE]

            # format data to produce scaled genomic/metagenoic data matrix
            if(gp_model == "GBLUP" || gp_model == "gGBLUP"){
              if (gwas_Gpred == TRUE){
                geno <- t(wdG[, -c(1:5)])
                colnames(geno) <- wdG$SNP
                geno <- geno[rownames(geno) %in% common_ids, ]; keeprownames <- rownames(geno)
                geno <- apply(geno, 2, function(x) as.numeric(as.character(x)))
                rownames(geno) <- keeprownames
                geno <- as.data.frame(geno)
                geno[] <- apply(geno, 2, function(col) {
                  if (anyNA(col)) {
                    col[is.na(col)] <- median(col, na.rm = TRUE)
                  }
                  col
                })
                geno <- as.matrix(geno)  # ensure matrix
                # Scale genotype matrix based on the specified genetic model
                if (gene_model == "Full" || gene_model == "All") {
                  geno_add <- scale(geno, center = TRUE, scale = TRUE)                 # Additive component
                  # geno_dom <- scale(ifelse(geno == 1, 1, 0), center = TRUE, scale = TRUE)  # Dominance component
                  geno_dom <- scale(geno * (ploidy - geno), center = TRUE, scale = TRUE) # model dominance as quadratic deviation from the linear additive dosage
                  geno.A_scaled <- geno_add[, colSums(is.na(geno_add)) == 0]; rownames(geno.A_scaled) <- common_ids
                  geno.D_scaled <- geno_dom[, colSums(is.na(geno_add)) == 0]; rownames(geno.D_scaled) <- common_ids
                } else {
                  if (gene_model == "Additive") {
                    geno_scaled <- scale(geno, center = TRUE, scale = TRUE)
                  }
                  if (gene_model == "Dominance") {
                    # geno_dom <- ifelse(geno == 1, 1, 0)  # Encode dominance (1 for heterozygotes)
                    geno_dom <- scale(geno * (ploidy - geno), center = TRUE, scale = TRUE) # model dominance as quadratic deviation from the linear additive dosage
                    geno_scaled <- scale(geno_dom, center = TRUE, scale = TRUE)
                    geno_scaled <- geno_scaled[, colSums(is.na(geno_scaled)) == 0]
                  }
                  if (gene_model == "Additive_Dominance") {
                    geno_add <- scale(geno, center = TRUE, scale = TRUE)                 # Additive component
                    # geno_dom <- scale(ifelse(geno == 1, 1, 0), center = TRUE, scale = TRUE)  # Dominance component
                    geno_dom <- scale(geno * (ploidy - geno), center = TRUE, scale = TRUE) # model dominance as quadratic deviation from the linear additive dosage
                    geno_scaled <- cbind(geno_add, geno_dom)
                    geno_scaled <- geno_scaled[, colSums(is.na(geno_scaled)) == 0]
                  }
                  rownames(geno_scaled) <- common_ids
                  geno_scaled <- geno_scaled[, colSums(is.na(geno_scaled)) == 0]
                }
              } else {
                # geno <- t(dG[, -c(1:5)])
                # colnames(geno) <- dG$SNP
                geno <- pop_data
                geno <- geno[rownames(geno) %in% common_ids, ]
                geno <- apply(geno, 2, function(x) as.numeric(as.character(x)))
                geno <- as.data.frame(geno)
                geno[] <- apply(geno, 2, function(col) {
                  if (anyNA(col)) {
                    col[is.na(col)] <- median(col, na.rm = TRUE)
                  }
                  col
                })
                geno <- as.matrix(geno)  # ensure matrix
                # Scale genotype matrix based on the specified genetic model
                if (gene_model == "Full" || gene_model == "All") {
                  geno_add <- scale(geno, center = TRUE, scale = TRUE)                 # Additive component
                  # geno_dom <- scale(ifelse(geno == 1, 1, 0), center = TRUE, scale = TRUE)  # Dominance component
                  geno_dom <- scale(geno * (ploidy - geno), center = TRUE, scale = TRUE) # model dominance as quadratic deviation from the linear additive dosage
                  geno.A_scaled <- geno_add[, colSums(is.na(geno_add)) == 0]; rownames(geno.A_scaled) <- common_ids
                  geno.D_scaled <- geno_dom[, colSums(is.na(geno_add)) == 0]; rownames(geno.D_scaled) <- common_ids
                  geno.A_scaled[is.na(geno.A_scaled)] <- 0
                  geno.D_scaled[is.na(geno.D_scaled)] <- 0
                } else {
                  if (gene_model == "Additive") {
                    geno_add <- scale(geno, center = TRUE, scale = TRUE)
                    geno_scaled <- scale(geno_add, center = TRUE, scale = TRUE)
                    geno_scaled <- geno_scaled[, colSums(is.na(geno_scaled)) == 0]
                  }
                  if (gene_model == "Dominance") {
                    # geno_dom <- ifelse(geno == 1, 1, 0)  # Encode dominance (1 for heterozygotes)
                    geno_dom <- scale(geno * (ploidy - geno), center = TRUE, scale = TRUE) # model dominance as quadratic deviation from the linear additive dosage
                    geno_scaled <- scale(geno_dom, center = TRUE, scale = TRUE)
                    geno_scaled <- geno_scaled[, colSums(is.na(geno_scaled)) == 0]
                  }
                  rownames(geno_scaled) <- common_ids
                }
              }
            }
            if(gp_model == "gBLUP" || gp_model == "gGBLUP"){
              mgeno <- metag_clr
              mgeno <- mgeno[rownames(mgeno) %in% common_ids, ]
              mgeno <- apply(mgeno, 2, function(x) as.numeric(as.character(x)))
              mgeno <- as.data.frame(mgeno)
              mgeno[] <- apply(mgeno, 2, function(col) {
                if (anyNA(col)) {
                  col[is.na(col)] <- median(col, na.rm = TRUE)
                }
                col
              })
              mgeno <- as.matrix(mgeno)  # ensure matrix
              mgeno_scaled <- scale(mgeno, center = TRUE, scale = TRUE)
              rownames(mgeno_scaled) <- common_ids
            }
            Y.masked <- as.numeric(Y.masked[[1]])

            # compute epistatic kernels
            if (gene_model == "Full" || gene_model == "All"){
              if(gp_model == "GBLUP"){
                K_A <- myKI.A
                K_D <- myKI.D
                geno.A_scaled <- as.matrix(geno.A_scaled); mode(geno.A_scaled) <- "numeric"
                geno.D_scaled <- as.matrix(geno.D_scaled); mode(geno.D_scaled) <- "numeric"
                # Function to compute Hadamard (element-wise) product for interactions
                interaction_kernel <- function(K1, K2) {(K1 * K2)}
                # Compute interaction kernels
                K_AxD   <- interaction_kernel(K_A, K_D)
                K_AxA   <- interaction_kernel(K_A, K_A)
                K_DxD   <- interaction_kernel(K_D, K_D)
                K_A <- (K_A + t(K_A)) / 2
                K_D <- (K_D + t(K_D)) / 2
                K_AxD <- (K_AxD + t(K_AxD)) / 2
                K_AxA <- (K_AxA + t(K_AxA)) / 2
                K_DxD <- (K_DxD + t(K_DxD)) / 2
                kernels <- list(A = K_A, D = K_D, AxD = K_AxD, AxA = K_AxA, DxD = K_DxD)
              }
              if(gp_model == "gBLUP"){
                K_M <- metagKIx
                kernels <- list(M=K_M)
              }
              if(gp_model == "gGBLUP"){
                K_M <- metagKIx
                K_A <- myKI.A
                K_D <- myKI.D
                geno.A_scaled <- as.matrix(geno.A_scaled); mode(geno.A_scaled) <- "numeric"
                geno.D_scaled <- as.matrix(geno.D_scaled); mode(geno.D_scaled) <- "numeric"
                # Function to compute Hadamard (element-wise) product for interactions
                interaction_kernel <- function(K1, K2) {(K1 * K2)}
                # Compute interaction kernels
                K_AxD   <- interaction_kernel(K_A, K_D)
                K_AxA   <- interaction_kernel(K_A, K_A)
                K_DxD   <- interaction_kernel(K_D, K_D)
                K_A <- (K_A + t(K_A)) / 2
                K_D <- (K_D + t(K_D)) / 2
                K_AxD <- (K_AxD + t(K_AxD)) / 2
                K_AxA <- (K_AxA + t(K_AxA)) / 2
                K_DxD <- (K_DxD + t(K_DxD)) / 2
                mgeno_scale <- as.matrix(mgeno_scaled); mode(mgeno_scaled) <- "numeric"
                K_AxAxM <- interaction_kernel(K_AxA, K_M)
                K_DxDxM <- interaction_kernel(K_DxD, K_M)
                K_AxDxM <- interaction_kernel(K_AxD, K_M)
                K_AxAxM <- (K_AxAxM + t(K_AxAxM)) / 2
                K_DxDxM <- (K_DxDxM + t(K_DxDxM)) / 2
                K_AxDxM <- (K_AxDxM + t(K_AxDxM)) / 2
                K_M <- (K_M + t(K_M)) / 2
                kernels <- list(A = K_A, D = K_D, AxD = K_AxD, AxA = K_AxA, DxD = K_DxD, AxAxM = K_AxAxM, DxDxM = K_DxDxM, AxDxM = K_AxDxM, M = K_M)
              }
              method_names <- c("GBLUP", "rrBLUP_markers", "RKHS_BGLR", "BRR", "BayesA", "BayesB", "BayesC", "BL")
            }

            if (multitraitGP == TRUE){
              if(gp_model == "GBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  prepare_covariates <- function(Xcov) {
                    if (is.list(Xcov) && !is.data.frame(Xcov)) {
                      Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                    } else {
                      Xcov_df <- as.data.frame(Xcov)
                    }

                    # remove intercept to avoid duplication
                    Xcov_mat <- model.matrix(~ . - 1, data = Xcov_df)

                    # drop duplicates
                    Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                    # drop near-zero variance
                    nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                    if (any(nzv)) {
                      message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                              paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                    }

                    # enforce full rank
                    qrX <- qr(Xcov_mat)
                    if (qrX$rank < ncol(Xcov_mat)) {
                      drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                      message("Dropping collinear covariates: ",
                              paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
                    }

                    if (ncol(Xcov_mat) == 0) return(NULL)
                    return(Xcov_mat)
                  }

                  # GBLUP  with rrBLUP package
                  pred_list <- list()
                  for (kernel_name in names(kernels)) {
                    K <- kernels[[kernel_name]]
                    covTraits <- colnames(Y.covmasked)
                    gebv_list <- list()
                    for (covt in covTraits){
                      y <- Y.covmasked[[covt]]
                      # remove NA individuals (mixed.solve handles but safer)
                      ok <- !is.na(y)
                      sol <- mixed.solve(y = y[ok], K = K[ok, ok, drop = FALSE])
                      # sol$u is GEBV vector for individuals with y; bring back into full n vector
                      u_full <- rep(NA, nrow(Y.covmasked))
                      u_full[which(ok)] <- sol$u[rownames(K)[ok]]   # name-align
                      gebv_list[[covt]] <- u_full
                    }
                    Y.tmasked <- as.data.frame(Y.masked)
                    for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                    Xcov <- Y.tmasked[, -1, drop = FALSE]
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      message("⚠️ No valid covariates left, setting Xcov_mat = NULL")
                      Xcov_mat <- NULL
                    } else {
                      # Optional: impute missing values
                      if (anyNA(Xcov_mat)) {
                        message("Replacing NA covariate values with column means")
                        Xcov_mat <- apply(Xcov_mat, 2, function(col) {
                          ifelse(is.na(col), mean(col, na.rm = TRUE), col)
                        })
                        Xcov_mat <- as.matrix(Xcov_mat)
                      }
                    }
                    # ---- Handle cases where covariates are dropped ----
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL   # rrBLUP accepts NULL
                    } else {
                      # make sure it's a matrix
                      Xcov_mat <- as.matrix(Xcov_mat)

                      # if rownames missing, assign from Y.tmasked
                      if (is.null(rownames(Xcov_mat))) {
                        rownames(Xcov_mat) <- rownames(Y.tmasked)
                      }

                      # if nrow mismatch, force alignment by merging
                      if (nrow(Xcov_mat) != nrow(Y.tmasked)) {
                        Xcov_mat <- Xcov_mat[match(rownames(Y.tmasked), rownames(Xcov_mat)), , drop = FALSE]
                      }

                      stopifnot(nrow(Xcov_mat) == nrow(Y.tmasked))
                    }

                    model_gblup <- rrBLUP::mixed.solve(y = Y.tmasked[,1], K = K, X=Xcov_mat)
                    pred_gblup <- model_gblup$u[test_ids]
                    pred_list[[kernel_name]] <- pred_gblup
                  }
                  pred_gblup_all <- do.call(cbind, pred_list)
                  rownames(pred_gblup_all) <- test_ids

                  #--- Stack GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, gblup  = pred_gblup_all)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- as.formula(paste(trait, "~", paste(colnames(gblup_stack)[-c(1:2)], collapse = " + ")))
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # rrBLUP marker effects model: K_A
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.A_scaled <- as.matrix(geno.A_scaled)
                    mode(geno.A_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.A_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.A_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.A_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_A <- as.vector(Z_test %*% model_rrblup$u)


                    # rrBLUP marker effects model: K_D
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.D_scaled <- as.matrix(geno.D_scaled)
                    mode(geno.D_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.D_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.D_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.D_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_D <- as.vector(Z_test %*% model_rrblup$u)

                    # --- C: Stacked rrBLUP outputs via glm ---
                    rrblup_stack <- data.frame(y = y_test, A  = pred_rrblup_A, D  = pred_rrblup_D)
                    colnames(rrblup_stack)[1:2] <- colnames(y_test)
                    formula_rrblup <- paste0(trait," ~ A + D")
                    stack_fit <- lm(formula_rrblup, data = rrblup_stack)
                    pred_rrblup <- predict(stack_fit, newdata = rrblup_stack)


                    # RHKs with BGLR package
                    pred_list <- list()
                    for (kernel_name in names(kernels)) {
                      K <- kernels[[kernel_name]]
                      covTraits <- colnames(Y.covmasked)
                      gebv_list <- list()
                      for(covt in covTraits){
                        y <- Y.covmasked[[covt]]
                        ok <- !is.na(y)
                        fm <- BGLR(y = y[ok],  ETA = list(list(K=K[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                        # Create full-length vector aligned with original y
                        gebv_full <- rep(NA, length(y))
                        gebv_full[ok] <- as.numeric(b)
                        gebv_list[[covt]] <- gebv_full
                      }
                      Y.tmasked <- as.data.frame(Y.masked)
                      for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      Xcov_mat <- as.matrix(Y.tmasked[, -1, drop = FALSE])
                      train_idx <- which(!is.na(Y.tmasked[,1]))
                      test_idx  <- which(is.na(Y.tmasked[,1]))
                      y_train <- Y.tmasked[train_idx, 1]
                      K_train <- K[train_idx, train_idx]
                      X_train <- Xcov_mat[train_idx, , drop = FALSE]
                      fit <- BGLR(y = y_train, ETA = list(list(K = K_train, model = "RKHS"), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      u_train <- as.numeric(fit$ETA[[1]]$u)
                      K_test_train <- K[test_idx, train_idx]
                      pred_test <- K_test_train %*% u_train
                      pred_list[[kernel_name]] <- as.numeric(pred_test)
                    }
                    pred_rkhs_all <- do.call(cbind, pred_list)
                    rownames(pred_rkhs_all) <- rownames(Y.tmasked)[test_idx]

                    #--- Stack RKHS models ---
                    y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                    rkhs_stack <- data.frame(y = y_test, rkhs  = pred_rkhs_all)
                    colnames(rkhs_stack)[1:2] <- colnames(y_test)
                    formula_rkhs <- as.formula(paste(trait, "~", paste(colnames(rkhs_stack)[-c(1:2)], collapse = " + ")))
                    fit_stack  <- lm(formula_rkhs, data = rkhs_stack)
                    pred_rkhs <- predict(fit_stack, newdata = rkhs_stack)


                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                      run_independent_bayes <- function(model_name, Y.masked, Y.covmasked, geno.A_scaled, geno.D_scaled, nIter, burnIn) {
                        covTraits <- colnames(Y.covmasked)
                        # ---- Additive step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = geno.A_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(geno.A_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit_A <- BGLR(y = y_train, ETA = list(list(X = geno.A_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all_A <- rep(NA, nrow(Y.tmasked))
                        yHat_all_A[train_idx] <- fit_A$yHat
                        b_markersA <- fit_A$ETA[[1]]$b
                        b_fixedA   <- fit_A$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markersA)) geno.A_scaled[test_idx, , drop = FALSE] %*% b_markersA else 0
                        fixed_part  <- if (!is.null(b_fixedA))   X_test %*% b_fixedA else 0
                        pred_A <- as.numeric(marker_part + fixed_part)
                        # pred_A <- as.numeric(geno.A_scaled[test_idx, , drop = FALSE] %*% b_markersA + X_test %*% b_fixedA)
                        # ---- Dominance step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = geno.D_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(geno.D_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit_D <- BGLR(y = y_train, ETA = list( list(X = geno.D_scaled[train_idx, , drop = FALSE], model = model_name), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all_D <- rep(NA, nrow(Y.tmasked))
                        yHat_all_D[train_idx] <- fit_D$yHat
                        b_markersD <- fit_D$ETA[[1]]$b
                        b_fixedD   <- fit_D$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markersD)) geno.D_scaled[test_idx, , drop = FALSE] %*% b_markersD else 0
                        fixed_part  <- if (!is.null(b_fixedD))   X_test %*% b_fixedD else 0
                        pred_D <- as.numeric(marker_part + fixed_part)
                        # pred_D <- as.numeric(geno.D_scaled[test_idx, , drop = FALSE] %*% b_markersD + X_test %*% b_fixedD)
                        # ---- Return predictions for test set ----
                        return(data.frame(pred_A = pred_A, pred_D = pred_D, row.names = rownames(Y.tmasked)[test_idx]))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, Y.covmasked, geno.A_scaled, geno.D_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(n.cores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "Y.covmasked", "geno.A_scaled", "geno.D_scaled", "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked, Y.covmasked = Y.covmasked,geno.A_scaled = geno.A_scaled, geno.D_scaled = geno.D_scaled, nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }

                      # Run all Bayesian models in parallel
                      preds_stack <- run_parallel_stack(Y.masked = Y.masked, Y.covmasked = Y.covmasked, geno.A_scaled = geno.A_scaled, geno.D_scaled = geno.D_scaled, nIter = nIter, burnIn = burnIn, n.cores = ncores)

                      # Stack predictions using lm for each model
                      stacked_preds <- lapply(preds_stack, function(pred_df) {
                        test_rows <- which(is.na(Y.masked))
                        y_test <- as.numeric(Y.raw[test_rows, 2])

                        stack_df <- data.frame(y = y_test, pred_A = pred_df$pred_A, pred_D = pred_df$pred_D)
                        stack_df <- stack_df[complete.cases(stack_df), ]

                        if (nrow(stack_df) == 0) {
                          warning("No valid rows to fit linear model.")
                          return(rep(NA, length(test_rows)))
                        }

                        # Stacking model: combine pred_A and pred_D
                        fit <- lm(y ~ pred_A + pred_D, data = stack_df)
                        predict(fit, newdata = stack_df)
                      })

                      # stacked_preds is a named list of prediction vectors for each bayes model
                      names(stacked_preds) <- bayes_models
                      pred_brr <- data.frame(ID = test_ids, Prediction = stacked_preds[["BRR"]])
                      pred_bayesA <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesA"]])
                      pred_bayesB <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesB"]])
                      pred_bayesC <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesC"]])
                      pred_bayesLasso <- data.frame(ID = test_ids, Prediction = stacked_preds[["BL"]])

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr$Prediction, BayesA = pred_bayesA$Prediction,
                                                      BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr$Prediction,
                                                      BayesA = pred_bayesA$Prediction, BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                } else {
                  prepare_covariates <- function(Xcov) {
                    if (is.list(Xcov) && !is.data.frame(Xcov)) {
                      Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                    } else {
                      Xcov_df <- as.data.frame(Xcov)
                    }

                    # remove intercept to avoid duplication
                    Xcov_mat <- model.matrix(~ . - 1, data = Xcov_df)

                    # drop duplicates
                    Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                    # drop near-zero variance
                    nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                    if (any(nzv)) {
                      message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                              paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                    }

                    # enforce full rank
                    qrX <- qr(Xcov_mat)
                    if (qrX$rank < ncol(Xcov_mat)) {
                      drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                      message("Dropping collinear covariates: ",
                              paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
                    }

                    if (ncol(Xcov_mat) == 0) return(NULL)
                    return(Xcov_mat)
                  }

                  # GBLUP with rrBLUP package
                  covTraits <- colnames(Y.covmasked)
                  gebv_list <- list()
                  for (covt in covTraits){
                    y <- Y.covmasked[[covt]]
                    # remove NA individuals (mixed.solve handles but safer)
                    ok <- !is.na(y)
                    sol <- mixed.solve(y = y[ok], K = myKIx[ok, ok, drop = FALSE])
                    # sol$u is GEBV vector for individuals with y; bring back into full n vector
                    u_full <- rep(NA, nrow(Y.covmasked))
                    u_full[which(ok)] <- sol$u[rownames(myKIx)[ok]]   # name-align
                    gebv_list[[covt]] <- u_full
                  }
                  Y.tmasked <- as.data.frame(Y.masked)
                  for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                  Xcov <- Y.tmasked[, -1, drop = FALSE]
                  Xcov_mat <- prepare_covariates(Xcov)
                  if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                    message("⚠️ No valid covariates left, setting Xcov_mat = NULL")
                    Xcov_mat <- NULL
                  } else {
                    # Optional: impute missing values
                    if (anyNA(Xcov_mat)) {
                      message("Replacing NA covariate values with column means")
                      Xcov_mat <- apply(Xcov_mat, 2, function(col) {
                        ifelse(is.na(col), mean(col, na.rm = TRUE), col)
                      })
                      Xcov_mat <- as.matrix(Xcov_mat)
                    }
                  }
                  # ---- Handle cases where covariates are dropped ----
                  if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                    Xcov_mat <- NULL   # rrBLUP accepts NULL
                  } else {
                    # make sure it's a matrix
                    Xcov_mat <- as.matrix(Xcov_mat)

                    # if rownames missing, assign from Y.tmasked
                    if (is.null(rownames(Xcov_mat))) {
                      rownames(Xcov_mat) <- rownames(Y.tmasked)
                    }

                    # if nrow mismatch, force alignment by merging
                    if (nrow(Xcov_mat) != nrow(Y.tmasked)) {
                      Xcov_mat <- Xcov_mat[match(rownames(Y.tmasked), rownames(Xcov_mat)), , drop = FALSE]
                    }

                    stopifnot(nrow(Xcov_mat) == nrow(Y.tmasked))
                  }

                  model_gblup <- rrBLUP::mixed.solve(y = Y.tmasked[,1], K = K, X=Xcov_mat)
                  pred_gblup <- model_gblup$u[test_ids]  # genomic breeding values

                  if (!is.null(Additional_models)){
                    # rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    # Ensure genotype matrix is numeric
                    geno_scaled <- as.matrix(geno_scaled)
                    mode(geno_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    # Align test set to training SNPs (safe check)
                    Z_test <- geno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test,
                                    center = colMeans(geno_scaled[train_ids, kept_snps]),
                                    scale = FALSE)
                    pred_rrblup <- as.vector(Z_test %*% model_rrblup$u)

                    # RHKs with BGLR package
                    covTraits <- colnames(Y.covmasked)
                    gebv_list <- list()
                    for(covt in covTraits){
                      y <- Y.covmasked[[covt]]
                      ok <- !is.na(y)
                      fm <- BGLR(y = y[ok],  ETA = list(list(K=myKIx[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                      # Create full-length vector aligned with original y
                      gebv_full <- rep(NA, length(y))
                      gebv_full[ok] <- as.numeric(b)
                      gebv_list[[covt]] <- gebv_full
                    }
                    Y.tmasked <- as.data.frame(Y.masked)
                    for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                    Xcov_mat <- as.matrix(Y.tmasked[, -1, drop = FALSE])
                    train_idx <- which(!is.na(Y.tmasked[,1]))
                    test_idx  <- which(is.na(Y.tmasked[,1]))
                    y_train <- Y.tmasked[train_idx, 1]
                    K_train <- myKIx[train_idx, train_idx]
                    X_train <- Xcov_mat[train_idx, , drop = FALSE]
                    fit <- BGLR(y = y_train, ETA = list(list(K = K_train, model = "RKHS"), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                    u_train <- as.numeric(fit$ETA[[1]]$u)
                    K_test_train <- myKIx[test_idx, train_idx]
                    pred_rkhs <-as.numeric( K_test_train %*% u_train)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                      run_independent_bayes <- function(model_name, Y.masked, Y.covmasked, geno_scaled, nIter, burnIn) {
                        covTraits <- colnames(Y.covmasked)
                        # ---- Additive step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = geno_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(geno_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit <- BGLR(y = y_train, ETA = list(list(X = geno_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all <- rep(NA, nrow(Y.tmasked))
                        yHat_all[train_idx] <- fit$yHat
                        b_markers <- fit$ETA[[1]]$b
                        b_fixed   <- fit$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markers)) geno_scaled[test_idx, , drop = FALSE] %*% b_markers else 0
                        fixed_part  <- if (!is.null(b_fixed))   X_test %*% b_fixed else 0
                        pred <- as.numeric(marker_part + fixed_part)
                        return(data.frame(pred = pred, row.names = rownames(Y.tmasked)[test_idx]))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, Y.covmasked, geno_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(n.cores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "Y.covmasked", "geno_scaled", "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked, Y.covmasked = Y.covmasked, geno_scaled = geno_scaled, nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }

                      # Run all Bayesian models in parallel
                      preds <- run_parallel_stack(Y.masked = Y.masked, Y.covmasked = Y.covmasked, geno_scaled = geno_scaled, nIter = nIter, burnIn = burnIn, n.cores = ncores)

                      names(preds) <- bayes_models
                      pred_brr <- preds$BRR; colnames(pred_brr)[1] <- "Prediction"
                      pred_bayesA <- preds$BayesA; colnames(pred_bayesA)[1] <- "Prediction"
                      pred_bayesB <- preds$BayesB; colnames(pred_bayesB)[1] <- "Prediction"
                      pred_bayesC <- preds$BayesC; colnames(pred_bayesC)[1] <- "Prediction"
                      pred_bayesLasso <- preds$BL; colnames(pred_bayesLasso)[1] <- "Prediction"

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr$Prediction, BayesA = pred_bayesA$Prediction,
                                                      BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr$Prediction,
                                                      BayesA = pred_bayesA$Prediction, BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                }
              }
              if(gp_model == "gBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  prepare_covariates <- function(Xcov) {
                    if (is.list(Xcov) && !is.data.frame(Xcov)) {
                      Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                    } else {
                      Xcov_df <- as.data.frame(Xcov)
                    }

                    # remove intercept to avoid duplication
                    Xcov_mat <- model.matrix(~ . - 1, data = Xcov_df)

                    # drop duplicates
                    Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                    # drop near-zero variance
                    nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                    if (any(nzv)) {
                      message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                              paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                    }

                    # enforce full rank
                    qrX <- qr(Xcov_mat)
                    if (qrX$rank < ncol(Xcov_mat)) {
                      drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                      message("Dropping collinear covariates: ",
                              paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
                    }

                    if (ncol(Xcov_mat) == 0) return(NULL)
                    return(Xcov_mat)
                  }

                  # GBLUP  with rrBLUP package
                  pred_list <- list()
                  for (kernel_name in names(kernels)) {
                    K <- kernels[[kernel_name]]
                    covTraits <- colnames(Y.covmasked)
                    gebv_list <- list()
                    for (covt in covTraits){
                      y <- Y.covmasked[[covt]]
                      # remove NA individuals (mixed.solve handles but safer)
                      ok <- !is.na(y)
                      sol <- mixed.solve(y = y[ok], K = K[ok, ok, drop = FALSE])
                      # sol$u is GEBV vector for individuals with y; bring back into full n vector
                      u_full <- rep(NA, nrow(Y.covmasked))
                      u_full[which(ok)] <- sol$u[rownames(K)[ok]]   # name-align
                      gebv_list[[covt]] <- u_full
                    }
                    Y.tmasked <- as.data.frame(Y.masked)
                    for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                    Xcov <- Y.tmasked[, -1, drop = FALSE]
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      message("⚠️ No valid covariates left, setting Xcov_mat = NULL")
                      Xcov_mat <- NULL
                    } else {
                      # Optional: impute missing values
                      if (anyNA(Xcov_mat)) {
                        message("Replacing NA covariate values with column means")
                        Xcov_mat <- apply(Xcov_mat, 2, function(col) {
                          ifelse(is.na(col), mean(col, na.rm = TRUE), col)
                        })
                        Xcov_mat <- as.matrix(Xcov_mat)
                      }
                    }
                    # ---- Handle cases where covariates are dropped ----
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL   # rrBLUP accepts NULL
                    } else {
                      # make sure it's a matrix
                      Xcov_mat <- as.matrix(Xcov_mat)

                      # if rownames missing, assign from Y.tmasked
                      if (is.null(rownames(Xcov_mat))) {
                        rownames(Xcov_mat) <- rownames(Y.tmasked)
                      }

                      # if nrow mismatch, force alignment by merging
                      if (nrow(Xcov_mat) != nrow(Y.tmasked)) {
                        Xcov_mat <- Xcov_mat[match(rownames(Y.tmasked), rownames(Xcov_mat)), , drop = FALSE]
                      }

                      stopifnot(nrow(Xcov_mat) == nrow(Y.tmasked))
                    }

                    model_gblup <- rrBLUP::mixed.solve(y = Y.tmasked[,1], K = K, X=Xcov_mat)
                    pred_gblup <- model_gblup$u[test_ids]
                    pred_list[[kernel_name]] <- pred_gblup
                  }
                  pred_gblup_all <- do.call(cbind, pred_list)
                  rownames(pred_gblup_all) <- test_ids

                  #--- Stack GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, gblup  = pred_gblup_all)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- as.formula(paste(trait, "~", paste(colnames(gblup_stack[-c(1:2)]), collapse = " + ")))
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # rrBLUP marker effects model: K_M
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_M <- as.vector(Z_test %*% model_rrblup$u)

                    # --- C: Stacked rrBLUP outputs via glm ---
                    rrblup_stack <- data.frame(y = y_test, M  = pred_rrblup_M)
                    colnames(rrblup_stack)[1:2] <- colnames(y_test)
                    formula_rrblup <- paste0(trait," ~ M")
                    stack_fit <- lm(formula_rrblup, data = rrblup_stack)
                    pred_rrblup <- predict(stack_fit, newdata = rrblup_stack)


                    # RHKs with BGLR package
                    pred_list <- list()
                    for (kernel_name in names(kernels)) {
                      K <- kernels[[kernel_name]]
                      covTraits <- colnames(Y.covmasked)
                      gebv_list <- list()
                      for(covt in covTraits){
                        y <- Y.covmasked[[covt]]
                        ok <- !is.na(y)
                        fm <- BGLR(y = y[ok],  ETA = list(list(K=K[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                        # Create full-length vector aligned with original y
                        gebv_full <- rep(NA, length(y))
                        gebv_full[ok] <- as.numeric(b)
                        gebv_list[[covt]] <- gebv_full
                      }
                      Y.tmasked <- as.data.frame(Y.masked)
                      for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      Xcov_mat <- as.matrix(Y.tmasked[, -1, drop = FALSE])
                      train_idx <- which(!is.na(Y.tmasked[,1]))
                      test_idx  <- which(is.na(Y.tmasked[,1]))
                      y_train <- Y.tmasked[train_idx, 1]
                      K_train <- K[train_idx, train_idx]
                      X_train <- Xcov_mat[train_idx, , drop = FALSE]
                      fit <- BGLR(y = y_train, ETA = list(list(K = K_train, model = "RKHS"), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      u_train <- as.numeric(fit$ETA[[1]]$u)
                      K_test_train <- K[test_idx, train_idx]
                      pred_test <- K_test_train %*% u_train
                      pred_list[[kernel_name]] <- as.numeric(pred_test)
                    }
                    pred_rkhs_all <- do.call(cbind, pred_list)
                    rownames(pred_rkhs_all) <- rownames(Y.tmasked)[test_idx]

                    #--- Stack RKHS models ---
                    y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                    rkhs_stack <- data.frame(y = y_test, rkhs  = pred_rkhs_all)
                    colnames(rkhs_stack)[1:2] <- colnames(y_test)
                    formula_rkhs <- as.formula(paste(trait, "~", paste(colnames(rkhs_stack)[-c(1:2)], collapse = " + ")))
                    fit_stack  <- lm(formula_rkhs, data = rkhs_stack)
                    pred_rkhs <- predict(fit_stack, newdata = rkhs_stack)


                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                      run_independent_bayes <- function(model_name, Y.masked, Y.covmasked, mgeno_scaled, nIter, burnIn) {
                        covTraits <- colnames(Y.covmasked)
                        # ---- Metagenome step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = mgeno_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(mgeno_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit_M <- BGLR(y = y_train, ETA = list(list(X = mgeno_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all_M <- rep(NA, nrow(Y.tmasked))
                        yHat_all_M[train_idx] <- fit_M$yHat
                        b_markersM <- fit_M$ETA[[1]]$b
                        b_fixedM   <- fit_M$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markersM)) mgeno_scaled[test_idx, , drop = FALSE] %*% b_markersM else 0
                        fixed_part  <- if (!is.null(b_fixedM))   X_test %*% b_fixedM else 0
                        pred_M <- as.numeric(marker_part + fixed_part)

                        # ---- Return predictions for test set ----
                        return(data.frame(pred_M = pred_M, row.names = rownames(Y.tmasked)[test_idx]))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, Y.covmasked, mgeno_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(n.cores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "Y.covmasked", "mgeno_scaled", "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked, Y.covmasked = Y.covmasked,mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }

                      # Run all Bayesian models in parallel
                      preds_stack <- run_parallel_stack(Y.masked = Y.masked, Y.covmasked = Y.covmasked, mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn, n.cores = ncores)

                      # Stack predictions using lm for each model
                      stacked_preds <- lapply(preds_stack, function(pred_df) {
                        test_rows <- which(is.na(Y.masked))
                        y_test <- as.numeric(Y.raw[test_rows, 2])

                        stack_df <- data.frame(y = y_test, pred_M = pred_df$pred_M)
                        stack_df <- stack_df[complete.cases(stack_df), ]

                        if (nrow(stack_df) == 0) {
                          warning("No valid rows to fit linear model.")
                          return(rep(NA, length(test_rows)))
                        }

                        # Stacking model: combine pred_M and pred_D
                        fit <- lm(y ~ pred_M, data = stack_df)
                        predict(fit, newdata = stack_df)
                      })

                      # stacked_preds is a named list of prediction vectors for each bayes model
                      names(stacked_preds) <- bayes_models
                      pred_brr <- data.frame(ID = test_ids, Prediction = stacked_preds[["BRR"]])
                      pred_bayesA <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesA"]])
                      pred_bayesB <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesB"]])
                      pred_bayesC <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesC"]])
                      pred_bayesLasso <- data.frame(ID = test_ids, Prediction = stacked_preds[["BL"]])

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr[,2], BayesA = pred_bayesA[,2],
                                                      BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr[,2],
                                                      BayesA = pred_bayesA[,2], BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                } else {
                  prepare_covariates <- function(Xcov) {
                    if (is.list(Xcov) && !is.data.frame(Xcov)) {
                      Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                    } else {
                      Xcov_df <- as.data.frame(Xcov)
                    }

                    # remove intercept to avoid duplication
                    Xcov_mat <- model.matrix(~ . - 1, data = Xcov_df)

                    # drop duplicates
                    Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                    # drop near-zero variance
                    nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                    if (any(nzv)) {
                      message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                              paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                    }

                    # enforce full rank
                    qrX <- qr(Xcov_mat)
                    if (qrX$rank < ncol(Xcov_mat)) {
                      drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                      message("Dropping collinear covariates: ",
                              paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
                    }

                    if (ncol(Xcov_mat) == 0) return(NULL)
                    return(Xcov_mat)
                  }

                  # GBLUP with rrBLUP package
                  covTraits <- colnames(Y.covmasked)
                  gebv_list <- list()
                  for (covt in covTraits){
                    y <- Y.covmasked[[covt]]
                    # remove NA individuals (mixed.solve handles but safer)
                    ok <- !is.na(y)
                    sol <- mixed.solve(y = y[ok], K = metagKIx[ok, ok, drop = FALSE])
                    # sol$u is GEBV vector for individuals with y; bring back into full n vector
                    u_full <- rep(NA, nrow(Y.covmasked))
                    u_full[which(ok)] <- sol$u[rownames(metagKIx)[ok]]   # name-align
                    gebv_list[[covt]] <- u_full
                  }
                  Y.tmasked <- as.data.frame(Y.masked)
                  for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                  Xcov <- Y.tmasked[, -1, drop = FALSE]
                  Xcov_mat <- prepare_covariates(Xcov)
                  if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                    message("⚠️ No valid covariates left, setting Xcov_mat = NULL")
                    Xcov_mat <- NULL
                  } else {
                    # Optional: impute missing values
                    if (anyNA(Xcov_mat)) {
                      message("Replacing NA covariate values with column means")
                      Xcov_mat <- apply(Xcov_mat, 2, function(col) {
                        ifelse(is.na(col), mean(col, na.rm = TRUE), col)
                      })
                      Xcov_mat <- as.matrix(Xcov_mat)
                    }
                  }
                  # ---- Handle cases where covariates are dropped ----
                  if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                    Xcov_mat <- NULL   # rrBLUP accepts NULL
                  } else {
                    # make sure it's a matrix
                    Xcov_mat <- as.matrix(Xcov_mat)

                    # if rownames missing, assign from Y.tmasked
                    if (is.null(rownames(Xcov_mat))) {
                      rownames(Xcov_mat) <- rownames(Y.tmasked)
                    }

                    # if nrow mismatch, force alignment by merging
                    if (nrow(Xcov_mat) != nrow(Y.tmasked)) {
                      Xcov_mat <- Xcov_mat[match(rownames(Y.tmasked), rownames(Xcov_mat)), , drop = FALSE]
                    }

                    stopifnot(nrow(Xcov_mat) == nrow(Y.tmasked))
                  }

                  model_gblup <- rrBLUP::mixed.solve(y = Y.tmasked[,1], K = K, X=Xcov_mat)
                  pred_gblup <- model_gblup$u[test_ids]  # genomic breeding values

                  if (!is.null(Additional_models)){
                    # rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    # Ensure genotype matrix is numeric
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (ensure full rank)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test,
                                    center = colMeans(mgeno_scaled[train_ids, kept_snps]),
                                    scale = FALSE)
                    pred_rrblup <- as.vector(Z_test %*% model_rrblup$u)

                    # RHKs with BGLR package
                    covTraits <- colnames(Y.covmasked)
                    gebv_list <- list()
                    for(covt in covTraits){
                      y <- Y.covmasked[[covt]]
                      ok <- !is.na(y)
                      fm <- BGLR(y = y[ok],  ETA = list(list(K=metagKIx[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                      # Create full-length vector aligned with original y
                      gebv_full <- rep(NA, length(y))
                      gebv_full[ok] <- as.numeric(b)
                      gebv_list[[covt]] <- gebv_full
                    }
                    Y.tmasked <- as.data.frame(Y.masked)
                    for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                    Xcov_mat <- as.matrix(Y.tmasked[, -1, drop = FALSE])
                    train_idx <- which(!is.na(Y.tmasked[,1]))
                    test_idx  <- which(is.na(Y.tmasked[,1]))
                    y_train <- Y.tmasked[train_idx, 1]
                    K_train <- metagKIx[train_idx, train_idx]
                    X_train <- Xcov_mat[train_idx, , drop = FALSE]
                    fit <- BGLR(y = y_train, ETA = list(list(K = K_train, model = "RKHS"), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                    u_train <- as.numeric(fit$ETA[[1]]$u)
                    K_test_train <- metagKIx[test_idx, train_idx]
                    pred_rkhs <-as.numeric( K_test_train %*% u_train)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                      run_independent_bayes <- function(model_name, Y.masked, Y.covmasked, mgeno_scaled, nIter, burnIn) {
                        covTraits <- colnames(Y.covmasked)
                        # ---- Metagenome step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = mgeno_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(mgeno_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit <- BGLR(y = y_train, ETA = list(list(X = mgeno_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all <- rep(NA, nrow(Y.tmasked))
                        yHat_all[train_idx] <- fit$yHat
                        b_markers <- fit$ETA[[1]]$b
                        b_fixed   <- fit$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markers)) mgeno_scaled[test_idx, , drop = FALSE] %*% b_markers else 0
                        fixed_part  <- if (!is.null(b_fixed))   X_test %*% b_fixed else 0
                        pred <- as.numeric(marker_part + fixed_part)
                        return(data.frame(pred = pred, row.names = rownames(Y.tmasked)[test_idx]))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, Y.covmasked, mgeno_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(n.cores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "Y.covmasked", "mgeno_scaled", "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked, Y.covmasked = Y.covmasked, mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }

                      # Run all Bayesian models in parallel
                      preds <- run_parallel_stack(Y.masked = Y.masked, Y.covmasked = Y.covmasked, mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn, n.cores = ncores)

                      names(preds) <- bayes_models
                      pred_brr <- preds$BRR; colnames(pred_brr)[1] <- "Prediction"
                      pred_bayesA <- preds$BayesA; colnames(pred_bayesA)[1] <- "Prediction"
                      pred_bayesB <- preds$BayesB; colnames(pred_bayesB)[1] <- "Prediction"
                      pred_bayesC <- preds$BayesC; colnames(pred_bayesC)[1] <- "Prediction"
                      pred_bayesLasso <- preds$BL; colnames(pred_bayesLasso)[1] <- "Prediction"

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr$Prediction, BayesA = pred_bayesA$Prediction,
                                                      BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr$Prediction,
                                                      BayesA = pred_bayesA$Prediction, BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                }
              }
              if(gp_model == "gGBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  prepare_covariates <- function(Xcov) {
                    if (is.list(Xcov) && !is.data.frame(Xcov)) {
                      Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                    } else {
                      Xcov_df <- as.data.frame(Xcov)
                    }

                    # remove intercept to avoid duplication
                    Xcov_mat <- model.matrix(~ . - 1, data = Xcov_df)

                    # drop duplicates
                    Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                    # drop near-zero variance
                    nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                    if (any(nzv)) {
                      message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                              paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                    }

                    # enforce full rank
                    qrX <- qr(Xcov_mat)
                    if (qrX$rank < ncol(Xcov_mat)) {
                      drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                      message("Dropping collinear covariates: ",
                              paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
                    }

                    if (ncol(Xcov_mat) == 0) return(NULL)
                    return(Xcov_mat)
                  }

                  # GBLUP  with rrBLUP package
                  pred_list <- list()
                  for (kernel_name in names(kernels)) {
                    K <- kernels[[kernel_name]]
                    covTraits <- colnames(Y.covmasked)
                    gebv_list <- list()
                    for (covt in covTraits){
                      y <- Y.covmasked[[covt]]
                      # remove NA individuals (mixed.solve handles but safer)
                      ok <- !is.na(y)
                      sol <- mixed.solve(y = y[ok], K = K[ok, ok, drop = FALSE])
                      # sol$u is GEBV vector for individuals with y; bring back into full n vector
                      u_full <- rep(NA, nrow(Y.covmasked))
                      u_full[which(ok)] <- sol$u[rownames(K)[ok]]   # name-align
                      gebv_list[[covt]] <- u_full
                    }
                    Y.tmasked <- as.data.frame(Y.masked)
                    for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                    Xcov <- Y.tmasked[, -1, drop = FALSE]
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      message("⚠️ No valid covariates left, setting Xcov_mat = NULL")
                      Xcov_mat <- NULL
                    } else {
                      # Optional: impute missing values
                      if (anyNA(Xcov_mat)) {
                        message("Replacing NA covariate values with column means")
                        Xcov_mat <- apply(Xcov_mat, 2, function(col) {
                          ifelse(is.na(col), mean(col, na.rm = TRUE), col)
                        })
                        Xcov_mat <- as.matrix(Xcov_mat)
                      }
                    }
                    # ---- Handle cases where covariates are dropped ----
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL   # rrBLUP accepts NULL
                    } else {
                      # make sure it's a matrix
                      Xcov_mat <- as.matrix(Xcov_mat)

                      # if rownames missing, assign from Y.tmasked
                      if (is.null(rownames(Xcov_mat))) {
                        rownames(Xcov_mat) <- rownames(Y.tmasked)
                      }

                      # if nrow mismatch, force alignment by merging
                      if (nrow(Xcov_mat) != nrow(Y.tmasked)) {
                        Xcov_mat <- Xcov_mat[match(rownames(Y.tmasked), rownames(Xcov_mat)), , drop = FALSE]
                      }

                      stopifnot(nrow(Xcov_mat) == nrow(Y.tmasked))
                    }

                    model_gblup <- rrBLUP::mixed.solve(y = Y.tmasked[,1], K = K, X=Xcov_mat)
                    pred_gblup <- model_gblup$u[test_ids]
                    pred_list[[kernel_name]] <- pred_gblup
                  }
                  pred_gblup_all <- do.call(cbind, pred_list)
                  rownames(pred_gblup_all) <- test_ids

                  #--- Stack GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, gblup  = pred_gblup_all)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- as.formula(paste(trait, "~", paste(colnames(gblup_stack)[-c(1:2)], collapse = " + ")))
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # rrBLUP marker effects model: K_A
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.A_scaled <- as.matrix(geno.A_scaled)
                    mode(geno.A_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.A_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.A_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.A_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_A <- as.vector(Z_test %*% model_rrblup$u)


                    # rrBLUP marker effects model: K_D
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.D_scaled <- as.matrix(geno.D_scaled)
                    mode(geno.D_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.D_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.D_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.D_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_D <- as.vector(Z_test %*% model_rrblup$u)

                    # rrBLUP marker effects model: K_M
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_M <- as.vector(Z_test %*% model_rrblup$u)

                    # --- C: Stacked rrBLUP outputs via glm ---
                    rrblup_stack <- data.frame(y = y_test, A  = pred_rrblup_A, D  = pred_rrblup_D, M  = pred_rrblup_M)
                    colnames(rrblup_stack)[1:2] <- colnames(y_test)
                    formula_rrblup <- paste0(trait," ~ A + D + M")
                    stack_fit <- lm(formula_rrblup, data = rrblup_stack)
                    pred_rrblup <- predict(stack_fit, newdata = rrblup_stack)


                    # RHKs with BGLR package
                    pred_list <- list()
                    for (kernel_name in names(kernels)) {
                      K <- kernels[[kernel_name]]
                      covTraits <- colnames(Y.covmasked)
                      gebv_list <- list()
                      for(covt in covTraits){
                        y <- Y.covmasked[[covt]]
                        ok <- !is.na(y)
                        fm <- BGLR(y = y[ok],  ETA = list(list(K=K[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                        # Create full-length vector aligned with original y
                        gebv_full <- rep(NA, length(y))
                        gebv_full[ok] <- as.numeric(b)
                        gebv_list[[covt]] <- gebv_full
                      }
                      Y.tmasked <- as.data.frame(Y.masked)
                      for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      Xcov_mat <- as.matrix(Y.tmasked[, -1, drop = FALSE])
                      train_idx <- which(!is.na(Y.tmasked[,1]))
                      test_idx  <- which(is.na(Y.tmasked[,1]))
                      y_train <- Y.tmasked[train_idx, 1]
                      K_train <- K[train_idx, train_idx]
                      X_train <- Xcov_mat[train_idx, , drop = FALSE]
                      fit <- BGLR(y = y_train, ETA = list(list(K = K_train, model = "RKHS"), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      u_train <- as.numeric(fit$ETA[[1]]$u)
                      K_test_train <- K[test_idx, train_idx]
                      pred_test <- K_test_train %*% u_train
                      pred_list[[kernel_name]] <- as.numeric(pred_test)
                    }
                    pred_rkhs_all <- do.call(cbind, pred_list)
                    rownames(pred_rkhs_all) <- rownames(Y.tmasked)[test_idx]

                    #--- Stack RKHS models ---
                    y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                    rkhs_stack <- data.frame(y = y_test, rkhs  = pred_rkhs_all)
                    colnames(rkhs_stack)[1:2] <- colnames(y_test)
                    formula_rkhs <- as.formula(paste(trait, "~", paste(colnames(rkhs_stack)[-c(1:2)], collapse = " + ")))
                    fit_stack  <- lm(formula_rkhs, data = rkhs_stack)
                    pred_rkhs <- predict(fit_stack, newdata = rkhs_stack)


                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                      run_independent_bayes <- function(model_name, Y.masked, Y.covmasked, geno.A_scaled, geno.D_scaled, mgeno_scaled, nIter, burnIn) {
                        covTraits <- colnames(Y.covmasked)
                        # ---- Additive step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = geno.A_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(geno.A_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit_A <- BGLR(y = y_train, ETA = list(list(X = geno.A_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all_A <- rep(NA, nrow(Y.tmasked))
                        yHat_all_A[train_idx] <- fit_A$yHat
                        b_markersA <- fit_A$ETA[[1]]$b
                        b_fixedA   <- fit_A$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markersA)) geno.A_scaled[test_idx, , drop = FALSE] %*% b_markersA else 0
                        fixed_part  <- if (!is.null(b_fixedA))   X_test %*% b_fixedA else 0
                        pred_A <- as.numeric(marker_part + fixed_part)
                        # pred_A <- as.numeric(geno.A_scaled[test_idx, , drop = FALSE] %*% b_markersA + X_test %*% b_fixedA)
                        # ---- Dominance step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = geno.D_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(geno.D_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit_D <- BGLR(y = y_train, ETA = list( list(X = geno.D_scaled[train_idx, , drop = FALSE], model = model_name), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all_D <- rep(NA, nrow(Y.tmasked))
                        yHat_all_D[train_idx] <- fit_D$yHat
                        b_markersD <- fit_D$ETA[[1]]$b
                        b_fixedD   <- fit_D$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markersD)) geno.D_scaled[test_idx, , drop = FALSE] %*% b_markersD else 0
                        fixed_part  <- if (!is.null(b_fixedD))   X_test %*% b_fixedD else 0
                        pred_D <- as.numeric(marker_part + fixed_part)
                        # pred_D <- as.numeric(geno.D_scaled[test_idx, , drop = FALSE] %*% b_markersD + X_test %*% b_fixedD)
                        # ---- Metagenome step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = mgeno_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(mgeno_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit_M <- BGLR(y = y_train, ETA = list(list(X = mgeno_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all_M <- rep(NA, nrow(Y.tmasked))
                        yHat_all_M[train_idx] <- fit_M$yHat
                        b_markersM <- fit_M$ETA[[1]]$b
                        b_fixedM   <- fit_M$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markersM)) mgeno_scaled[test_idx, , drop = FALSE] %*% b_markersM else 0
                        fixed_part  <- if (!is.null(b_fixedM))   X_test %*% b_fixedM else 0
                        pred_M <- as.numeric(marker_part + fixed_part)
                        # ---- Return predictions for test set ----
                        return(data.frame(pred_A = pred_A, pred_D = pred_D, pred_M = pred_M, row.names = rownames(Y.tmasked)[test_idx]))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, Y.covmasked, geno.A_scaled, geno.D_scaled, mgeno_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(n.cores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "Y.covmasked", "geno.A_scaled", "geno.D_scaled", "mgeno_scaled", "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked, Y.covmasked = Y.covmasked,geno.A_scaled = geno.A_scaled, geno.D_scaled = geno.D_scaled, mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }

                      # Run all Bayesian models in parallel
                      preds_stack <- run_parallel_stack(Y.masked = Y.masked, Y.covmasked = Y.covmasked, geno.A_scaled = geno.A_scaled, geno.D_scaled = geno.D_scaled, mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn, n.cores = ncores)

                      # Stack predictions using lm for each model
                      stacked_preds <- lapply(preds_stack, function(pred_df) {
                        test_rows <- which(is.na(Y.masked))
                        y_test <- as.numeric(Y.raw[test_rows, 2])

                        stack_df <- data.frame(y = y_test, pred_A = pred_df$pred_A, pred_D = pred_df$pred_D, pred_M = pred_df$pred_M)
                        stack_df <- stack_df[complete.cases(stack_df), ]

                        if (nrow(stack_df) == 0) {
                          warning("No valid rows to fit linear model.")
                          return(rep(NA, length(test_rows)))
                        }

                        # Stacking model: combine pred_A and pred_D
                        fit <- lm(y ~ pred_A + pred_D + pred_M, data = stack_df)
                        predict(fit, newdata = stack_df)
                      })

                      # stacked_preds is a named list of prediction vectors for each bayes model
                      names(stacked_preds) <- bayes_models
                      pred_brr <- data.frame(ID = test_ids, Prediction = stacked_preds[["BRR"]])
                      pred_bayesA <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesA"]])
                      pred_bayesB <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesB"]])
                      pred_bayesC <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesC"]])
                      pred_bayesLasso <- data.frame(ID = test_ids, Prediction = stacked_preds[["BL"]])

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr$Prediction, BayesA = pred_bayesA$Prediction,
                                                      BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr$Prediction,
                                                      BayesA = pred_bayesA$Prediction, BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                } else {
                  prepare_covariates <- function(Xcov) {
                    if (is.list(Xcov) && !is.data.frame(Xcov)) {
                      Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                    } else {
                      Xcov_df <- as.data.frame(Xcov)
                    }

                    # remove intercept to avoid duplication
                    Xcov_mat <- model.matrix(~ . - 1, data = Xcov_df)

                    # drop duplicates
                    Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                    # drop near-zero variance
                    nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                    if (any(nzv)) {
                      message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                              paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                    }

                    # enforce full rank
                    qrX <- qr(Xcov_mat)
                    if (qrX$rank < ncol(Xcov_mat)) {
                      drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                      message("Dropping collinear covariates: ",
                              paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                      Xcov_mat <- Xcov_mat[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
                    }

                    if (ncol(Xcov_mat) == 0) return(NULL)
                    return(Xcov_mat)
                  }

                  # GBLUP with rrBLUP package
                  covTraits <- colnames(Y.covmasked)
                  gebv_list <- list()
                  for (covt in covTraits){
                    y <- Y.covmasked[[covt]]
                    # remove NA individuals (mixed.solve handles but safer)
                    ok <- !is.na(y)
                    sol <- mixed.solve(y = y[ok], K = myKIx[ok, ok, drop = FALSE])
                    # sol$u is GEBV vector for individuals with y; bring back into full n vector
                    u_full <- rep(NA, nrow(Y.covmasked))
                    u_full[which(ok)] <- sol$u[rownames(myKIx)[ok]]   # name-align
                    gebv_list[[covt]] <- u_full
                  }
                  Y.tmasked <- as.data.frame(Y.masked)
                  for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                  Xcov <- Y.tmasked[,-1]
                  model_gblup_g <- rrBLUP::mixed.solve(y = Y.tmasked[,1], K = myKIx, X=Xcov)
                  pred_gblup_g <- model_gblup_g$u[test_ids]  # genomic breeding values
                  gebv_list <- list()
                  for (covt in covTraits){
                    y <- Y.covmasked[[covt]]
                    # remove NA individuals (mixed.solve handles but safer)
                    ok <- !is.na(y)
                    sol <- mixed.solve(y = y[ok], K =metagKIx[ok, ok, drop = FALSE])
                    # sol$u is GEBV vector for individuals with y; bring back into full n vector
                    u_full <- rep(NA, nrow(Y.covmasked))
                    u_full[which(ok)] <- sol$u[rownames(metagKIx)[ok]]   # name-align
                    gebv_list[[covt]] <- u_full
                  }
                  Y.tmasked <- as.data.frame(Y.masked)
                  for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                  Xcov <- Y.tmasked[, -1, drop = FALSE]
                  Xcov_mat <- prepare_covariates(Xcov)
                  if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                    message("⚠️ No valid covariates left, setting Xcov_mat = NULL")
                    Xcov_mat <- NULL
                  } else {
                    # Optional: impute missing values
                    if (anyNA(Xcov_mat)) {
                      message("Replacing NA covariate values with column means")
                      Xcov_mat <- apply(Xcov_mat, 2, function(col) {
                        ifelse(is.na(col), mean(col, na.rm = TRUE), col)
                      })
                      Xcov_mat <- as.matrix(Xcov_mat)
                    }
                  }
                  # ---- Handle cases where covariates are dropped ----
                  if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                    Xcov_mat <- NULL   # rrBLUP accepts NULL
                  } else {
                    # make sure it's a matrix
                    Xcov_mat <- as.matrix(Xcov_mat)

                    # if rownames missing, assign from Y.tmasked
                    if (is.null(rownames(Xcov_mat))) {
                      rownames(Xcov_mat) <- rownames(Y.tmasked)
                    }

                    # if nrow mismatch, force alignment by merging
                    if (nrow(Xcov_mat) != nrow(Y.tmasked)) {
                      Xcov_mat <- Xcov_mat[match(rownames(Y.tmasked), rownames(Xcov_mat)), , drop = FALSE]
                    }

                    stopifnot(nrow(Xcov_mat) == nrow(Y.tmasked))
                  }

                  model_gblup_m <- rrBLUP::mixed.solve(y = Y.tmasked[,1], K = metagKIx, X=Xcov_mat)
                  pred_gblup_m <- model_gblup_m$u[test_ids]  # genomic breeding values

                  #--- Define stacking function for two GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, geno  = pred_gblup_g, micro = pred_gblup_m)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- paste0(trait," ~ geno + micro")
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # metagenomic RHKs with BGLR package
                    covTraits <- colnames(Y.covmasked)
                    gebv_list <- list()
                    for(covt in covTraits){
                      y <- Y.covmasked[[covt]]
                      ok <- !is.na(y)
                      fm <- BGLR(y = y[ok],  ETA = list(list(K=metagKIx[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                      # Create full-length vector aligned with original y
                      gebv_full <- rep(NA, length(y))
                      gebv_full[ok] <- as.numeric(b)
                      gebv_list[[covt]] <- gebv_full
                    }
                    Y.tmasked <- as.data.frame(Y.masked)
                    for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                    Xcov_mat <- as.matrix(Y.tmasked[, -1, drop = FALSE])
                    train_idx <- which(!is.na(Y.tmasked[,1]))
                    test_idx  <- which(is.na(Y.tmasked[,1]))
                    y_train <- Y.tmasked[train_idx, 1]
                    K_train <- metagKIx[train_idx, train_idx]
                    X_train <- Xcov_mat[train_idx, , drop = FALSE]
                    fit <- BGLR(y = y_train, ETA = list(list(K = K_train, model = "RKHS"), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                    u_train <- as.numeric(fit$ETA[[1]]$u)
                    K_test_train <- metagKIx[test_idx, train_idx]
                    pred_rkhs_m <-as.numeric( K_test_train %*% u_train)

                    # genomic RHKs with BGLR package
                    gebv_list <- list()
                    for(covt in covTraits){
                      y <- Y.covmasked[[covt]]
                      ok <- !is.na(y)
                      fm <- BGLR(y = y[ok],  ETA = list(list(K=myKIx[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                      # Create full-length vector aligned with original y
                      gebv_full <- rep(NA, length(y))
                      gebv_full[ok] <- as.numeric(b)
                      gebv_list[[covt]] <- gebv_full
                    }
                    Y.tmasked <- as.data.frame(Y.masked)
                    for (covt in covTraits) { Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                    Xcov_mat <- as.matrix(Y.tmasked[, -1, drop = FALSE])
                    train_idx <- which(!is.na(Y.tmasked[,1]))
                    test_idx  <- which(is.na(Y.tmasked[,1]))
                    y_train <- Y.tmasked[train_idx, 1]
                    K_train <- myKIx[train_idx, train_idx]
                    X_train <- Xcov_mat[train_idx, , drop = FALSE]
                    fit <- BGLR(y = y_train, ETA = list(list(K = K_train, model = "RKHS"), list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                    u_train <- as.numeric(fit$ETA[[1]]$u)
                    K_test_train <- myKIx[test_idx, train_idx]
                    pred_rkhs_g <-as.numeric( K_test_train %*% u_train)
                    #--- Stack RKHS models ---
                    y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                    rkhs_stack <- data.frame(y = y_test, geno  = pred_rkhs_g, micro = pred_rkhs_m)
                    colnames(rkhs_stack)[1:2] <- colnames(y_test)
                    formula_rkhs <- paste0(trait," ~ geno + micro")
                    fit_stack  <- lm(formula_rkhs, data = rkhs_stack)
                    pred_rkhs <- predict(fit_stack, newdata = rkhs_stack)

                    # metagenomic rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    if (is.list(Xcov) && !is.data.frame(Xcov)) {
                      Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                    } else {
                      Xcov_df <- as.data.frame(Xcov)
                    }
                    Xcov_mat <- model.matrix(~ ., data = Xcov_df)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X=Xcov_mat)
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_m <- as.vector(Z_test %*% model_rrblup$u)

                    # genomic rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno_scaled <- as.matrix(geno_scaled)
                    mode(geno_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_g <- as.vector(Z_test %*% model_rrblup$u)

                    # metagenomic rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    # --- SNP preparation ---
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))
                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    # --- Covariate preparation ---
                    prepare_covariates <- function(Xcov) {
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }

                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)

                      # Drop duplicate columns
                      Xcov_mat <- Xcov_mat[, !duplicated(colnames(Xcov_mat)), drop = FALSE]

                      # Drop near-zero variance columns
                      nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                      if (any(nzv)) {
                        message("Dropping ", sum(nzv), " near-zero variance covariates: ",
                                paste(colnames(Xcov_mat)[nzv], collapse = ", "))
                        Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                      }

                      # Drop collinear columns (full rank check)
                      qrX <- qr(Xcov_mat)
                      if (qrX$rank < ncol(Xcov_mat)) {
                        drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                        message("Dropping collinear covariates: ",
                                paste(colnames(Xcov_mat)[drop_idx], collapse = ", "))
                        keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                        Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                      }

                      return(Xcov_mat)
                    }
                    Xcov_mat <- prepare_covariates(Xcov)
                    if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                      Xcov_mat <- NULL
                    }
                    # --- rrBLUP model ---
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                    # --- Prediction ---
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_m <- as.vector(Z_test %*% model_rrblup$u)


                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                      run_independent_bayes <- function(model_name, Y.masked, Y.covmasked, geno_scaled, mgeno_scaled, nIter, burnIn) {
                        covTraits <- colnames(Y.covmasked)
                        # ---- Genomic step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = geno_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(geno_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit <- BGLR(y = y_train, ETA = list(list(X = geno_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all <- rep(NA, nrow(Y.tmasked))
                        yHat_all[train_idx] <- fit$yHat
                        b_markers <- fit$ETA[[1]]$b
                        b_fixed   <- fit$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markers)) geno_scaled[test_idx, , drop = FALSE] %*% b_markers else 0
                        fixed_part  <- if (!is.null(b_fixed))   X_test %*% b_fixed else 0
                        pred_g <- as.numeric(marker_part + fixed_part)
                        # ---- Metagenomic step ----
                        gebv_list <- list()
                        for (covt in covTraits) {
                          y <- Y.covmasked[[covt]]
                          ok <- !is.na(y)
                          fm <- BGLR(y = y[ok], ETA = list(list(X = mgeno_scaled[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$b
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(mgeno_scaled[ok, , drop = FALSE] %*% b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.tmasked <- as.data.frame(Y.masked)
                        for (covt in covTraits) {Y.tmasked[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                        train_idx <- which(!is.na(Y.tmasked[,1]))
                        test_idx  <- which(is.na(Y.tmasked[,1]))
                        y_train <- Y.tmasked[train_idx, 1]
                        X_train <- as.matrix(Y.tmasked[train_idx, -1, drop = FALSE])
                        X_test  <- as.matrix(Y.tmasked[test_idx, -1, drop = FALSE])
                        fit <- BGLR(y = y_train, ETA = list(list(X = mgeno_scaled[train_idx, , drop = FALSE], model = model_name),list(X = X_train, model = "FIXED")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        yHat_all <- rep(NA, nrow(Y.tmasked))
                        yHat_all[train_idx] <- fit$yHat
                        b_markers <- fit$ETA[[1]]$b
                        b_fixed   <- fit$ETA[[2]]$b
                        X_test[is.na(X_test)] <- 0
                        marker_part <- if (!is.null(b_markers)) mgeno_scaled[test_idx, , drop = FALSE] %*% b_markers else 0
                        fixed_part  <- if (!is.null(b_fixed))   X_test %*% b_fixed else 0
                        pred_m <- as.numeric(marker_part + fixed_part)
                        return(data.frame(pred_g = pred_g, pred_m = pred_m, row.names = rownames(Y.tmasked)[test_idx]))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, Y.covmasked, geno_scaled, mgeno_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(n.cores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "Y.covmasked", "geno_scaled","mgeno_scaled", "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked, Y.covmasked = Y.covmasked, geno_scaled = geno_scaled, mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }

                      # Run all Bayesian models in parallel
                      preds_stack <- run_parallel_stack(Y.masked = Y.masked, Y.covmasked = Y.covmasked, geno_scaled = geno_scaled, mgeno_scaled = mgeno_scaled, nIter = nIter, burnIn = burnIn, n.cores = ncores)

                      # Stack predictions using lm for each model
                      stacked_preds <- lapply(preds_stack, function(pred_df) {
                        test_rows <- which(is.na(Y.masked))
                        y_test <- as.numeric(Y.raw[test_rows, 2])

                        stack_df <- data.frame(y = y_test, pred_g = pred_df$pred_g, pred_m = pred_df$pred_m)
                        stack_df <- stack_df[complete.cases(stack_df), ]

                        if (nrow(stack_df) == 0) {
                          warning("No valid rows to fit linear model.")
                          return(rep(NA, length(test_rows)))
                        }

                        # Stacking model: combine pred_g and pred_m
                        fit <- lm(y ~ pred_g + pred_m, data = stack_df)
                        predict(fit, newdata = stack_df)
                      })

                      names(stacked_preds) <- bayes_models
                      pred_brr <- as.data.frame(stacked_preds$BRR); colnames(pred_brr)[1] <- "Prediction"
                      pred_bayesA <- as.data.frame(stacked_preds$BayesA); colnames(pred_bayesA)[1] <- "Prediction"
                      pred_bayesB <- as.data.frame(stacked_preds$BayesB); colnames(pred_bayesB)[1] <- "Prediction"
                      pred_bayesC <- as.data.frame(stacked_preds$BayesC); colnames(pred_bayesC)[1] <- "Prediction"
                      pred_bayesLasso <- as.data.frame(stacked_preds$BL); colnames(pred_bayesLasso)[1] <- "Prediction"

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr$Prediction, BayesA = pred_bayesA$Prediction,
                                                      BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr$Prediction,
                                                      BayesA = pred_bayesA$Prediction, BayesB = pred_bayesB$Prediction, BayesC = pred_bayesC$Prediction, BayesLasso = pred_bayesLasso$Prediction)
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                }
              }
            } else {
              if(gp_model == "GBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  # GBLUP  with rrBLUP package
                  pred_list <- list()
                  for (kernel_name in names(kernels)) {
                    K <- kernels[[kernel_name]]
                    model_gblup <- rrBLUP::mixed.solve(y = Y.masked, K = K)
                    pred_gblup <- model_gblup$u[test_ids]
                    pred_list[[kernel_name]] <- pred_gblup
                  }
                  pred_gblup_all <- do.call(cbind, pred_list)
                  rownames(pred_gblup_all) <- test_ids
                  #--- Stack GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, gblup  = pred_gblup_all)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- as.formula(paste(trait, "~", paste(colnames(gblup_stack[,-c(1:2)]), collapse = " + ")))
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # RHKs with BGLR package
                    pred_list <- list()
                    for (kernel_name in names(kernels)) {
                      K <- kernels[[kernel_name]]
                      fit <- BGLR(y = Y.masked, ETA = list(list(K = K, model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      pred_rkhs <- rep(NA, length(Y.masked))
                      pred_rkhs[is.na(Y.masked)] <- fit$yHat[is.na(Y.masked)]
                      pred_rkhs <- as.numeric(pred_rkhs[is.na(Y.masked)])
                      pred_list[[kernel_name]] <- pred_rkhs
                    }
                    pred_rkhs_all <- do.call(cbind, pred_list)
                    rownames(pred_rkhs_all) <- test_ids
                    #--- Stack GBLUP models ---
                    y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                    rkhs_stack <- data.frame(y = y_test, rkhs  = pred_rkhs_all)
                    colnames(rkhs_stack)[1:2] <- colnames(y_test)
                    formula_rkhs <- as.formula(paste(trait, "~", paste(colnames(rkhs_stack[,-c(1:2)]), collapse = " + ")))
                    fit_stack  <- lm(formula_rkhs, data = rkhs_stack)
                    pred_rkhs <- predict(fit_stack, newdata = rkhs_stack)

                    # rrBLUP marker effects model: K_A
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.A_scaled <- as.matrix(geno.A_scaled)
                    mode(geno.A_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.A_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.A_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.A_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_A <- as.vector(Z_test %*% model_rrblup$u)

                    # rrBLUP marker effects model: K_D
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.D_scaled <- as.matrix(geno.D_scaled)
                    mode(geno.D_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.D_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.D_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.D_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_D <- as.vector(Z_test %*% model_rrblup$u)

                    # --- C: Stacked rrBLUP outputs via glm ---
                    rrblup_stack <- data.frame(y = y_test, A  = pred_rrblup_A, D  = pred_rrblup_D)
                    colnames(rrblup_stack)[1:2] <- colnames(y_test)
                    formula_rrblup <- paste0(trait," ~ A + D")
                    stack_fit <- lm(formula_rrblup, data = rrblup_stack)
                    pred_rrblup <- predict(stack_fit, newdata = rrblup_stack)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to run one stacked Bayesian
                      run_independent_bayes <- function(model_name, Y.masked, geno.A_scaled, geno.D_scaled, nIter, burnIn) {
                        fit_A <- BGLR(y = Y.masked, ETA = list(list(X = geno.A_scaled, model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        pred_A <- geno.A_scaled %*% fit_A$ETA[[1]]$b
                        fit_D <- BGLR(y = Y.masked, ETA = list(list(X = geno.D_scaled, model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        pred_D <- geno.D_scaled %*% fit_D$ETA[[1]]$b
                        return(data.frame(pred_A = as.numeric(pred_A), pred_D = as.numeric(pred_D)))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, geno.A_scaled, geno.D_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(ncores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "geno.A_scaled", "geno.D_scaled",
                                                      "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked,
                                                geno.A_scaled = geno.A_scaled, geno.D_scaled = geno.D_scaled,
                                                nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }
                      # Run all Bayesian models in parallel
                      preds_stack <- run_parallel_stack(Y.masked = Y.masked, geno.A_scaled = geno.A_scaled, geno.D_scaled = geno.D_scaled,
                                                        nIter = nIter, burnIn = burnIn, n.cores = ncores)
                      # Stack predictions using lm for each model
                      stacked_preds <- lapply(preds_stack, function(pred_list) {
                        pred_vec <- Reduce(`+`, pred_list)   # element-wise sum
                        test_rows <- which(is.na(Y.masked))
                        y_test <- as.numeric(Y.raw[test_rows, 2])
                        pred_test <- pred_vec[test_rows]
                        stack_df <- data.frame(y = y_test, pred = pred_test)
                        stack_df <- stack_df[complete.cases(stack_df), ]
                        if (nrow(stack_df) == 0) {
                          warning("No valid rows to fit linear model.")
                          return(rep(NA, length(test_rows)))
                        }
                        fit <- lm(y ~ pred, data = stack_df)
                        predict(fit, newdata = stack_df)
                      })
                      # stacked_preds is a named list of prediction vectors for each bayes model
                      names(stacked_preds) <- bayes_models
                      pred_brr <- data.frame(ID = test_ids, Prediction = stacked_preds[["BRR"]])
                      pred_bayesA <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesA"]])
                      pred_bayesB <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesB"]])
                      pred_bayesC <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesC"]])
                      pred_bayesLasso <- data.frame(ID = test_ids, Prediction = stacked_preds[["BL"]])

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr[,2], BayesA = pred_bayesA[,2],
                                                      BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr[,2],
                                                      BayesA = pred_bayesA[,2], BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                } else {
                  # GBLUP with rrBLUP package
                  model_gblup <- rrBLUP::mixed.solve(y = Y.masked, K = myKIx)
                  pred_gblup <- model_gblup$u[test_ids]  # genomic breeding values

                  if (!is.null(Additional_models)){
                    # RHKs with BGLR package
                    fit <- BGLR(y = Y.masked, ETA = list(list(K = myKIx, model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                    pred_rkhs <- rep(NA, length(Y.masked))
                    pred_rkhs[is.na(Y.masked)] <- fit$yHat[is.na(Y.masked)]
                    pred_rkhs <- as.numeric(pred_rkhs[is.na(Y.masked)])

                    # rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno_scaled <- as.matrix(geno_scaled)
                    mode(geno_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup <- as.vector(Z_test %*% model_rrblup$u)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      run_model <- function(model_name) {
                        fit <- BGLR(y = Y.masked,
                                    ETA = list(list(X = geno_scaled, model = model_name)),
                                    nIter = nIter,
                                    burnIn = burnIn,
                                    verbose = FALSE)
                        pred <- as.numeric(geno_scaled %*% fit$ETA[[1]]$b)
                        return(pred)
                      }

                      cl <- makeCluster(ncores)
                      clusterExport(cl, c("BGLR", "Y.masked", "geno_scaled", "nIter", "burnIn", "run_model"))
                      preds <- parLapply(cl, bayes_models, run_model)
                      stopCluster(cl)

                      names(preds) <- bayes_models
                      pred_brr <- data.frame(ID = common_ids, Prediction = preds[["BRR"]]); pred_brr <- pred_brr[pred_brr$ID %in% test_ids, ]
                      pred_bayesA <- data.frame(ID = common_ids, Prediction = preds[["BayesA"]]); pred_bayesA <- pred_bayesA[pred_bayesA$ID %in% test_ids, ]
                      pred_bayesB <- data.frame(ID = common_ids, Prediction = preds[["BayesB"]]); pred_bayesB <- pred_bayesB[pred_bayesB$ID %in% test_ids, ]
                      pred_bayesC <- data.frame(ID = common_ids, Prediction = preds[["BayesC"]]); pred_bayesC <- pred_bayesC[pred_bayesC$ID %in% test_ids, ]
                      pred_bayesLasso <- data.frame(ID = common_ids, Prediction = preds[["BL"]]); pred_bayesLasso <- pred_bayesLasso[pred_bayesLasso$ID %in% test_ids, ]

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr[,2], BayesA = pred_bayesA[,2],
                                                      BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr[,2],
                                                      BayesA = pred_bayesA[,2], BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                }
              }
              if(gp_model == "gBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  # GBLUP with rrBLUP package
                  pred_list <- list()
                  for (kernel_name in names(kernels)) {
                    K <- kernels[[kernel_name]]
                    model_gblup <- rrBLUP::mixed.solve(y = Y.masked, K = K)
                    pred_gblup <- model_gblup$u[test_ids]
                    pred_list[[kernel_name]] <- pred_gblup
                  }
                  pred_gblup_all <- do.call(cbind, pred_list)
                  rownames(pred_gblup_all) <- test_ids
                  #--- Stack GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, gblup  = pred_gblup_all)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- as.formula(paste(trait, "~ M"))
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # RHKs with BGLR package
                    pred_list <- list()
                    for (kernel_name in names(kernels)) {
                      K <- kernels[[kernel_name]]
                      fit <- BGLR(y = Y.masked, ETA = list(list(K = K, model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      pred_rkhs <- rep(NA, length(Y.masked))
                      pred_rkhs[is.na(Y.masked)] <- fit$yHat[is.na(Y.masked)]
                      pred_rkhs <- as.numeric(pred_rkhs[is.na(Y.masked)])
                      pred_list[[kernel_name]] <- pred_rkhs
                    }
                    pred_rkhs_all <- do.call(cbind, pred_list)
                    rownames(pred_rkhs_all) <- test_ids
                    #--- Stack GBLUP models ---
                    y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                    rkhs_stack <- data.frame(y = y_test, rkhs  = pred_rkhs_all)
                    colnames(rkhs_stack)[1:2] <- colnames(y_test)
                    formula_rkhs <- as.formula(paste(trait, "~ M"))
                    fit_stack  <- lm(formula_rkhs, data = rkhs_stack)
                    pred_rkhs <- predict(fit_stack, newdata = rkhs_stack)

                    # rrBLUP marker effects model: m
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_M <- as.vector(Z_test %*% model_rrblup$u)

                    # --- C: Stacked rrBLUP outputs via glm ---
                    rrblup_stack <- data.frame(y = y_test, pred_rrblup_M)
                    colnames(rrblup_stack)[1:2] <- colnames(y_test)
                    formula_rrblup <- paste0(trait," ~ pred_rrblup_M")
                    stack_fit <- lm(formula_rrblup, data = rrblup_stack)
                    pred_rrblup <- predict(stack_fit, newdata = rrblup_stack)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to run one stacked Bayesian
                      run_independent_bayes <- function(model_name, Y.masked, mgeno_scaled, nIter, burnIn) {
                        fit_M <- BGLR(y = Y.masked, ETA = list(list(X = mgeno_scaled, model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        pred_M <- mgeno_scaled %*% fit_M$ETA[[1]]$b
                        return(data.frame(pred_M = as.numeric(pred_M)))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, mgeno_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(ncores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "mgeno_scaled",
                                                      "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked, mgeno_scaled = mgeno_scaled,
                                                nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }
                      # Run all Bayesian models in parallel
                      preds_stack <- run_parallel_stack(Y.masked = Y.masked, mgeno_scaled = mgeno_scaled,
                                                        nIter = nIter, burnIn = burnIn, n.cores = ncores)
                      # Stack predictions using lm for each model
                      stacked_preds <- lapply(preds_stack, function(pred_list) {
                        pred_vec <- Reduce(`+`, pred_list)   # element-wise sum
                        test_rows <- which(is.na(Y.masked))
                        y_test <- as.numeric(Y.raw[test_rows, 2])
                        pred_test <- pred_vec[test_rows]
                        stack_df <- data.frame(y = y_test, pred = pred_test)
                        stack_df <- stack_df[complete.cases(stack_df), ]
                        if (nrow(stack_df) == 0) {
                          warning("No valid rows to fit linear model.")
                          return(rep(NA, length(test_rows)))
                        }
                        fit <- lm(y ~ pred, data = stack_df)
                        predict(fit, newdata = stack_df)
                      })
                      # stacked_preds is a named list of prediction vectors for each bayes model
                      names(stacked_preds) <- bayes_models
                      pred_brr <- data.frame(ID = test_ids, Prediction = stacked_preds[["BRR"]])
                      pred_bayesA <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesA"]])
                      pred_bayesB <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesB"]])
                      pred_bayesC <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesC"]])
                      pred_bayesLasso <- data.frame(ID = test_ids, Prediction = stacked_preds[["BL"]])

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr[,2], BayesA = pred_bayesA[,2],
                                                      BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr[,2],
                                                      BayesA = pred_bayesA[,2], BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                } else {
                  # GBLUP with rrBLUP package
                  model_gblup <- rrBLUP::mixed.solve(y = Y.masked, K = metagKIx)
                  pred_gblup <- model_gblup$u[test_ids]  # genomic breeding values

                  if (!is.null(Additional_models)){
                    # RHKs with BGLR package
                    fit <- BGLR(y = Y.masked, ETA = list(list(K = metagKIx, model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                    pred_rkhs <- rep(NA, length(Y.masked))
                    pred_rkhs[is.na(Y.masked)] <- fit$yHat[is.na(Y.masked)]
                    pred_rkhs <- as.numeric(pred_rkhs[is.na(Y.masked)])

                    # rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup <- as.vector(Z_test %*% model_rrblup$u)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      run_model <- function(model_name) {
                        fit <- BGLR(y = Y.masked,
                                    ETA = list(list(X = mgeno_scaled, model = model_name)),
                                    nIter = nIter,
                                    burnIn = burnIn,
                                    verbose = FALSE)
                        pred <- as.numeric(mgeno_scaled %*% fit$ETA[[1]]$b)
                        return(pred)
                      }

                      cl <- makeCluster(ncores)
                      clusterExport(cl, c("BGLR", "Y.masked", "mgeno_scaled", "nIter", "burnIn", "run_model"))
                      preds <- parLapply(cl, bayes_models, run_model)
                      stopCluster(cl)

                      names(preds) <- bayes_models
                      pred_brr <- data.frame(ID = common_ids, Prediction = preds[["BRR"]]); pred_brr <- pred_brr[pred_brr$ID %in% test_ids, ]
                      pred_bayesA <- data.frame(ID = common_ids, Prediction = preds[["BayesA"]]); pred_bayesA <- pred_bayesA[pred_bayesA$ID %in% test_ids, ]
                      pred_bayesB <- data.frame(ID = common_ids, Prediction = preds[["BayesB"]]); pred_bayesB <- pred_bayesB[pred_bayesB$ID %in% test_ids, ]
                      pred_bayesC <- data.frame(ID = common_ids, Prediction = preds[["BayesC"]]); pred_bayesC <- pred_bayesC[pred_bayesC$ID %in% test_ids, ]
                      pred_bayesLasso <- data.frame(ID = common_ids, Prediction = preds[["BL"]]); pred_bayesLasso <- pred_bayesLasso[pred_bayesLasso$ID %in% test_ids, ]

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr[,2], BayesA = pred_bayesA[,2],
                                                      BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr[,2],
                                                      BayesA = pred_bayesA[,2], BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                }
              }
              if(gp_model == "gGBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  # GBLUP with rrBLUP package
                  pred_list <- list()
                  for (kernel_name in names(kernels)) {
                    K <- kernels[[kernel_name]]
                    model_gblup <- rrBLUP::mixed.solve(y = Y.masked, K = K)
                    pred_gblup <- model_gblup$u[test_ids]
                    pred_list[[kernel_name]] <- pred_gblup
                  }
                  pred_gblup_all <- do.call(cbind, pred_list)
                  rownames(pred_gblup_all) <- test_ids
                  #--- Stack GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, gblup  = pred_gblup_all)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- as.formula(paste(trait, "~", paste(colnames(gblup_stack[,-c(1:2)]), collapse = " + ")))
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # RHKs with BGLR package
                    pred_list <- list()
                    for (kernel_name in names(kernels)) {
                      K <- kernels[[kernel_name]]
                      fit <- BGLR(y = Y.masked, ETA = list(list(K = K, model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      pred_rkhs <- rep(NA, length(Y.masked))
                      pred_rkhs[is.na(Y.masked)] <- fit$yHat[is.na(Y.masked)]
                      pred_rkhs <- as.numeric(pred_rkhs[is.na(Y.masked)])
                      pred_list[[kernel_name]] <- pred_rkhs
                    }
                    pred_rkhs_all <- do.call(cbind, pred_list)
                    rownames(pred_rkhs_all) <- test_ids
                    #--- Stack GBLUP models ---
                    y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                    rkhs_stack <- data.frame(y = y_test, rkhs  = pred_rkhs_all)
                    colnames(rkhs_stack)[1:2] <- colnames(y_test)
                    formula_rkhs <- as.formula(paste(trait, "~", paste(colnames(rkhs_stack[,-c(1:2)]), collapse = " + ")))
                    fit_stack  <- lm(formula_rkhs, data = rkhs_stack)
                    pred_rkhs <- predict(fit_stack, newdata = rkhs_stack)

                    # rrBLUP marker effects model: K_A
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.A_scaled <- as.matrix(geno.A_scaled)
                    mode(geno.A_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.A_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.A_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.A_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_A <- as.vector(Z_test %*% model_rrblup$u)

                    # rrBLUP marker effects model: K_D
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    geno.D_scaled <- as.matrix(geno.D_scaled)
                    mode(geno.D_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(geno.D_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- geno.D_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(geno.D_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_D <- as.vector(Z_test %*% model_rrblup$u)

                    # rrBLUP marker effects model: m
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup_M <- as.vector(Z_test %*% model_rrblup$u)

                    # --- C: Stacked rrBLUP outputs via glm ---
                    rrblup_stack <- data.frame(y = y_test, A  = pred_rrblup_A, D  = pred_rrblup_D, M = pred_rrblup_M)
                    colnames(rrblup_stack)[1:2] <- colnames(y_test)
                    formula_rrblup <- paste0(trait," ~ A + D + M")
                    stack_fit <- lm(formula_rrblup, data = rrblup_stack)
                    pred_rrblup <- predict(stack_fit, newdata = rrblup_stack)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      # Function to run one stacked Bayesian
                      run_independent_bayes <- function(model_name, Y.masked, geno.A_scaled, geno.D_scaled, mgeno_scaled, nIter, burnIn) {
                        fit_A <- BGLR(y = Y.masked, ETA = list(list(X = geno.A_scaled, model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        pred_A <- geno.A_scaled %*% fit_A$ETA[[1]]$b
                        fit_D <- BGLR(y = Y.masked, ETA = list(list(X = geno.D_scaled, model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        pred_D <- geno.D_scaled %*% fit_D$ETA[[1]]$b
                        fit_M <- BGLR(y = Y.masked, ETA = list(list(X = mgeno_scaled, model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        pred_M <- mgeno_scaled %*% fit_M$ETA[[1]]$b
                        return(data.frame(pred_A = as.numeric(pred_A), pred_D = as.numeric(pred_D), pred_M = as.numeric(pred_M)))
                      }

                      # Wrapper to run all models in parallel
                      run_parallel_stack <- function(Y.masked, geno.A_scaled, geno.D_scaled, mgeno_scaled, nIter, burnIn, n.cores = ncores) {
                        cl <- makeCluster(ncores)
                        clusterEvalQ(cl, library(BGLR))
                        clusterExport(cl, varlist = c("Y.masked", "geno.A_scaled", "geno.D_scaled", "mgeno_scaled",
                                                      "nIter", "burnIn", "run_independent_bayes"), envir = environment())
                        preds_list <- parLapply(cl, bayes_models, function(model) {
                          run_independent_bayes(model, Y.masked = Y.masked,
                                                geno.A_scaled = geno.A_scaled, geno.D_scaled = geno.D_scaled, mgeno_scaled = mgeno_scaled,
                                                nIter = nIter, burnIn = burnIn)
                        })
                        stopCluster(cl)
                        names(preds_list) <- bayes_models
                        return(preds_list)
                      }
                      # Run all Bayesian models in parallel
                      preds_stack <- run_parallel_stack(Y.masked = Y.masked, geno.A_scaled = geno.A_scaled,
                                                        geno.D_scaled = geno.D_scaled, mgeno_scaled = mgeno_scaled,
                                                        nIter = nIter, burnIn = burnIn, n.cores = ncores)
                      # Stack predictions using lm for each model
                      stacked_preds <- lapply(preds_stack, function(pred_list) {
                        pred_vec <- Reduce(`+`, pred_list)   # element-wise sum
                        test_rows <- which(is.na(Y.masked))
                        y_test <- as.numeric(Y.raw[test_rows, 2])
                        pred_test <- pred_vec[test_rows]
                        stack_df <- data.frame(y = y_test, pred = pred_test)
                        stack_df <- stack_df[complete.cases(stack_df), ]
                        if (nrow(stack_df) == 0) {
                          warning("No valid rows to fit linear model.")
                          return(rep(NA, length(test_rows)))
                        }
                        fit <- lm(y ~ pred, data = stack_df)
                        predict(fit, newdata = stack_df)
                      })
                      # stacked_preds is a named list of prediction vectors for each bayes model
                      names(stacked_preds) <- bayes_models
                      pred_brr <- data.frame(ID = test_ids, Prediction = stacked_preds[["BRR"]])
                      pred_bayesA <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesA"]])
                      pred_bayesB <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesB"]])
                      pred_bayesC <- data.frame(ID = test_ids, Prediction = stacked_preds[["BayesC"]])
                      pred_bayesLasso <- data.frame(ID = test_ids, Prediction = stacked_preds[["BL"]])

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr[,2], BayesA = pred_bayesA[,2],
                                                      BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr[,2],
                                                      BayesA = pred_bayesA[,2], BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                } else {
                  # GBLUP with rrBLUP package: g
                  model_gblup_g <- rrBLUP::mixed.solve(y = Y.masked, K = myKIx)
                  pred_gblup_g <- model_gblup_g$u[test_ids]  # genomic breeding values

                  # GBLUP with rrBLUP package: m
                  model_gblup_m <- rrBLUP::mixed.solve(y = Y.masked, K = metagKIx)
                  pred_gblup_m <- model_gblup_m$u[test_ids]  # genomic breeding values

                  #--- Define stacking function for two GBLUP models ---
                  y_test <- Y.raw[rownames(Y.raw) %in% test_ids, ]
                  gblup_stack <- data.frame(y = y_test, geno  = pred_gblup_g, micro = pred_gblup_m)
                  colnames(gblup_stack)[1:2] <- colnames(y_test)
                  formula_gblup <- paste0(trait," ~ geno + micro")
                  fit_stack  <- lm(formula_gblup, data = gblup_stack)
                  pred_gblup <- predict(fit_stack, newdata = gblup_stack)

                  if (!is.null(Additional_models)){
                    # RHKs with BGLR package
                    fit <- BGLR(y = Y.masked, ETA = list(list(K = metagKIx, model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                    pred_rkhs <- rep(NA, length(Y.masked))
                    pred_rkhs[is.na(Y.masked)] <- fit$yHat[is.na(Y.masked)]
                    pred_rkhs <- as.numeric(pred_rkhs[is.na(Y.masked)])

                    # rrBLUP marker effects model
                    train_ids <- which(!is.na(Y.masked))
                    y_train <- Y.masked[train_ids]
                    mgeno_scaled <- as.matrix(mgeno_scaled)
                    mode(mgeno_scaled) <- "numeric"
                    prepare_rrblup_matrix <- function(Z, y) {
                      Z <- as.matrix(Z)
                      mode(Z) <- "numeric"
                      keep_cols <- apply(Z, 2, function(col) {
                        all(is.finite(col)) && var(col, na.rm = TRUE) > 1e-8
                      })
                      Z_clean <- Z[, keep_cols, drop = FALSE]
                      Z_clean <- scale(Z_clean, center = TRUE, scale = FALSE)
                      stopifnot(nrow(Z_clean) == length(y))
                      stopifnot(!anyNA(Z_clean), !any(is.infinite(Z_clean)))

                      return(Z_clean)
                    }
                    Z_train <- prepare_rrblup_matrix(mgeno_scaled[train_ids, ], y_train)
                    model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
                    kept_snps <- colnames(Z_train)
                    Z_test <- mgeno_scaled[test_ids, kept_snps, drop = FALSE]
                    Z_test <- scale(Z_test, center = colMeans(mgeno_scaled[train_ids, kept_snps]), scale = FALSE)
                    pred_rrblup <- as.vector(Z_test %*% model_rrblup$u)

                    # Bayesian-based genomic predictions
                    bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                    if(any(grepl("Bayes",unlist(Additional_models)))){
                      run_model <- function(model_name) {
                        fit <- BGLR(y = Y.masked,
                                    ETA = list(list(X = mgeno_scaled, model = model_name)),
                                    nIter = nIter,
                                    burnIn = burnIn,
                                    verbose = FALSE)
                        pred <- as.numeric(mgeno_scaled %*% fit$ETA[[1]]$b)
                        return(pred)
                      }

                      cl <- makeCluster(ncores)
                      clusterExport(cl, c("BGLR", "Y.masked", "mgeno_scaled", "nIter", "burnIn", "run_model"))
                      preds <- parLapply(cl, bayes_models, run_model)
                      stopCluster(cl)

                      names(preds) <- bayes_models
                      pred_brr <- data.frame(ID = common_ids, Prediction = preds[["BRR"]]); pred_brr <- pred_brr[pred_brr$ID %in% test_ids, ]
                      pred_bayesA <- data.frame(ID = common_ids, Prediction = preds[["BayesA"]]); pred_bayesA <- pred_bayesA[pred_bayesA$ID %in% test_ids, ]
                      pred_bayesB <- data.frame(ID = common_ids, Prediction = preds[["BayesB"]]); pred_bayesB <- pred_bayesB[pred_bayesB$ID %in% test_ids, ]
                      pred_bayesC <- data.frame(ID = common_ids, Prediction = preds[["BayesC"]]); pred_bayesC <- pred_bayesC[pred_bayesC$ID %in% test_ids, ]
                      pred_bayesLasso <- data.frame(ID = common_ids, Prediction = preds[["BL"]]); pred_bayesLasso <- pred_bayesLasso[pred_bayesLasso$ID %in% test_ids, ]

                      if(!("rrBLUP" %in% Additional_models)){
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], BRR = pred_brr[,2], BayesA = pred_bayesA[,2],
                                                      BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      } else {
                        stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup, RKHS = pred_rkhs, BRR = pred_brr[,2],
                                                      BayesA = pred_bayesA[,2], BayesB = pred_bayesB[,2], BayesC = pred_bayesC[,2], BayesLasso = pred_bayesLasso[,2])
                      }
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      if(!("rrBLUP" %in% Additional_models)){
                        fit_stack <- lm(y ~ BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      } else {
                        fit_stack <- lm(y ~ GBLUP + rrBLUP + RKHS + BRR + BayesA + BayesB + BayesC + BayesLasso, data = stack_allmodels)
                      }
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    } else {
                      stack_allmodels <- data.frame(y = Y.raw[rownames(Y.raw) %in% test_ids,2], GBLUP = pred_gblup, rrBLUP = pred_rrblup)
                      stack_allmodels_scaled <- as.data.frame(scale(stack_allmodels[, -1]))
                      stack_allmodels_scaled$y <- stack_allmodels$y
                      stack_allmodels <-  stack_allmodels_scaled
                      fit_stack <- lm(y ~ GBLUP + rrBLUP, data = stack_allmodels)
                      stacked_prediction <- predict(fit_stack, newdata = stack_allmodels)
                      stack_allmodels$StackedPrediction <- stacked_prediction
                      cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                    }
                  }
                }
              }
            }

            #Calculate correlation and store them
            r.GBLUP <- cor(pred_gblup, Y.raw[rownames(Y.raw) %in% test_ids,2])
            storage.GBLUP[rep,1]=r.GBLUP

            if (!is.null(Additional_models)){
              r.rrBLUP <- cor(pred_rrblup, Y.raw[rownames(Y.raw) %in% test_ids,2])
              storage.rrBLUP[rep,1]=r.rrBLUP

              r.RKHS <- cor(pred_rkhs, Y.raw[rownames(Y.raw) %in% test_ids,2])
              storage.RKHS[rep,1]=r.RKHS

              if(any(grepl("Bayes",unlist(Additional_models)))){
                r.BRR <- cor(pred_brr$Prediction, Y.raw[rownames(Y.raw) %in% test_ids,2])
                r.BayesA <- cor(pred_bayesA$Prediction, Y.raw[rownames(Y.raw) %in% test_ids,2])
                r.BayesB <- cor(pred_bayesB$Prediction, Y.raw[rownames(Y.raw) %in% test_ids,2])
                r.BayesC <- cor(pred_bayesC$Prediction, Y.raw[rownames(Y.raw) %in% test_ids,2])
                r.BayesL <- cor(pred_bayesLasso$Prediction, Y.raw[rownames(Y.raw) %in% test_ids,2])
                r.mStacked <- cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                storage.BRR[rep,1]=r.BRR
                storage.BayesA[rep,1]=r.BayesA
                storage.BayesB[rep,1]=r.BayesB
                storage.BayesC[rep,1]=r.BayesC
                storage.BayesL[rep,1]=r.BayesL
                storage.mStacked[rep,1]=r.mStacked
              } else {
                r.BRR <- NA
                r.BayesA <- NA
                r.BayesB <- NA
                r.BayesC <- NA
                r.BayesL <- NA
                r.mStacked <- NA
                r.mStacked <- cor(stack_allmodels$y, stack_allmodels$StackedPrediction)
                storage.BRR[rep,1]=NA
                storage.BayesA[rep,1]=NA
                storage.BayesB[rep,1]=NA
                storage.BayesC[rep,1]=NA
                storage.BayesL[rep,1]=NA
                storage.mStacked[rep,1]=r.mStacked
              }
            } else {
              r.BRR <- NA
              r.BayesA <- NA
              r.BayesB <- NA
              r.BayesC <- NA
              r.BayesL <- NA
              r.mStacked <- NA
              storage.rrBLUP[rep,1]=NA
              storage.RKHS[rep,1]=NA
              storage.BRR[rep,1]=NA
              storage.BayesA[rep,1]=NA
              storage.BayesB[rep,1]=NA
              storage.BayesC[rep,1]=NA
              storage.BayesL[rep,1]=NA
              storage.mStacked[rep,1]=NA
            }

            if(!exists("pop_data")){pop_data <- data.frame(matrix(nrow = 2, ncol = 0))}
            if (gwas_Gpred == TRUE) {
              if (metagenome_covariate == TRUE) {
                time <- timestamp()
                repi <- sprintf("%04d", rep)
                rep_iteration <- paste("completed rep-",repi,"\t", round(r.GBLUP,digits=2),"\t",round(r.rrBLUP,digits=2),"\t",round(r.RKHS,digits=2),"\t",round(r.BRR,digits=2),"\t",round(r.BayesA,digits=2),"\t",
                                       round(r.BayesB,digits=2),"\t",round(r.BayesC,digits=2),"\t",round(r.BayesL,digits=2),"\t",round(r.mStacked,digits=2),"\t", nrow(Y.raw),"\t",(ncol(pop_data)),"\t",time,sep="")
                write(rep_iteration,file=paste(trait,"_progress_gwas_",gene_model,"_metag.txt",sep=""),append=TRUE)
              } else {
                time <- timestamp()
                repi <- sprintf("%04d", rep)
                rep_iteration <- paste("completed rep-",repi,"\t", round(r.GBLUP,digits=2),"\t",round(r.rrBLUP,digits=2),"\t",round(r.RKHS,digits=2),"\t",round(r.BRR,digits=2),"\t",round(r.BayesA,digits=2),"\t",
                                       round(r.BayesB,digits=2),"\t",round(r.BayesC,digits=2),"\t",round(r.BayesL,digits=2),"\t", round(r.mStacked,digits=2),"\t", nrow(Y.raw),"\t",(ncol(pop_data)),"\t",time,sep="")
                write(rep_iteration,file=paste(trait,"_progress_gwas_",gene_model,".txt",sep=""),append=TRUE)
              }
            }
            if (gwas_Gpred == FALSE) {
              if (metagenome_covariate == TRUE) {
                time <- timestamp()
                repi <- sprintf("%04d", rep)
                rep_iteration <- paste("completed rep-",repi,"\t", round(r.GBLUP,digits=2),"\t",round(r.rrBLUP,digits=2),"\t",round(r.RKHS,digits=2),"\t",round(r.BRR,digits=2),"\t",round(r.BayesA,digits=2),"\t",
                                       round(r.BayesB,digits=2),"\t",round(r.BayesC,digits=2),"\t",round(r.BayesL,digits=2),"\t", round(r.mStacked,digits=2),"\t", nrow(Y.raw),"\t",ncol(pop_data),"\t",time,sep="")
                write(rep_iteration,file=paste(trait,"_progress_",gene_model,"_metag.txt",sep=""),append=TRUE)
              } else {
                time <- timestamp()
                repi <- sprintf("%04d", rep)
                rep_iteration <- paste("completed rep-",repi,"\t",round(r.GBLUP,digits=2),"\t",round(r.rrBLUP,digits=2),"\t",round(r.RKHS,digits=2),"\t",round(r.BRR,digits=2),"\t",round(r.BayesA,digits=2),"\t",
                                       round(r.BayesB,digits=2),"\t",round(r.BayesC,digits=2),"\t",round(r.BayesL,digits=2),"\t", round(r.mStacked,digits=2),"\t", nrow(Y.raw),"\t",ncol(pop_data),"\t",time,sep="")
                write(rep_iteration,file=paste(trait,"_progress_",gene_model,".txt",sep=""),append=TRUE)
              }
            }
            gc()
          }#End of for (rep in 1:t)
          stderror <- function(x) sd(x)/sqrt(length(x))
          storage=cbind(storage.GBLUP,storage.rrBLUP,storage.RKHS,storage.BRR,storage.BayesA,storage.BayesB,storage.BayesC,storage.BayesL,storage.mStacked)
          colnames(storage)=c("PA_GBLUP","PA_rrBLUP","PA_RKHS","PA_BRR","PA_BayesA","PA_BayesB","PA_BayesC","PA_BayesL","PA_mStacked")
          write.table(storage, "GAPIT.Cross.Validation.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
          Final_summary <- data.frame(matrix(nrow = 0, ncol = 6))
          colnames(Final_summary) <- c("Traits","Median","Mean","StdErr","PA_method")
          PAmethod=1
          for (PAm in c("GBLUP","rrBLUP","RKHS","BRR","BayesA","BayesB","BayesC","BayesL","mStacked")){
            median <-c(outdir,median(storage[,PAmethod], na.rm = TRUE),
                       mean(storage[,PAmethod], na.rm = TRUE), stderror(storage[,PAmethod]), PAm)
            median <- as.data.frame(t(median))
            colnames(median) <- c("Traits","Median","Mean","StdErr","PA_method")
            Final_summary <- rbind(Final_summary,median)
            PAmethod=PAmethod+1
          }
          write.table(Final_summary, paste("../",GP_run_title,"/",outdir,"gp_summary_stats.txt",sep=""), col.names=TRUE, row.names=FALSE, quote = FALSE, sep = "\t", append=FALSE)
          setwd("../")
          file.rename(paste0(getwd(),"/",outdir),paste0(getwd(),"/",GP_run_title,"/",outdir))
        } else {setwd("../")}
      }
    }
  }
}
