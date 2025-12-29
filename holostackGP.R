#!/usr/bin/env Rscript

holostackGP <- function(
    wdir = "./",
    projname="proj_GP",
    MTME=FALSE,
    phenofile=NULL,
    genofile=NULL,
    metagenomefile=NULL,
    ploidy=2,
    traits=NULL,
    covariate=NULL,
    # "metagenomic or microbiome", genomic","holobiont", "metagenomic+genomic or microbiome+genomic"
    kernel=NULL,
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    gwas_pred=FALSE,
    subsample_markers=NULL,
    topHits=100,
    gene_model="Full",           # "Additive", "Dominance", "metagenome or microbiome", "Full"
    R_libpath=NULL
) {
  load_packages <- function(pkgs, lib = Sys.getenv("R_LIBS_USER")) {

    if (nzchar(lib)) {
      .libPaths(c(lib, .libPaths()))
    }

    for (p in pkgs) {

      if (!requireNamespace(p, quietly = TRUE)) {
        message("Installing missing package: ", p)

        tryCatch(
          install.packages(
            p,
            repos = "https://cloud.r-project.org",
            dependencies = TRUE
          ),
          error = function(e) {
            stop("Failed to install package '", p, "': ", conditionMessage(e))
          }
        )
      }

      ok <- suppressPackageStartupMessages(
        require(p, character.only = TRUE, quietly = TRUE)
      )

      if (!ok) {
        stop("Package '", p, "' is installed but failed to load")
      }
    }

    invisible(TRUE)
  }
  pkgs <- c("data.table", "dplyr", "mice", "AGHmatrix", "vegan", "compositions", "ggcorrplot", "glmnet",
            "nlme", "lsmeans", "agricolae", "parallel", "doParallel", "foreach", "rrBLUP", "BGLR", "GWASpoly")



  #############################################################################################################################################################################
  # Specify parameters
  ####################
  # Step 1: Specify covariate if required and import data.
  setwd(wdir)
  GP_run_title <- projname
  MTME <- MTME
  traits <- traits
  if (is.null(covariate)) {
    covariate <- NULL
  } else {
    covariate <- covariate
  }
  if(tolower(kernel) == "metagenomic") {kernel <- "gBLUP"}
  if(tolower(kernel) == "genomic") {kernel <- "GBLUP"}
  if (tolower(kernel) %in% c("holobiont", "metagenomic+genomic")) {kernel <- "gGBLUP"}
  if (tolower(gene_model) %in% c("metagenome", "microbiome")) {kernel <- "gBLUP"}
  if(kernel == "GBLUP"){metagenomefile <- NULL}
  myY <- read.table(phenofile, head = TRUE, sep="\t", check.names=FALSE)
  myG <- NULL
  if (!is.null(genofile) && !is.na(genofile) && nzchar(genofile) && tolower(genofile) != "null") {
    myG <- read.table(genofile, header = FALSE, sep = "\t", check.names = FALSE)
  }
  metagenome_data <- metagenomefile
  gp_model <- kernel
  gwas_Gpred <- gwas_pred
  gwas_model <- "MLM"
  Additional_models <- c("rrBLUP","BRR","BayesA","BayesB","BayesC","BayesianLasso","RHKS")     # "rrBLUP","BRR","BayesA","BayesB","BayesC","BayesianLasso" or NULL

  number_reps <- as.numeric(CVrep)
  nfold_CV <- as.numeric(k_fold)
  maf_threshold <- as.numeric(maf)
  perc_missing <- as.numeric(geno_missing_rate)
  if (is.null(subsample_markers)) {
    subsample_markers <- NULL
  } else {
    subsample_markers <- as.numeric(subsample_markers)
  }
  gene_models <- gene_model
  select_gwasGPmodel <- NULL                                       #c("2-dom-ref","3-dom-alt","3-dom-ref")
  ploidy_levels <- as.numeric(ploidy)

  weight_by <- "scores"   # use effects, scores (i.e. -log10(pvalues) or pvalues.
  ntop_hits <- 100
  mtraits <- NULL
  ncores <- 5
  LOCO <- FALSE
  drop_ind <- ""

  #metagenome-based parameters
  mincorr <- 0                                 # minimum correlation coefficient (trait vs subsets of metagenome)
  maxcorr <- 0.8                                 # maximum correlation coefficient (trait vs subsets of metagenome)
  pvalue <- 1
  clr_transform <- TRUE
  pheno_zero_inflated <- FALSE
  metag_zero_inflated <- TRUE
  impute_zeroinflation_metagcov <- FALSE
  corr_coeff <- NULL
  entire_metagenome <- TRUE
  crosstrait <- NULL
  options(warn=0)


  delete_short_files <- function(path = ".") {
    files <- list.files(path, full.names = TRUE, recursive = TRUE)
    for (f in files) {
      if (file.info(f)$isdir) next  # skip directories

      n_lines <- length(readLines(f, n = 2))
      if (n_lines < 2) {
        message("Deleting short file: ", f)
        file.remove(f)
      }
    }
  }
  delete_empty_dirs <- function(path = ".") {
    dirs <- list.dirs(path, full.names = TRUE, recursive = TRUE)

    # sort longest path first so we delete deepest directories first
    dirs <- dirs[order(nchar(dirs), decreasing = TRUE)]

    for (d in dirs) {
      if (length(list.files(d, recursive = FALSE)) == 0) {
        message("Deleting empty directory: ", d)
        unlink(d, recursive = TRUE)
      }
    }
  }
  delete_short_files("./")
  delete_empty_dirs("./")



  for (trait in c(traits)) {
    for (ploidy in c(ploidy_levels)) {
      for (gene_model in c(gene_models)) {
        model_selection <- TRUE
        if (gene_model == "full" || gene_model == "FULL" ){ gene_model <- "Full"}
        if (gene_model == "all" || gene_model == "ALL" ){ gene_model <- "All"}
        if (gene_model == "additive" || gene_model == "ADDITIVE" ){ gene_model <- "Additive"}
        if (gene_model == "dominance" || gene_model == "DOMINANCE" ){ gene_model <- "Dominance"}
        if(gp_model == "GBLUP"){metagenome_covariate <- FALSE; metagenome_data <- NULL}
        if(gp_model == "gGBLUP"){metagenome_covariate <- TRUE; metag_method <- "Aitchison"}
        if(gp_model == "gBLUP"){gene_model <- "metagenome"; metag_method <- "Aitchison"; metagenome_covariate <- FALSE}
        metag_pca <- TRUE
        dir.create(GP_run_title, showWarnings=FALSE, recursive=TRUE)
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
        if (!is.null(myG)){
          dG <- myG; colnames(dG) <- dG[1,]; dG <- dG[-1,]
          if("pvalue" %in% colnames(dG)){dG <- subset(dG, select=-c(pvalue))}
        } else {dG <- NULL}
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
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","RKHS","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_gwas_",gene_model,"_metag.txt",sep=""),append=TRUE)
            } else {
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","RKHS","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_gwas_",gene_model,".txt",sep=""),append=TRUE)
            }
          }
          if (gwas_Gpred == FALSE) {
            if (metagenome_covariate == TRUE) {
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","RKHS","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_",gene_model,"_metag.txt",sep=""),append=TRUE)
            } else {
              rep_iteration <- paste("completed repplication","\t","GBLUP","\t","rrBLUP","\t","RKHS","\t","BRR","\t","BayesA","\t","BayesB","\t","BayesC","\t","BayesL","\t","mStacked","\t",
                                     "#_samples","\t","#_markers","\t","time_stamp",sep="")
              write(rep_iteration,file=paste(trait,"_progress_",gene_model,".txt",sep=""),append=TRUE)
            }
          }

          top_percentile <- 0.95; low_percentile <- 0.05
          topP <- subset(Y.raw, Y.raw[,2] >= quantile(Y.raw[,2],top_percentile)); topP$class <- "high_value"
          lowP <- subset(Y.raw, Y.raw[,2] <= quantile(Y.raw[,2],low_percentile)); lowP$class <- "low_value"
          Precision_in <- NULL; Precision_out <- data.frame(matrix(nrow = 0, ncol = 6))
          colnames(Precision_out) <- c("Trait","Rep","Precision_topP","topP_n","Precision_lowP","lowP_n")

          if (gp_model == "gBLUP" || gp_model == "gGBLUP") {
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
              metag_clr <- as.matrix(compositions::clr(metag + 1e-6))
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
          if (gp_model == "GBLUP" || gp_model == "gGBLUP") {
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
            if (is.numeric(subsample_markers)) {
              set.seed(123)
              pop_data <- dplyr::sample_n(pop_data, subsample_markers)
            }
            pop_data <- as.matrix(t(pop_data))
            #Computing the full-autopolyploid matrix based on Slater 2016 (Eq. 8 and 9)
            if (gene_model == "Additive"){
              G_matrix <- AGHmatrix::Gmatrix(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA,
                                             maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                             pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                             ratio.check = FALSE)
              myKI <- normalize_kinmat(as.matrix(G_matrix))
            }
            if (gene_model == "Dominance"){
              if (ploidy == 2){Gmethod <- "Vitezica"}
              if (ploidy > 2){Gmethod <- "Slater"}
              G_matrix <- AGHmatrix::Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
                                             maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                             pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                             ratio.check = FALSE)
              myKI <- normalize_kinmat(as.matrix(G_matrix))
            }
            if (gene_model == "Full" || gene_model == "All" ){
              G_matrix <- AGHmatrix::Gmatrix(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA,
                                             maf = maf_threshold, thresh.missing = 1, verify.posdef = FALSE, ploidy = ploidy,
                                             pseudo.diploid = FALSE, integer = TRUE, ratio = FALSE, impute.method = "mode",
                                             ratio.check = FALSE)
              myKI.Add <- normalize_kinmat(as.matrix(G_matrix))
              if (ploidy == 2){Gmethod <- "Vitezica"}
              if (ploidy > 2){Gmethod <- "Slater"}
              G_matrix <- AGHmatrix::Gmatrix(SNPmatrix = pop_data, method = Gmethod, missingValue = NA,
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
            imp <- mice::mice(Y.cov, m = 5, method = "pmm", seed = 123)
            Y.cov <- mice::complete(imp, 1)   # first imputed dataset
            Y.cov <- Y.cov[,c("Taxa",covariate)]
          } else {Y.cov <- NULL}

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
          #BGLR parameters
          nIter <- 15000           # number of iterations: Bayesian methods using BGLR package
          burnIn <- 5000          # number of burnin iterations: Bayesian methods using BGLR package
          verboseBGLR <- FALSE


          # format genomic/metagenomic relationship matrix
          if(gp_model == "GBLUP" || gp_model == "gGBLUP"){
            if (gene_model == "Full" || gene_model == "All"){
              myKI.A <- t(myKI.Add)
              myKI.A <- as.matrix(myKI.A); myKI.A <- apply(myKI.A, 2, as.numeric); storage.mode(myKI.A) <- "numeric"
              rownames(myKI.A) <- colnames(myKI.A)
              common_ids <- intersect(colnames(myKI.A), rownames(Y.raw))
              myKI.A <- myKI.A[common_ids, common_ids]

              myKI.D <- t(myKI.Dom)
              myKI.D <- as.matrix(myKI.D); myKI.D <- apply(myKI.D, 2, as.numeric); storage.mode(myKI.D) <- "numeric"
              rownames(myKI.D) <- colnames(myKI.D)
              common_ids <- intersect(colnames(myKI.D), rownames(Y.raw))
              myKI.D <- myKI.D[common_ids, common_ids]
            } else {
              myKIx <- t(myKI)
              myKIx <- as.matrix(myKIx); myKIx <- apply(myKIx, 2, as.numeric); storage.mode(myKIx) <- "numeric"
              rownames(myKIx) <- colnames(myKIx)
              common_ids <- intersect(colnames(myKIx), rownames(Y.raw))
              myKIx <- myKIx[common_ids, common_ids]
            }
          }
          if(gp_model == "gBLUP" || gp_model == "gGBLUP"){
            metagKIx <- t(metagKI); colnames(metagKIx) <- metagKIx[1,]; metagKIx <- metagKIx[-1,]; rownames(metagKIx) <- colnames(metagKIx)
            metagKIx <- as.matrix(metagKIx); metagKIx <- apply(metagKIx, 2, as.numeric); storage.mode(metagKIx) <- "numeric"
            rownames(metagKIx) <- colnames(metagKIx)
            common_ids <- intersect(colnames(metagKIx), rownames(Y.raw))
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
          Y.raw <- Y.raw[rownames(Y.raw) %in% common_ids, , drop = FALSE]

          # format data to produce scaled genomic/metagenomic data matrix
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
          Y.raw <- as.numeric(Y.raw[[1]])

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
          }

          for(rep in seq(from=1, to=t, by=1)){
            rep <- as.integer(rep)
            # create empty dataframe for out-of-fold (OOF) predictions
            if (MTME == TRUE){
              print("run MTME a false to generate predictions, which you can use for stacking in a separate R script")
            } else {
              if(gp_model == "GBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  col_names <- c("GBLUP.A", "GBLUP.D", "GBLUP.AxD", "GBLUP.AxA", "GBLUP.DxD")
                  pred_gblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_gblup_OOF) <- col_names
                  col_names <- c("rrBLUP.A", "rrBLUP.D")
                  pred_rrblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rrblup_OOF) <- col_names
                  col_names <- c("RKHS.A", "RKHS.D", "RKHS.AxD", "RKHS.AxA", "RKHS.DxD")
                  pred_rkhs_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rkhs_OOF) <- col_names
                  col_names <- c("BRR.A", "BRR.D", "BayesA.A", "BayesA.D", "BayesB.A", "BayesB.D", "BayesC.A", "BayesC.D", "BL.A", "BL.D")
                  pred_bayes_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_bayes_OOF) <- col_names
                } else {
                  col_names <- c("GBLUP")
                  pred_gblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_gblup_OOF) <- col_names
                  col_names <- c("rrBLUP")
                  pred_rrblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rrblup_OOF) <- col_names
                  col_names <- c("RKHS")
                  pred_rkhs_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rkhs_OOF) <- col_names
                  col_names <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
                  pred_bayes_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_bayes_OOF) <- col_names
                }
              }
              if(gp_model == "gBLUP"){
                col_names <- c("GBLUP")
                pred_gblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_gblup_OOF) <- col_names
                col_names <- c("rrBLUP")
                pred_rrblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rrblup_OOF) <- col_names
                col_names <- c("RKHS")
                pred_rkhs_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rkhs_OOF) <- col_names
                col_names <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
                pred_bayes_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_bayes_OOF) <- col_names
              }
              if(gp_model == "gGBLUP"){
                if (gene_model == "Full" || gene_model == "All"){
                  col_names <- c("GBLUP.A","GBLUP.D","GBLUP.AxD","GBLUP.AxA","GBLUP.DxD","GBLUP.AxAxM","GBLUP.DxDxM","GBLUP.AxDxM","GBLUP.M")
                  pred_gblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_gblup_OOF) <- col_names
                  col_names <- c("rrBLUP.A", "rrBLUP.D", "rrBLUP.M")
                  pred_rrblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rrblup_OOF) <- col_names
                  col_names <- c("RKHS.A","RKHS.D","RKHS.AxD","RKHS.AxA","RKHS.DxD","RKHS.AxAxM","RKHS.DxDxM","RKHS.AxDxM","RKHS.M")
                  pred_rkhs_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rkhs_OOF) <- col_names
                  col_names <- c("BRR.A", "BRR.D", "BRR.M", "BayesA.A", "BayesA.D", "BayesA.M","BayesB.A", "BayesB.D", "BayesB.M",
                                 "BayesC.A", "BayesC.D", "BayesC.DM", "BL.A", "BL.D", "BL.M")
                  pred_bayes_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_bayes_OOF) <- col_names
                } else {
                  col_names <- c("GBLUP", "GBLUP.M")
                  pred_gblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_gblup_OOF) <- col_names
                  col_names <- c("rrBLUP", "rrBLUP.M")
                  pred_rrblup_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rrblup_OOF) <- col_names
                  col_names <- c("RKHS", "RKHS.M")
                  pred_rkhs_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_rkhs_OOF) <- col_names
                  col_names <- c("BRR", "BRR.M", "BayesA", "BayesA.M", "BayesB", "BayesB.M", "BayesC", "BayesC.M", "BL", "BL.M")
                  pred_bayes_OOF <- data.frame(matrix(ncol = length(col_names), nrow = 0)); colnames(pred_bayes_OOF) <- col_names
                }
              }
            }
            parse_Xcov <- function(Y.cov, id_col = "Taxa") {
              if (is.null(Y.cov)) return(NULL)
              if (!id_col %in% colnames(Y.cov))
                stop("Y.cov must contain column: ", id_col)
              covar_cols <- setdiff(colnames(Y.cov), id_col)
              if (length(covar_cols) == 0)
                return(NULL)
              Xcov <- as.matrix(Y.cov[, covar_cols, drop = FALSE])
              rownames(Xcov) <- Y.cov[[id_col]]
              # Drop zero-variance covariates globally
              keep <- apply(Xcov, 2, function(z) sd(z, na.rm = TRUE) > 0)
              Xcov <- Xcov[, keep, drop = FALSE]
              if (ncol(Xcov) == 0) return(NULL)
              Xcov
            }
            Xcov_all <- parse_Xcov(Y.cov)
            split_Xcov <- function(Xcov, train_ids, test_ids) {
              if (is.null(Xcov))
                return(list(X_train = NULL, X_test = NULL))
              X_train <- Xcov[train_ids, , drop = FALSE]
              X_test  <- Xcov[test_ids,  , drop = FALSE]
              # Center using training data only
              mu <- colMeans(X_train, na.rm = TRUE)
              X_train <- sweep(X_train, 2, mu, "-")
              X_test  <- sweep(X_test,  2, mu, "-")
              list(
                X_train = X_train,
                X_test  = X_test
              )
            }

            ## (1) Create folds per replicate
            set.seed(1000+rep)
            fold_id <- sample(rep(1:nfold_CV, length.out = n))
            names(fold_id) <- Y.raw$Taxa
            gc()
            ## (2) Storage for OOF predictions

            # Computing predictions
            for(kfold in 1:nfold_CV){
              sample.missing=sample(1:n,n.missing)
              if(n.missing>0){
                test_ids  <- Y.raw$Taxa[ fold_id == kfold ]
                train_ids <- Y.raw$Taxa[ fold_id != kfold ]

                # Subset phenotype
                Y_train <- Y.raw[Y.raw$Taxa %in% train_ids, trait]
                Y_test  <- Y.raw[Y.raw$Taxa %in% test_ids, trait]

                # Subset covariates
                cov_split <- split_Xcov(Xcov_all, train_ids, test_ids)

                ## ============================================================
                ## Subset relationship matrices (n × n)
                ## ============================================================
                if (gp_model %in% c("GBLUP", "gGBLUP")) {
                  if (gene_model %in% c("Full", "All")) {
                    myKI.A_train <- myKI.A[train_ids, train_ids, drop = FALSE]
                    myKI.A_test  <- myKI.A[test_ids,  train_ids, drop = FALSE]
                    myKI.D_train <- myKI.D[train_ids, train_ids, drop = FALSE]
                    myKI.D_test  <- myKI.D[test_ids,  train_ids, drop = FALSE]
                  } else {
                    myKIx_train <- myKIx[train_ids, train_ids, drop = FALSE]
                    myKIx_test  <- myKIx[test_ids,  train_ids, drop = FALSE]
                  }
                }
                ## Metagenomic kernel (gBLUP / gGBLUP)
                if (gp_model %in% c("gBLUP", "gGBLUP")) {
                  metagKIx_train <- metagKIx[train_ids, train_ids, drop = FALSE]
                  metagKIx_test  <- metagKIx[test_ids,  train_ids, drop = FALSE]
                }

                ## ============================================================
                ## Subset genotype matrices (n × p)
                ## ============================================================
                if (gp_model %in% c("GBLUP", "gGBLUP")) {
                  if (gene_model %in% c("Full", "All")) {
                    geno.A_scaled_train <- geno.A_scaled[train_ids, , drop = FALSE]
                    geno.A_scaled_test  <- geno.A_scaled[test_ids,  , drop = FALSE]
                    geno.D_scaled_train <- geno.D_scaled[train_ids, , drop = FALSE]
                    geno.D_scaled_test  <- geno.D_scaled[test_ids,  , drop = FALSE]
                  } else {
                    geno_scaled_train <- geno_scaled[train_ids, , drop = FALSE]
                    geno_scaled_test  <- geno_scaled[test_ids,  , drop = FALSE]
                  }
                }
                ## Metagenomic genotype matrix (n × p)
                if (gp_model %in% c("gBLUP", "gGBLUP")) {
                  mgeno_scaled_train <- mgeno_scaled[train_ids, , drop = FALSE]
                  mgeno_scaled_test  <- mgeno_scaled[test_ids,  , drop = FALSE]
                }

                ## ============================================================
                ## Subset epistatic kernels (n × n)
                ## ============================================================
                kernels_train <- lapply(kernels, function(K) {
                  K[train_ids, train_ids, drop = FALSE]
                })
                kernels_test <- lapply(kernels, function(K) {
                  K[test_ids, train_ids, drop = FALSE]
                })

                ## ------------------------------------------------------------
                ## Sanity checks
                stopifnot(
                  nrow(myKI.A_train) == length(train_ids),
                  ncol(myKI.A_train) == length(train_ids),
                  nrow(myKI.A_test)  == length(test_ids),
                  ncol(myKI.A_test)  == length(train_ids)
                )
                stopifnot(
                  nrow(geno.A_scaled_train) == length(train_ids),
                  nrow(geno.A_scaled_test)  == length(test_ids)
                )

                if (!is.null(covariate)){
                  cov_split <- split_Xcov(Xcov_all, train_ids, test_ids)
                  X_train <- cov_split$X_train
                  X_test  <- cov_split$X_test
                } else {X_train <- NULL; X_test <- NULL}

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

              method_names <- c("GBLUP", "rrBLUP_markers", "RKHS_BGLR", "BRR", "BayesA", "BayesB", "BayesC", "BL")
              # --- Covariate preparation ---
              prepare_covariates <- function(Xcov) {

                # Return NULL immediately if no covariates
                if (is.null(Xcov)) {
                  return(NULL)
                }

                # Coerce to data.frame safely
                if (is.list(Xcov) && !is.data.frame(Xcov)) {
                  Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                } else {
                  Xcov_df <- as.data.frame(Xcov)
                }

                # Empty data.frame guard
                if (ncol(Xcov_df) == 0) {
                  return(NULL)
                }

                # Drop all-NA or constant columns
                keep <- vapply(Xcov_df, function(x) {
                  !(all(is.na(x)) || length(unique(na.omit(x))) <= 1)
                }, logical(1))
                Xcov_df <- Xcov_df[, keep, drop = FALSE]

                if (ncol(Xcov_df) == 0) {
                  return(NULL)
                }

                # Build design matrix (no intercept)
                Xcov_mat <- model.matrix(~ . - 1, data = Xcov_df)

                # Drop duplicate columns after factor expansion
                dup <- duplicated(colnames(Xcov_mat))
                if (any(dup)) {
                  Xcov_mat <- Xcov_mat[, !dup, drop = FALSE]
                }

                # Drop near-zero variance columns
                nzv <- apply(Xcov_mat, 2, function(x) var(x, na.rm = TRUE) < 1e-8)
                if (any(nzv)) {
                  message(
                    "Dropping ", sum(nzv), " near-zero variance covariates: ",
                    paste(colnames(Xcov_mat)[nzv], collapse = ", ")
                  )
                  Xcov_mat <- Xcov_mat[, !nzv, drop = FALSE]
                }

                # Enforce full rank by dropping collinear columns
                if (ncol(Xcov_mat) > 1) {
                  qrX <- qr(Xcov_mat)
                  if (qrX$rank < ncol(Xcov_mat)) {
                    drop_idx <- setdiff(seq_len(ncol(Xcov_mat)), qrX$pivot[seq_len(qrX$rank)])
                    message(
                      "Dropping collinear covariates: ",
                      paste(colnames(Xcov_mat)[drop_idx], collapse = ", ")
                    )
                    keep_idx <- qrX$pivot[seq_len(qrX$rank)]
                    Xcov_mat <- Xcov_mat[, keep_idx, drop = FALSE]
                  }
                }

                # Final check: return NULL if no valid covariates left
                if (ncol(Xcov_mat) == 0) {
                  return(NULL)
                }

                return(Xcov_mat)
              }

              # Genomic predictions
              if (MTME == TRUE){
                print("run MTME a false to generate predictions, which you can use for stacking in a separate R script")
              } else {
                if(gp_model == "GBLUP"){
                  if (gene_model == "Full" || gene_model == "All"){
                    # GBLUP  with rrBLUP package
                    pred_list <- list()
                    for (kernel_name in names(kernels_train)) {
                      K_train <- kernels_train[[kernel_name]]   # train × train
                      K <- K_train
                      K_test  <- kernels_test[[kernel_name]]    # test  × train
                      covTraits <- character(0)
                      if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                      gebv_list <- list()
                      if (length(covTraits) > 0) {
                        for (covt in covTraits){
                          y <- X_train[[covt]]
                          # remove NA individuals (rrBLUP::mixed.solve handles but safer)
                          ok <- !is.na(y)
                          sol <- rrBLUP::mixed.solve(y = y[ok], K = K[ok, ok, drop = FALSE])
                          # sol$u is GEBV vector for individuals with y; bring back into full n vector
                          u_full <- rep(NA, nrow(X_train))
                          u_full[which(ok)] <- sol$u[rownames(K)[ok]]   # name-align
                          gebv_list[[covt]] <- u_full
                        }
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
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

                        # if rownames missing, assign from Y.ttrain
                        if (is.null(rownames(Xcov_mat))) {
                          rownames(Xcov_mat) <- rownames(Y.ttrain)
                        }

                        # if nrow mismatch, force alignment by merging
                        if (nrow(Xcov_mat) != nrow(Y.ttrain)) {
                          Xcov_mat <- Xcov_mat[match(rownames(Y.ttrain), rownames(Xcov_mat)), , drop = FALSE]
                        }

                        stopifnot(nrow(Xcov_mat) == nrow(Y.ttrain))
                      }

                      model_gblup <- rrBLUP::mixed.solve(y = Y.ttrain[,1], K = K, X=Xcov_mat)
                      pred_gblup <- as.numeric(K_test %*% model_gblup$u)
                      pred_gblup <- data.frame(pred = pred_gblup, row.names = test_ids, stringsAsFactors = FALSE)
                      pred_list[[kernel_name]] <- pred_gblup
                    }

                    pred_gblup_all <- do.call(cbind, lapply(pred_list, function(x) x[, 1, drop = FALSE]))
                    colnames(pred_gblup_all) <- paste0("GBLUP.", names(pred_list))
                    rownames(pred_gblup_all) <- test_ids
                    pred_gblup_OOF <- rbind(pred_gblup_OOF,pred_gblup_all)

                    if (!is.null(Additional_models)){
                      # rrBLUP marker effects model: K_A
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      geno.A_scaled_train <- as.matrix(geno.A_scaled_train)
                      mode(geno.A_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(geno.A_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      Z_test <- geno.A_scaled_test[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(geno.A_scaled_test[test_ids, kept_snps]), scale = FALSE)
                      pred_rrblup_A <- as.vector(Z_test %*% model_rrblup$u)

                      # rrBLUP marker effects model: K_D
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      geno.D_scaled_train <- as.matrix(geno.D_scaled_train)
                      mode(geno.D_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(geno.D_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      Z_test <- geno.D_scaled_test[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(geno.D_scaled_test[test_ids, kept_snps]), scale = FALSE)
                      pred_rrblup_D <- as.vector(Z_test %*% model_rrblup$u)

                      pred_rrblup_all <- data.frame(rrBLUP.A  = pred_rrblup_A, rrBLUP.D  = pred_rrblup_D)
                      rownames(pred_rrblup_all) <- test_ids
                      pred_rrblup_OOF <- rbind(pred_rrblup_OOF,pred_rrblup_all)

                      # RHKs with BGLR package
                      pred_list <- list()
                      for (kernel_name in names(kernels_train)) {
                        K_train <- kernels_train[[kernel_name]]   # train × train
                        K <- K_train
                        K_test  <- kernels_test[[kernel_name]]    # test  × train
                        covTraits <- character(0)
                        if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                        gebv_list <- list()
                        for(covt in covTraits){
                          y <- X_train[[covt]]
                          ok <- !is.na(y)
                          library(BGLR)
                          fm <- BGLR(y = y[ok],  ETA = list(list(K=K[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                          # Create full-length vector aligned with original y
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.ttrain <- as.data.frame(Y_train)
                        if (length(covTraits) > 0) {
                          for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                        }
                        Xcov_mat <- as.matrix(Y.ttrain[, -1, drop = FALSE])
                        train_idx <- which(!is.na(Y.ttrain[,1]))
                        test_idx  <- which(is.na(Y.ttrain[,1]))
                        y_train <- Y.ttrain[train_idx, 1]
                        K_train <- K[train_idx, train_idx]
                        X_train <- Xcov_mat[train_idx, , drop = FALSE]
                        library(BGLR)
                        if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                        ## Step 4: build ETA conditionally
                        ETA <- list(list(K = K_train, model = "RKHS"))
                        if (!is.null(X_train) && ncol(X_train) > 0) {
                          X_train <- Xcov_mat[train_idx, , drop = FALSE]
                          ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                        }
                        ## Step 5: fit final model
                        fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        u_train <- as.numeric(fit$ETA[[1]]$u)
                        names(u_train) <- train_ids
                        pred_test <- as.vector(K_test %*% u_train)
                        names(pred_test) <- test_ids
                        pred_list[[kernel_name]] <- as.numeric(pred_test)
                      }

                      pred_rkhs_all <- do.call(cbind, pred_list)
                      colnames(pred_rkhs_all) <- paste0("RKHS.", names(pred_list))
                      rownames(pred_rkhs_all) <- test_ids
                      pred_rkhs_OOF <- rbind(pred_rkhs_OOF,pred_rkhs_all)

                      # Bayesian-based genomic predictions
                      bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                      gc()
                      if(any(grepl("Bayes",unlist(Additional_models)))){
                        # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                        run_independent_bayes <- function(model_name, Y_train, X_train = NULL, geno.A_scaled_train, geno.D_scaled_train, nIter, burnIn) {
                          suppressPackageStartupMessages(library(BGLR))
                          covTraits <- character(0)
                          if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                          # ---- Additive step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = geno.A_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(geno.A_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = geno.A_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit_A <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all_A <- rep(NA, nrow(Y.ttrain))
                          yHat_all_A[train_idx] <- fit_A$yHat
                          b_markersA <- fit_A$ETA[[1]]$b
                          b_fixedA <- if (length(fit_A$ETA) >= 2) fit_A$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markersA)) geno.A_scaled_test %*% b_markersA else 0
                          fixed_part  <- if (!is.null(b_fixedA))   X_test %*% b_fixedA else 0
                          pred_A <- as.numeric(marker_part + fixed_part)
                          # ---- Dominance step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              library(BGLR)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = geno.D_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(geno.D_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = geno.D_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit_D <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all_D <- rep(NA, nrow(Y.ttrain))
                          yHat_all_D[train_idx] <- fit_D$yHat
                          b_markersD <- fit_D$ETA[[1]]$b
                          b_fixedD <- if (length(fit_D$ETA) >= 2) fit_D$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markersD)) geno.D_scaled_test %*% b_markersD else 0
                          fixed_part  <- if (!is.null(b_fixedD))   X_test %*% b_fixedD else 0
                          pred_D <- as.numeric(marker_part + fixed_part)
                          # pred_D <- as.numeric(geno.D_scaled_train[test_idx, , drop = FALSE] %*% b_markersD + X_test %*% b_fixedD)
                          # ---- Return predictions for test set ----
                          return(data.frame(A = pred_A, D = pred_D, row.names = rownames(Y.ttrain)[test_idx]))
                        }

                        # Wrapper to run all models in parallel
                        run_parallel_stack <- function(Y_train, X_train, covariate = NULL, geno.A_scaled_train, geno.D_scaled_train, nIter, burnIn, n.cores = ncores) {
                          # If covariate is NULL, make sure X_train is also NULL
                          if (is.null(covariate)) {X_train <- NULL}
                          # If X_train exists but has zero columns, set it to NULL
                          if (!is.null(X_train) && ncol(X_train) == 0) {X_train <- NULL}
                          cl <- parallel::makeCluster(n.cores, type = "PSOCK")
                          on.exit(parallel::stopCluster(cl), add = TRUE)
                          parallel::clusterSetRNGStream(cl, iseed = 12345)
                          parallel::clusterEvalQ(cl, {
                            suppressPackageStartupMessages(library(BGLR))
                            NULL
                          })
                          parallel::clusterExport(cl, varlist = c("Y_train", "X_train", "geno.A_scaled_train", "geno.D_scaled_train", "nIter", "burnIn", "run_independent_bayes", "prepare_covariates"), envir = environment())
                          preds_list <- parallel::parLapply(cl, bayes_models, function(model) {
                            run_independent_bayes(model, Y_train = Y_train, X_train = X_train, geno.A_scaled_train = geno.A_scaled_train, geno.D_scaled_train = geno.D_scaled_train, nIter = nIter, burnIn = burnIn)
                          })
                          names(preds_list) <- bayes_models
                          return(preds_list)
                        }
                        # Run all Bayesian models in parallel
                        preds_stack <- run_parallel_stack(Y_train = Y_train, X_train = X_train, geno.A_scaled_train = geno.A_scaled_train, geno.D_scaled_train = geno.D_scaled_train, nIter = nIter, burnIn = burnIn, n.cores = ncores)
                        pred_bayes_all <- as.data.frame(preds_stack)
                        rownames(pred_bayes_all) <- test_ids
                        pred_bayes_OOF <- rbind( pred_bayes_OOF,pred_bayes_all)
                      }
                    }

                  } else {
                    # GBLUP with rrBLUP package
                    covTraits <- character(0)
                    {if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                      gebv_list <- list()
                      if (length(covTraits) > 0) {
                        for (covt in covTraits){
                          y <- X_train[[covt]]
                          # remove NA individuals (rrBLUP::mixed.solve handles but safer)
                          ok <- !is.na(y)
                          sol <- rrBLUP::mixed.solve(y = y[ok], K = myKIx_train[ok, ok, drop = FALSE])
                          # sol$u is GEBV vector for individuals with y; bring back into full n vector
                          u_full <- rep(NA, nrow(X_train))
                          u_full[which(ok)] <- sol$u[rownames(myKIx_train)[ok]]   # name-align
                          gebv_list[[covt]] <- u_full
                        }
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
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

                        # if rownames missing, assign from Y.ttrain
                        if (is.null(rownames(Xcov_mat))) {
                          rownames(Xcov_mat) <- rownames(Y.ttrain)
                        }

                        # if nrow mismatch, force alignment by merging
                        if (nrow(Xcov_mat) != nrow(Y.ttrain)) {
                          Xcov_mat <- Xcov_mat[match(rownames(Y.ttrain), rownames(Xcov_mat)), , drop = FALSE]
                        }

                        stopifnot(nrow(Xcov_mat) == nrow(Y.ttrain))
                      }}

                    model_gblup <- rrBLUP::mixed.solve(y = Y.ttrain[,1], K = myKIx_train, X=Xcov_mat)
                    pred_gblup <- as.numeric(K_test %*% model_gblup$u)
                    pred_gblup <- data.frame(pred = pred_gblup, row.names = test_ids, stringsAsFactors = FALSE)
                    colnames(pred_gblup) <- "GBLUP"
                    pred_gblup_OOF <- rbind(pred_gblup_OOF,pred_gblup)

                    if (!is.null(Additional_models)){
                      # rrBLUP marker effects model
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      # Ensure genotype matrix is numeric
                      geno_scaled_train <- as.matrix(geno_scaled_train)
                      mode(geno_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(geno_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      # Align test set to training SNPs (safe check)
                      Z_test <- geno_scaled_train[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test,
                                      center = colMeans(geno_scaled_train[train_ids, kept_snps]),
                                      scale = FALSE)
                      pred_rrblup <- as.data.frame(as.vector(Z_test %*% model_rrblup$u))
                      rownames(pred_rrblup) <- test_ids
                      colnames(pred_rrblup) <- "rrBLUP"
                      pred_rrblup_OOF <- rbind(pred_rrblup_OOF,pred_rrblup)

                      # RHKs with BGLR package
                      covTraits <- character(0)
                      if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                      gebv_list <- list()
                      for(covt in covTraits){
                        y <- X_train[[covt]]
                        ok <- !is.na(y)
                        library(BGLR)
                        fm <- BGLR(y = y[ok],  ETA = list(list(K=myKIx_train[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                        # Create full-length vector aligned with original y
                        gebv_full <- rep(NA, length(y))
                        gebv_full[ok] <- as.numeric(b)
                        gebv_list[[covt]] <- gebv_full
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- as.matrix(Y.ttrain[, -1, drop = FALSE])
                      train_idx <- which(!is.na(Y.ttrain[,1]))
                      test_idx  <- which(is.na(Y.ttrain[,1]))
                      y_train <- Y.ttrain[train_idx, 1]
                      K_train <- myKIx_train[train_idx, train_idx]
                      X_train <- Xcov_mat[train_idx, , drop = FALSE]
                      library(BGLR)
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                      ## Step 4: build ETA conditionally
                      ETA <- list(list(K = K_train, model = "RKHS"))
                      if (!is.null(X_train) && ncol(X_train) > 0) {
                        X_train <- Xcov_mat[train_idx, , drop = FALSE]
                        ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                      }
                      ## Step 5: fit final model
                      fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      u_train <- as.numeric(fit$ETA[[1]]$u)
                      K_test_train <- myKIx_train[test_idx, train_idx]
                      pred_rkhs <- as.data.frame(as.numeric( K_test_train %*% u_train))
                      rownames(pred_rkhs) <- test_ids
                      colnames(pred_rkhs) <- "RKHS"
                      pred_rkhs_OOF <- rbind(pred_rkhs_OOF,pred_rkhs)

                      # Bayesian-based genomic predictions
                      bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                      gc()
                      if(any(grepl("Bayes",unlist(Additional_models)))){
                        # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                        run_independent_bayes <- function(model_name, Y_train, X_train = NULL, geno_scaled_train, nIter, burnIn) {
                          suppressPackageStartupMessages(library(BGLR))
                          covTraits <- character(0)
                          if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                          # ---- Additive step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = geno_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(geno_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = geno_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all <- rep(NA, nrow(Y.ttrain))
                          yHat_all[train_idx] <- fit$yHat
                          b_markers <- fit$ETA[[1]]$b
                          b_fixed <- if (length(fit$ETA) >= 2) fit$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markers)) geno_scaled_test %*% b_markers else 0
                          fixed_part  <- if (!is.null(b_fixed))   X_test %*% b_fixed else 0
                          pred <- as.numeric(marker_part + fixed_part)
                          return(data.frame(pred = pred, row.names = rownames(Y.ttrain)[test_idx]))
                        }

                        # Wrapper to run all models in parallel
                        run_parallel_stack <- function(Y_train, X_train, covariate = NULL, geno_scaled_train, nIter, burnIn, n.cores = ncores) {
                          # If covariate is NULL, make sure X_train is also NULL
                          if (is.null(covariate)) {X_train <- NULL}
                          # If X_train exists but has zero columns, set it to NULL
                          if (!is.null(X_train) && ncol(X_train) == 0) {X_train <- NULL}
                          cl <- parallel::makeCluster(n.cores, type = "PSOCK")
                          on.exit(parallel::stopCluster(cl), add = TRUE)
                          parallel::clusterSetRNGStream(cl, iseed = 12345)
                          parallel::clusterEvalQ(cl, {
                            suppressPackageStartupMessages(library(BGLR))
                            NULL
                          })
                          parallel::clusterExport(cl, varlist = c("Y_train", "X_train", "geno_scaled_train", "nIter", "burnIn", "run_independent_bayes", "prepare_covariates"), envir = environment())
                          preds_list <- parallel::parLapply(cl, bayes_models, function(model) {
                            run_independent_bayes(model, Y_train = Y_train, X_train = X_train, geno_scaled_train = geno_scaled_train, nIter = nIter, burnIn = burnIn)
                          })
                          names(preds_list) <- bayes_models
                          return(preds_list)
                        }

                        # Run all Bayesian models in parallel
                        preds <- run_parallel_stack(Y_train = Y_train, X_train = X_train, geno_scaled_train = geno_scaled_train, nIter = nIter, burnIn = burnIn, n.cores = ncores)
                        pred_bayes_all <- as.data.frame(preds)
                        rownames(pred_bayes_all) <- test_ids
                        pred_bayes_OOF <- rbind(pred_bayes_OOF,pred_bayes_all)
                      }
                    }
                  }
                }
                if(gp_model == "gBLUP"){
                    # GBLUP  with rrBLUP package
                    pred_list <- list()
                    for (kernel_name in names(kernels_train)) {
                      K_train <- kernels_train[[kernel_name]]   # train × train
                      K <- K_train
                      K_test  <- kernels_test[[kernel_name]]    # test  × train
                      covTraits <- character(0)
                      if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                      gebv_list <- list()
                      if (length(covTraits) > 0) {
                        for (covt in covTraits){
                          y <- X_train[[covt]]
                          # remove NA individuals (rrBLUP::mixed.solve handles but safer)
                          ok <- !is.na(y)
                          sol <- rrBLUP::mixed.solve(y = y[ok], K = K[ok, ok, drop = FALSE])
                          # sol$u is GEBV vector for individuals with y; bring back into full n vector
                          u_full <- rep(NA, nrow(X_train))
                          u_full[which(ok)] <- sol$u[rownames(K)[ok]]   # name-align
                          gebv_list[[covt]] <- u_full
                        }
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
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

                        # if rownames missing, assign from Y.ttrain
                        if (is.null(rownames(Xcov_mat))) {
                          rownames(Xcov_mat) <- rownames(Y.ttrain)
                        }

                        # if nrow mismatch, force alignment by merging
                        if (nrow(Xcov_mat) != nrow(Y.ttrain)) {
                          Xcov_mat <- Xcov_mat[match(rownames(Y.ttrain), rownames(Xcov_mat)), , drop = FALSE]
                        }

                        stopifnot(nrow(Xcov_mat) == nrow(Y.ttrain))
                      }

                      model_gblup <- rrBLUP::mixed.solve(y = Y.ttrain[,1], K = K, X=Xcov_mat)
                      pred_gblup <- as.numeric(K_test %*% model_gblup$u)
                      pred_gblup <- data.frame(pred = pred_gblup, row.names = test_ids, stringsAsFactors = FALSE)
                      pred_list[[kernel_name]] <- pred_gblup
                    }
                    pred_gblup_all <- do.call(cbind, lapply(pred_list, function(x) x[, 1, drop = FALSE]))
                    colnames(pred_gblup_all) <- names(pred_list)
                    rownames(pred_gblup) <- "GBLUP"
                    rownames(pred_gblup_all) <- test_ids
                    pred_gblup_OOF <- rbind(pred_gblup_OOF,pred_gblup_all)

                    if (!is.null(Additional_models)){
                      # rrBLUP marker effects model: K_M
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      mgeno_scaled_train <- as.matrix(mgeno_scaled_train)
                      mode(mgeno_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(mgeno_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      Z_test <- mgeno_scaled_test[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(mgeno_scaled_test[test_ids, kept_snps]), scale = FALSE)
                      pred_rrblup <- as.data.frame(as.vector(Z_test %*% model_rrblup$u))
                      colnames(pred_rrblup) <- "rrBLUP"
                      rownames(pred_rrblup) <- test_ids
                      pred_rrblup_OOF <- rbind(pred_rrblup_OOF,pred_rrblup)

                      # RHKs with BGLR package
                      pred_list <- list()
                      for (kernel_name in names(kernels_train)) {
                        K_train <- kernels_train[[kernel_name]]   # train × train
                        K <- K_train
                        K_test  <- kernels_test[[kernel_name]]    # test  × train
                        covTraits <- character(0)
                        if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                        gebv_list <- list()
                        for(covt in covTraits){
                          y <- X_train[[covt]]
                          ok <- !is.na(y)
                          library(BGLR)
                          fm <- BGLR(y = y[ok],  ETA = list(list(K=K[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                          # Create full-length vector aligned with original y
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.ttrain <- as.data.frame(Y_train)
                        if (length(covTraits) > 0) {
                          for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                        }
                        Xcov_mat <- as.matrix(Y.ttrain[, -1, drop = FALSE])
                        train_idx <- which(!is.na(Y.ttrain[,1]))
                        test_idx  <- which(is.na(Y.ttrain[,1]))
                        y_train <- Y.ttrain[train_idx, 1]
                        K_train <- K[train_idx, train_idx]
                        X_train <- Xcov_mat[train_idx, , drop = FALSE]
                        library(BGLR)
                        if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                        ## Step 4: build ETA conditionally
                        ETA <- list(list(K = K_train, model = "RKHS"))
                        if (!is.null(X_train) && ncol(X_train) > 0) {
                          X_train <- Xcov_mat[train_idx, , drop = FALSE]
                          ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                        }
                        ## Step 5: fit final model
                        fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        u_train <- as.numeric(fit$ETA[[1]]$u)
                        names(u_train) <- train_ids
                        pred_test <- as.vector(K_test %*% u_train)
                        names(pred_test) <- test_ids
                        pred_list[[kernel_name]] <- as.numeric(pred_test)
                      }
                      pred_rkhs_all <- do.call(cbind, pred_list)
                      rownames(pred_rkhs_all) <- rownames(Y.ttrain)[test_idx]
                      colnames(pred_rkhs_all) <- "RKHS"
                      rownames(pred_rkhs_all) <- test_ids
                      pred_rkhs_OOF <- rbind(pred_rkhs_OOF,pred_rkhs_all)

                      # Bayesian-based genomic predictions
                      bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                      gc()
                      if(any(grepl("Bayes",unlist(Additional_models)))){
                        # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                        run_independent_bayes <- function(model_name, Y_train, X_train = NULL, mgeno_scaled_train, nIter, burnIn) {
                          suppressPackageStartupMessages(library(BGLR))
                          covTraits <- character(0)
                          if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                          # ---- Metagenome step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = mgeno_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(mgeno_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = mgeno_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit_M <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all_M <- rep(NA, nrow(Y.ttrain))
                          yHat_all_M[train_idx] <- fit_M$yHat
                          b_markersM <- fit_M$ETA[[1]]$b
                          b_fixedM <- if (length(fit_M$ETA) >= 2) fit_M$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markersM)) mgeno_scaled_test %*% b_markersM else 0
                          fixed_part  <- if (!is.null(b_fixedM))   X_test %*% b_fixedM else 0
                          pred_M <- as.numeric(marker_part + fixed_part)

                          # ---- Return predictions for test set ----
                          return(data.frame(pred_M = pred_M, row.names = rownames(Y.ttrain)[test_idx]))
                        }

                        # Wrapper to run all models in parallel
                        run_parallel_stack <- function(Y_train, X_train, covariate = NULL, mgeno_scaled_train, nIter, burnIn, n.cores = ncores) {
                          # If covariate is NULL, make sure X_train is also NULL
                          if (is.null(covariate)) {X_train <- NULL}
                          # If X_train exists but has zero columns, set it to NULL
                          if (!is.null(X_train) && ncol(X_train) == 0) {X_train <- NULL}
                          cl <- parallel::makeCluster(n.cores, type = "PSOCK")
                          on.exit(parallel::stopCluster(cl), add = TRUE)
                          parallel::clusterSetRNGStream(cl, iseed = 12345)
                          parallel::clusterEvalQ(cl, {
                            suppressPackageStartupMessages(library(BGLR))
                            NULL
                          })
                          parallel::clusterExport(cl, varlist = c("Y_train", "X_train", "mgeno_scaled_train", "nIter", "burnIn", "run_independent_bayes", "prepare_covariates"), envir = environment())
                          preds_list <- parallel::parLapply(cl, bayes_models, function(model) {
                            run_independent_bayes(model, Y_train = Y_train, X_train = X_train,mgeno_scaled_train = mgeno_scaled_train, nIter = nIter, burnIn = burnIn)
                          })
                          names(preds_list) <- bayes_models
                          return(preds_list)
                        }

                        # Run all Bayesian models in parallel
                        preds <- run_parallel_stack(Y_train = Y_train, X_train = X_train, mgeno_scaled_train = mgeno_scaled_train, nIter = nIter, burnIn = burnIn, n.cores = ncores)
                        pred_bayes_all <- as.data.frame(preds)
                        rownames(pred_bayes_all) <- test_ids
                        pred_bayes_OOF <- rbind(pred_bayes_OOF,pred_bayes_all)
                      }
                    }
                  }
                if(gp_model == "gGBLUP"){
                  if (gene_model == "Full" || gene_model == "All"){
                    # GBLUP  with rrBLUP package
                    pred_list <- list()
                    for (kernel_name in names(kernels_train)) {
                      K_train <- kernels_train[[kernel_name]]   # train × train
                      K <- K_train
                      K_test  <- kernels_test[[kernel_name]]    # test  × train
                      covTraits <- character(0)
                      if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                      gebv_list <- list()
                      if (length(covTraits) > 0) {
                        for (covt in covTraits){
                          y <- X_train[[covt]]
                          # remove NA individuals (rrBLUP::mixed.solve handles but safer)
                          ok <- !is.na(y)
                          sol <- rrBLUP::mixed.solve(y = y[ok], K = K[ok, ok, drop = FALSE])
                          # sol$u is GEBV vector for individuals with y; bring back into full n vector
                          u_full <- rep(NA, nrow(X_train))
                          u_full[which(ok)] <- sol$u[rownames(K)[ok]]   # name-align
                          gebv_list[[covt]] <- u_full
                        }
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
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

                        # if rownames missing, assign from Y.ttrain
                        if (is.null(rownames(Xcov_mat))) {
                          rownames(Xcov_mat) <- rownames(Y.ttrain)
                        }

                        # if nrow mismatch, force alignment by merging
                        if (nrow(Xcov_mat) != nrow(Y.ttrain)) {
                          Xcov_mat <- Xcov_mat[match(rownames(Y.ttrain), rownames(Xcov_mat)), , drop = FALSE]
                        }

                        stopifnot(nrow(Xcov_mat) == nrow(Y.ttrain))
                      }

                      model_gblup <- rrBLUP::mixed.solve(y = Y.ttrain[,1], K = K, X=Xcov_mat)
                      pred_gblup <- as.numeric(K_test %*% model_gblup$u)
                      pred_gblup <- data.frame(pred = pred_gblup, row.names = test_ids, stringsAsFactors = FALSE)
                      pred_list[[kernel_name]] <- pred_gblup
                    }
                    pred_gblup_all <- do.call(cbind, lapply(pred_list, function(x) x[, 1, drop = FALSE]))
                    colnames(pred_gblup_all) <- names(pred_list)
                    colnames(pred_gblup_all) <- paste0("GBLUP.", names(pred_list))
                    rownames(pred_gblup_all) <- test_ids
                    pred_gblup_OOF <- rbind(pred_gblup_OOF,pred_gblup_all)

                    if (!is.null(Additional_models)){
                      # rrBLUP marker effects model: K_A
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      geno.A_scaled_train <- as.matrix(geno.A_scaled_train)
                      mode(geno.A_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(geno.A_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      Z_test <- geno.A_scaled_test[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(geno.A_scaled_test[test_ids, kept_snps]), scale = FALSE)
                      pred_rrblup_A <- as.vector(Z_test %*% model_rrblup$u)

                      # rrBLUP marker effects model: K_D
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      geno.D_scaled_train <- as.matrix(geno.D_scaled_train)
                      mode(geno.D_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(geno.D_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      Z_test <- geno.D_scaled_test[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(geno.D_scaled_test[test_ids, kept_snps]), scale = FALSE)
                      pred_rrblup_D <- as.vector(Z_test %*% model_rrblup$u)

                      # rrBLUP marker effects model: K_M
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      mgeno_scaled_train <- as.matrix(mgeno_scaled_train)
                      mode(mgeno_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(mgeno_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      Z_test <- mgeno_scaled_train[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(mgeno_scaled_test[test_ids, kept_snps]), scale = FALSE)
                      pred_rrblup_M <- as.vector(Z_test %*% model_rrblup$u)

                      pred_rrblup_all <- data.frame(rrBLUP.A  = pred_rrblup_A, rrBLUP.D  = pred_rrblup_D, rrBLUP.M  = pred_rrblup_M)
                      rownames(pred_rrblup_all) <- test_ids
                      pred_rrblup_OOF <- rbind(pred_rrblup_OOF,pred_rrblup_all)


                      # RHKs with BGLR package
                      pred_list <- list()
                      for (kernel_name in names(kernels_train)) {
                        K_train <- kernels_train[[kernel_name]]   # train × train
                        K <- K_train
                        K_test  <- kernels_test[[kernel_name]]    # test  × train
                        covTraits <- character(0)
                        if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                        gebv_list <- list()
                        for(covt in covTraits){
                          y <- X_train[[covt]]
                          ok <- !is.na(y)
                          library(BGLR)
                          fm <- BGLR(y = y[ok],  ETA = list(list(K=K[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                          # Create full-length vector aligned with original y
                          gebv_full <- rep(NA, length(y))
                          gebv_full[ok] <- as.numeric(b)
                          gebv_list[[covt]] <- gebv_full
                        }
                        Y.ttrain <- as.data.frame(Y_train)
                        if (length(covTraits) > 0) {
                          for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                        }
                        Xcov_mat <- as.matrix(Y.ttrain[, -1, drop = FALSE])
                        train_idx <- which(!is.na(Y.ttrain[,1]))
                        test_idx  <- which(is.na(Y.ttrain[,1]))
                        y_train <- Y.ttrain[train_idx, 1]
                        K_train <- K[train_idx, train_idx]
                        X_train <- Xcov_mat[train_idx, , drop = FALSE]
                        library(BGLR)
                        if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                        ## Step 4: build ETA conditionally
                        ETA <- list(list(K = K_train, model = "RKHS"))
                        if (!is.null(X_train) && ncol(X_train) > 0) {
                          X_train <- Xcov_mat[train_idx, , drop = FALSE]
                          ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                        }
                        ## Step 5: fit final model
                        fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        u_train <- as.numeric(fit$ETA[[1]]$u)
                        names(u_train) <- train_ids
                        pred_test <- as.vector(K_test %*% u_train)
                        pred_list[[kernel_name]] <- as.numeric(pred_test)
                      }
                      pred_rkhs_all <- do.call(cbind, pred_list)
                      rownames(pred_rkhs_all) <- rownames(Y.ttrain)[test_idx]
                      colnames(pred_rkhs_all) <- paste0("RKHS.", names(pred_list))
                      rownames(pred_rkhs_all) <- test_ids
                      pred_rkhs_OOF <- rbind(pred_rkhs_OOF,pred_rkhs_all)


                      # Bayesian-based genomic predictions
                      bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                      gc()
                      if(any(grepl("Bayes",unlist(Additional_models)))){
                        # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                        run_independent_bayes <- function(model_name, Y_train, X_train = NULL, geno.A_scaled_train, geno.D_scaled_train, mgeno_scaled_train, nIter, burnIn) {
                          suppressPackageStartupMessages(library(BGLR))
                          covTraits <- character(0)
                          if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                          # ---- Additive step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = geno.A_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(geno.A_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = geno.A_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit_A <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all_A <- rep(NA, nrow(Y.ttrain))
                          yHat_all_A[train_idx] <- fit_A$yHat
                          b_markersA <- fit_A$ETA[[1]]$b
                          b_fixedA <- if (length(fit_A$ETA) >= 2) fit_A$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markersA)) geno.A_scaled_test %*% b_markersA else 0
                          fixed_part  <- if (!is.null(b_fixedA))   X_test %*% b_fixedA else 0
                          pred_A <- as.numeric(marker_part + fixed_part)
                          # pred_A <- as.numeric(geno.A_scaled_train[test_idx, , drop = FALSE] %*% b_markersA + X_test %*% b_fixedA)
                          # ---- Dominance step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              library(BGLR)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = geno.D_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(geno.D_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = geno.D_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit_D <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all_D <- rep(NA, nrow(Y.ttrain))
                          yHat_all_D[train_idx] <- fit_D$yHat
                          b_markersD <- fit_D$ETA[[1]]$b
                          b_fixedD <- if (length(fit_D$ETA) >= 2) fit_D$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markersD)) geno.D_scaled_test %*% b_markersD else 0
                          fixed_part  <- if (!is.null(b_fixedD))   X_test %*% b_fixedD else 0
                          pred_D <- as.numeric(marker_part + fixed_part)
                          # pred_D <- as.numeric(geno.D_scaled_train[test_idx, , drop = FALSE] %*% b_markersD + X_test %*% b_fixedD)
                          # ---- Metagenome step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              library(BGLR)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = mgeno_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(mgeno_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = mgeno_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit_M <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all_M <- rep(NA, nrow(Y.ttrain))
                          yHat_all_M[train_idx] <- fit_M$yHat
                          b_markersM <- fit_M$ETA[[1]]$b
                          b_fixedM <- if (length(fit_M$ETA) >= 2) fit_M$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markersM)) mgeno_scaled_test %*% b_markersM else 0
                          fixed_part  <- if (!is.null(b_fixedM))   X_test %*% b_fixedM else 0
                          pred_M <- as.numeric(marker_part + fixed_part)
                          # ---- Return predictions for test set ----
                          return(data.frame(A = pred_A, D = pred_D, M = pred_M, row.names = rownames(Y.ttrain)[test_idx]))
                        }

                        # Wrapper to run all models in parallel
                        run_parallel_stack <- function(Y_train, X_train, covariate = NULL, geno.A_scaled_train, geno.D_scaled_train, mgeno_scaled_train, nIter, burnIn, n.cores = ncores) {
                          # If covariate is NULL, make sure X_train is also NULL
                          if (is.null(covariate)) {X_train <- NULL}
                          # If X_train exists but has zero columns, set it to NULL
                          if (!is.null(X_train) && ncol(X_train) == 0) {X_train <- NULL}
                          cl <- parallel::makeCluster(n.cores, type = "PSOCK")
                          on.exit(parallel::stopCluster(cl), add = TRUE)
                          parallel::clusterSetRNGStream(cl, iseed = 12345)
                          parallel::clusterEvalQ(cl, {
                            suppressPackageStartupMessages(library(BGLR))
                            NULL
                          })
                          parallel::clusterExport(cl, varlist = c("Y_train", "X_train", "geno.A_scaled_train", "geno.D_scaled_train", "mgeno_scaled_train", "nIter", "burnIn", "run_independent_bayes", "prepare_covariates"), envir = environment())
                          preds_list <- parallel::parLapply(cl, bayes_models, function(model) {
                            run_independent_bayes(model, Y_train = Y_train, X_train = X_train,geno.A_scaled_train = geno.A_scaled_train, geno.D_scaled_train = geno.D_scaled_train, mgeno_scaled_train = mgeno_scaled_train, nIter = nIter, burnIn = burnIn)
                          })
                          names(preds_list) <- bayes_models
                          return(preds_list)
                        }

                        # Run all Bayesian models in parallel
                        preds_stack <- run_parallel_stack(Y_train = Y_train, X_train = X_train, geno.A_scaled_train = geno.A_scaled_train, geno.D_scaled_train = geno.D_scaled_train, mgeno_scaled_train = mgeno_scaled_train, nIter = nIter, burnIn = burnIn, n.cores = ncores)
                        pred_bayes_all <- as.data.frame(preds_stack)
                        rownames(pred_bayes_all) <- test_ids
                        pred_bayes_OOF <- rbind( pred_bayes_OOF,pred_bayes_all)
                      }

                    }
                  } else {
                    # GBLUP with rrBLUP package
                    {
                      covTraits <- character(0)
                      if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                      gebv_list <- list()
                      if (length(covTraits) > 0) {
                        for (covt in covTraits){
                          y <- X_train[[covt]]
                          # remove NA individuals (rrBLUP::mixed.solve handles but safer)
                          ok <- !is.na(y)
                          sol <- rrBLUP::mixed.solve(y = y[ok], K = myKIx_train[ok, ok, drop = FALSE])
                          # sol$u is GEBV vector for individuals with y; bring back into full n vector
                          u_full <- rep(NA, nrow(X_train))
                          u_full[which(ok)] <- sol$u[rownames(myKIx_train)[ok]]   # name-align
                          gebv_list[[covt]] <- u_full
                        }
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov <- Y.ttrain[,-1]
                      model_gblup_g <- rrBLUP::mixed.solve(y = Y.ttrain[,1], K = myKIx_train, X=Xcov)
                      pred_gblup_g <- as.numeric(K_test %*% model_gblup_g$u)

                      gebv_list <- list()
                      if (length(covTraits) > 0) {
                        for (covt in covTraits){
                          y <- X_train[[covt]]
                          # remove NA individuals (rrBLUP::mixed.solve handles but safer)
                          ok <- !is.na(y)
                          sol <- rrBLUP::mixed.solve(y = y[ok], K =metagKIx_train[ok, ok, drop = FALSE])
                          # sol$u is GEBV vector for individuals with y; bring back into full n vector
                          u_full <- rep(NA, nrow(X_train))
                          u_full[which(ok)] <- sol$u[rownames(metagKIx_train)[ok]]   # name-align
                          gebv_list[[covt]] <- u_full
                        }
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
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

                        # if rownames missing, assign from Y.ttrain
                        if (is.null(rownames(Xcov_mat))) {
                          rownames(Xcov_mat) <- rownames(Y.ttrain)
                        }

                        # if nrow mismatch, force alignment by merging
                        if (nrow(Xcov_mat) != nrow(Y.ttrain)) {
                          Xcov_mat <- Xcov_mat[match(rownames(Y.ttrain), rownames(Xcov_mat)), , drop = FALSE]
                        }

                        stopifnot(nrow(Xcov_mat) == nrow(Y.ttrain))
                      }
                      model_gblup_m <- rrBLUP::mixed.solve(y = Y.ttrain[,1], K = metagKIx_train, X=Xcov_mat)
                      pred_gblup_m <- as.numeric(K_test %*% model_gblup_m$u)
                    }
                    pred_gblup_all <- data.frame(GBLUP.G  = pred_gblup_g, GBLUP.M  = pred_gblup_m)
                    rownames(pred_gblup_all) <- test_ids
                    pred_gblup_OOF <- rbind(pred_gblup_OOF,pred_gblup_all)


                    if (!is.null(Additional_models)){
                      # metagenomic RHKs with BGLR package
                      covTraits <- character(0)
                      if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                      gebv_list <- list()
                      for(covt in covTraits){
                        y <- X_train[[covt]]
                        ok <- !is.na(y)
                        library(BGLR)
                        fm <- BGLR(y = y[ok],  ETA = list(list(K=metagKIx_train[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                        # Create full-length vector aligned with original y
                        gebv_full <- rep(NA, length(y))
                        gebv_full[ok] <- as.numeric(b)
                        gebv_list[[covt]] <- gebv_full
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- as.matrix(Y.ttrain[, -1, drop = FALSE])
                      train_idx <- which(!is.na(Y.ttrain[,1]))
                      test_idx  <- which(is.na(Y.ttrain[,1]))
                      y_train <- Y.ttrain[train_idx, 1]
                      K_train <- metagKIx_train[train_idx, train_idx]
                      X_train <- Xcov_mat[train_idx, , drop = FALSE]
                      library(BGLR)
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                      ## Step 4: build ETA conditionally
                      ETA <- list(list(K = K_train, model = "RKHS"))
                      if (!is.null(X_train) && ncol(X_train) > 0) {
                        X_train <- Xcov_mat[train_idx, , drop = FALSE]
                        ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                      }
                      ## Step 5: fit final model
                      fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      u_train <- as.numeric(fit$ETA[[1]]$u)
                      pred_rkhs_m <- as.numeric( K_test_train %*% u_train)

                      # genomic RHKs with BGLR package
                      gebv_list <- list()
                      for(covt in covTraits){
                        y <- X_train[[covt]]
                        ok <- !is.na(y)
                        library(BGLR)
                        fm <- BGLR(y = y[ok],  ETA = list(list(K=myKIx_train[ok, ok], model = "RKHS")), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                        b <- fm$ETA[[1]]$u  # 'u' contains the random effects / GEBVs for the n_obs
                        # Create full-length vector aligned with original y
                        gebv_full <- rep(NA, length(y))
                        gebv_full[ok] <- as.numeric(b)
                        gebv_list[[covt]] <- gebv_full
                      }
                      Y.ttrain <- as.data.frame(Y_train)
                      if (length(covTraits) > 0) {
                        for (covt in covTraits) { Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]]) }
                      }
                      Xcov_mat <- as.matrix(Y.ttrain[, -1, drop = FALSE])
                      train_idx <- which(!is.na(Y.ttrain[,1]))
                      test_idx  <- which(is.na(Y.ttrain[,1]))
                      y_train <- Y.ttrain[train_idx, 1]
                      K_train <- myKIx_train[train_idx, train_idx]
                      X_train <- Xcov_mat[train_idx, , drop = FALSE]
                      library(BGLR)
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                      ## Step 4: build ETA conditionally
                      ETA <- list(list(K = K_train, model = "RKHS"))
                      if (!is.null(X_train) && ncol(X_train) > 0) {
                        X_train <- Xcov_mat[train_idx, , drop = FALSE]
                        ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                      }
                      ## Step 5: fit final model
                      fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                      u_train <- as.numeric(fit$ETA[[1]]$u)
                      K_test_train <- myKIx_test[test_idx, train_idx]
                      pred_rkhs_g <-as.numeric( K_test_train %*% u_train)

                      pred_rkhs_all <- data.frame(RKHS.G  = pred_rkhs_g, RKHS.M  = pred_rkhs_m)
                      rownames(pred_rkhs_all) <- test_ids
                      pred_rkhs_OOF <- rbind(pred_rkhs_OOF,pred_rkhs_all)

                      # metagenomic rrBLUP marker effects model
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      geno_scaled_train <- as.matrix(mgeno_scaled_train)
                      mode(mgeno_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(mgeno_scaled_train[train_ids, ], y_train)
                      if (is.list(Xcov) && !is.data.frame(Xcov)) {
                        Xcov_df <- as.data.frame(do.call(cbind, Xcov))
                      } else {
                        Xcov_df <- as.data.frame(Xcov)
                      }
                      Xcov_mat <- model.matrix(~ ., data = Xcov_df)
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X=Xcov_mat)
                      kept_snps <- colnames(Z_train)
                      Z_test <- mgeno_scaled_test[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(mgeno_scaled_test), scale = FALSE)
                      pred_rrblup_m <- as.vector(Z_test %*% model_rrblup$u)

                      # genomic rrBLUP marker effects model
                      train_ids <- which(!is.na(Y_train))
                      y_train <- Y_train[train_ids]
                      geno_scaled_train <- as.matrix(geno_scaled_train)
                      mode(geno_scaled_train) <- "numeric"
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
                      Z_train <- prepare_rrblup_matrix(geno_scaled_train[train_ids, ], y_train)

                      # Usage: safely assign or NULL
                      Xcov_mat <- prepare_covariates(Y.ttrain[, -1, drop = FALSE])
                      if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {
                        Xcov_mat <- NULL
                      }
                      # --- rrBLUP model ---
                      model_rrblup <- rrBLUP::mixed.solve(y = y_train, Z = Z_train, X = Xcov_mat)
                      # --- Prediction ---
                      kept_snps <- colnames(Z_train)
                      Z_test <- geno_scaled_test[test_ids, kept_snps, drop = FALSE]
                      Z_test <- scale(Z_test, center = colMeans(geno_scaled_test[test_ids, kept_snps]), scale = FALSE)
                      pred_rrblup_g <- as.vector(Z_test %*% model_rrblup$u)

                      pred_rrblup_all <- data.frame(RKHS.G  = pred_rrblup_g, RKHS.M  = pred_rrblup_m)
                      rownames(pred_rrblup_all) <- test_ids
                      pred_rrblup_OOF <- rbind(pred_rrblup_OOF,pred_rrblup_all)


                      # Bayesian-based genomic predictions
                      bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")  # BL = Bayesian Lasso
                      gc()
                      if(any(grepl("Bayes",unlist(Additional_models)))){
                        # Function to fit a single Bayesian model for one phenotype and one or more multiple covariates
                        run_independent_bayes <- function(model_name, Y_train, X_train = NULL, geno_scaled_train, mgeno_scaled_train, nIter, burnIn) {
                          suppressPackageStartupMessages(library(BGLR))
                          covTraits <- character(0)
                          if (!is.null(X_train) && ncol(X_train) > 0) {covTraits <- colnames(X_train)}
                          # ---- Genomic step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = geno_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(geno_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = geno_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit <- BGLR(y = y_train, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
                          yHat_all <- rep(NA, nrow(Y.ttrain))
                          yHat_all[train_idx] <- fit$yHat
                          b_markers <- fit$ETA[[1]]$b
                          b_fixed <- if (length(fit$ETA) >= 2) fit$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markers)) geno_scaled_test %*% b_markers else 0
                          fixed_part  <- if (!is.null(b_fixed))   X_test %*% b_fixed else 0
                          pred_g <- as.numeric(marker_part + fixed_part)

                          # ---- Metagenomic step ----
                          gebv_list <- list()
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {
                              y <- X_train[[covt]]
                              ok <- !is.na(y)
                              library(BGLR)
                              fm <- BGLR(y = y[ok], ETA = list(list(X = mgeno_scaled_train[ok, , drop = FALSE], model = model_name)), nIter = nIter, burnIn = burnIn, verbose = FALSE)
                              b <- fm$ETA[[1]]$b
                              gebv_full <- rep(NA, length(y))
                              gebv_full[ok] <- as.numeric(mgeno_scaled_train[ok, , drop = FALSE] %*% b)
                              gebv_list[[covt]] <- gebv_full
                            }
                          }
                          Y.ttrain <- as.data.frame(Y_train)
                          if (length(covTraits) > 0) {
                            for (covt in covTraits) {Y.ttrain[[paste0("gebv_", covt)]] <- as.vector(gebv_list[[covt]])}
                          }
                          train_idx <- which(!is.na(Y.ttrain[,1]))
                          test_idx  <- which(is.na(Y.ttrain[,1]))
                          y_train <- Y.ttrain[train_idx, 1]
                          X_train <- as.matrix(Y.ttrain[train_idx, -1, drop = FALSE])
                          X_test  <- as.matrix(Y.ttrain[test_idx, -1, drop = FALSE])
                          library(BGLR)
                          if (is.null(Xcov_mat) || ncol(Xcov_mat) == 0) {Xcov_mat <- NULL}
                          ## Step 4: build ETA conditionally
                          ETA <- list(list(X = mgeno_scaled_train[train_idx, , drop = FALSE], model = model_name))
                          if (!is.null(X_train) && ncol(X_train) > 0) {
                            X_train <- Xcov_mat[train_idx, , drop = FALSE]
                            ETA <- c(ETA, list(list(X = X_train, model = "FIXED")))
                          }
                          ## Step 5: fit final model
                          fit <- BGLR(y = y_train, ETA = ETA, burnIn = burnIn, verbose = FALSE)
                          yHat_all <- rep(NA, nrow(Y.ttrain))
                          yHat_all[train_idx] <- fit$yHat
                          b_markers <- fit$ETA[[1]]$b
                          b_fixed <- if (length(fit$ETA) >= 2) fit$ETA[[2]]$b else NULL
                          X_test[is.na(X_test)] <- 0
                          marker_part <- if (!is.null(b_markers)) mgeno_scaled_test %*% b_markers else 0
                          fixed_part  <- if (!is.null(b_fixed))   X_test %*% b_fixed else 0
                          pred_m <- as.numeric(marker_part + fixed_part)

                          return(data.frame(G = pred_g, M = pred_m, row.names = rownames(Y.ttrain)[test_idx]))
                        }

                        # Wrapper to run all models in parallel
                        run_parallel_stack <- function(Y_train, X_train, covariate = NULL, geno_scaled_train, mgeno_scaled_train, nIter, burnIn, n.cores = ncores) {
                          # If covariate is NULL, make sure X_train is also NULL
                          if (is.null(covariate)) {X_train <- NULL}
                          # If X_train exists but has zero columns, set it to NULL
                          if (!is.null(X_train) && ncol(X_train) == 0) {X_train <- NULL}
                          cl <- parallel::makeCluster(n.cores, type = "PSOCK")
                          on.exit(parallel::stopCluster(cl), add = TRUE)
                          parallel::clusterSetRNGStream(cl, iseed = 12345)
                          parallel::clusterEvalQ(cl, {
                            suppressPackageStartupMessages(library(BGLR))
                            NULL
                          })
                          parallel::clusterExport(cl, varlist = c("Y_train", "X_train", "geno_scaled_train","mgeno_scaled_train", "nIter", "burnIn", "run_independent_bayes", "prepare_covariates"), envir = environment())
                          preds_list <- parallel::parLapply(cl, bayes_models, function(model) {
                            run_independent_bayes(model, Y_train = Y_train, X_train = X_train, geno_scaled_train = geno_scaled_train, mgeno_scaled_train = mgeno_scaled_train, nIter = nIter, burnIn = burnIn)
                          })
                          names(preds_list) <- bayes_models
                          return(preds_list)
                        }

                        # Run all Bayesian models in parallel
                        preds_stack <- run_parallel_stack(Y_train = Y_train, X_train = X_train, geno_scaled_train = geno_scaled_train, mgeno_scaled_train = mgeno_scaled_train, nIter = nIter, burnIn = burnIn, n.cores = ncores)
                        pred_bayes_all <- as.data.frame(preds_stack)
                        rownames(pred_bayes_all) <- test_ids
                        pred_bayes_OOF <- rbind( pred_bayes_OOF,pred_bayes_all)
                      }
                    }
                    gc()
                  }
                  gc()
                }
                gc()
              }
              gc()
            }
            gc(); rownames(Y.raw) <- Y.raw$Taxa

            # Model stacking
            if (MTME == TRUE){
              print("run MTME a false to generate predictions, which you can use for stacking in a separate R script")
            } else {
              # ---------------------------------------------
              # Two-step Elastic Net stacking for multi-kernel, multi-GP predictions
              # ---------------------------------------------

              stack_predictions_cv <- function(
                trait,
                gp_model = c("gBLUP", "GBLUP", "gGBLUP"),
                Y.raw,
                pred_gblup_OOF = NULL,
                pred_rrblup_OOF = NULL,
                pred_rkhs_OOF  = NULL,
                pred_bayes_OOF = NULL,
                id_col = "Taxa",
                seed = 123
              ) {

                set.seed(seed)
                gp_model <- match.arg(gp_model)

                ## -----------------------------
                ## 1. Prepare phenotype
                ## -----------------------------
                if (!id_col %in% colnames(Y.raw))
                  stop("Y.raw must contain ID column: ", id_col)

                rownames(Y.raw) <- Y.raw[[id_col]]

                if (!trait %in% colnames(Y.raw))
                  stop("Trait not found in Y.raw: ", trait)

                y_all <- Y.raw[[trait]]
                names(y_all) <- rownames(Y.raw)

                ## -----------------------------
                ## 2. Utilities
                ## -----------------------------
                cor_safe <- function(a, b) {
                  ok <- complete.cases(a, b)
                  if (sum(ok) < 5) return(NA_real_)
                  if (sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
                  cor(a[ok], b[ok])
                }

                stage1_stack <- function(X, y) {
                  X <- as.matrix(X)

                  keep <- apply(X, 2, function(z) sd(z, na.rm = TRUE) > 0)
                  X <- X[, keep, drop = FALSE]

                  if (ncol(X) == 0)
                    return(rep(NA_real_, length(y)))

                  if (ncol(X) == 1)
                    return(as.numeric(X[, 1]))

                  ok <- complete.cases(X, y)
                  if (sum(ok) < 5)
                    return(rep(NA_real_, length(y)))

                  fit <- tryCatch(
                    lm(y[ok] ~ X[ok, , drop = FALSE]),
                    error = function(e) NULL
                  )

                  yhat <- rep(NA_real_, length(y))

                  if (!is.null(fit)) {
                    yhat[ok] <- predict(fit, newdata = as.data.frame(X[ok, , drop = FALSE]))
                  }

                  if (all(is.na(yhat[ok]))) {
                    yhat[ok] <- rowMeans(X[ok, , drop = FALSE], na.rm = TRUE)
                  }

                  yhat
                }

                elastic_net_stack <- function(X, y) {
                  X <- as.matrix(X)
                  ok <- complete.cases(X, y)
                  if (sum(ok) < 5) return(rep(NA_real_, length(y)))

                  fit <- tryCatch(
                    glmnet::cv.glmnet(X[ok, , drop = FALSE], y[ok], alpha = 0.5),
                    error = function(e) NULL
                  )

                  yhat <- rep(NA_real_, length(y))
                  if (!is.null(fit))
                    yhat[ok] <- as.numeric(
                      predict(fit, newx = X[ok, , drop = FALSE], s = "lambda.min")
                    )

                  yhat
                }

                ## -----------------------------
                ## 3. Collect all OOF columns
                ## -----------------------------
                oof_all <- list(
                  pred_gblup_OOF,
                  pred_rrblup_OOF,
                  pred_rkhs_OOF,
                  pred_bayes_OOF
                )

                oof_all <- Filter(Negate(is.null), oof_all)
                if (length(oof_all) == 0)
                  stop("No OOF predictions supplied")

                oof_all <- do.call(cbind, oof_all)

                ## -----------------------------
                ## 4. Parse GP method from column names
                ## -----------------------------
                method_map <- sub("\\..*$", "", colnames(oof_all))
                methods <- unique(method_map)

                ## -----------------------------
                ## 5. Align IDs
                ## -----------------------------
                common_ids <- intersect(rownames(oof_all), names(y_all))
                if (length(common_ids) < 5)
                  stop("Insufficient overlap between phenotype and OOF predictions")

                y <- y_all[common_ids]
                oof_all <- oof_all[common_ids, , drop = FALSE]

                ## -----------------------------
                ## 6. Stage-1: within-method stacking
                ## -----------------------------
                method_preds <- list()
                method_cor   <- numeric(length(methods))
                names(method_cor) <- methods

                for (m in methods) {
                  cols <- which(method_map == m)
                  X <- oof_all[, cols, drop = FALSE]

                  yhat <- stage1_stack(X, y)
                  method_preds[[m]] <- yhat
                  method_cor[m] <- cor_safe(y, yhat)
                }

                ## -----------------------------
                ## 7. Stage-2: across-method stacking
                ## -----------------------------
                X2 <- do.call(cbind, method_preds)
                pred_stack <- elastic_net_stack(X2, y)
                pred_stack_cor <- cor_safe(y, pred_stack)

                ## -----------------------------
                ## 8. Return
                ## -----------------------------
                list(
                  pred_stack     = pred_stack,
                  pred_stack_cor = pred_stack_cor,
                  method_preds   = method_preds,
                  method_cor     = method_cor,
                  seed           = seed
                )
              }


              stack_result <- stack_predictions_cv(
                trait = trait,
                gp_model = gp_model,
                Y.raw = Y.raw,
                pred_gblup_OOF = pred_gblup_OOF,
                pred_rrblup_OOF = pred_rrblup_OOF,
                pred_rkhs_OOF = pred_rkhs_OOF,
                pred_bayes_OOF = pred_bayes_OOF,
                alpha = 0.5,
                seed = 123
              )
          }


            # Extract results
            pred_all       <- stack_result$pred_all      # all individual GP method predictions merged
            pred_stack     <- stack_result$pred_stack    # final stacked prediction
            pred_stack_cor <- stack_result$pred_stack_cor
            base_cor       <- stack_result$base_cor      # correlations for base models (GBLUP, rrBLUP, RKHS)
            bayes_cor      <- stack_result$bayes_cor     # correlations for Bayesian models (BRR, BayesA/B/C, BL)

            # ------------------------
            # Store correlations in replication matrices
            # ------------------------
            r.GBLUP <- if ("GBLUP" %in% names(base_cor)) base_cor[["GBLUP"]] else NA
            r.rrBLUP <- if ("rrBLUP" %in% names(base_cor)) base_cor[["rrBLUP"]] else NA
            r.RKHS  <- if ("RKHS" %in% names(base_cor)) base_cor[["RKHS"]] else NA

            r.BRR    <- if (!is.null(bayes_cor) && "BRR" %in% names(bayes_cor)) bayes_cor["BRR"] else NA
            r.BayesA <- if (!is.null(bayes_cor) && "BayesA" %in% names(bayes_cor)) bayes_cor["BayesA"] else NA
            r.BayesB <- if (!is.null(bayes_cor) && "BayesB" %in% names(bayes_cor)) bayes_cor["BayesB"] else NA
            r.BayesC <- if (!is.null(bayes_cor) && "BayesC" %in% names(bayes_cor)) bayes_cor["BayesC"] else NA
            r.BayesL <- if (!is.null(bayes_cor) && "BL" %in% names(bayes_cor)) bayes_cor["BL"] else NA

            r.mStacked <- pred_stack_cor

            # ------------------------
            # Store correlations in replication matrices
            # ------------------------
            storage.GBLUP[rep, 1] <- r.GBLUP
            storage.rrBLUP[rep, 1] <- r.rrBLUP
            storage.RKHS[rep, 1] <- r.RKHS
            storage.BRR[rep, 1] <- r.BRR
            storage.BayesA[rep, 1] <- r.BayesA
            storage.BayesB[rep, 1] <- r.BayesB
            storage.BayesC[rep, 1] <- r.BayesC
            storage.BayesL[rep, 1] <- r.BayesL
            storage.mStacked[rep, 1] <- r.mStacked


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
            outfile_predall <- paste("../",GP_run_title,"/",outdir,"ALL_PREDICTIONS.txt",sep="")
            write.table(pred_all, outfile_predall, col.names = !file.exists(outfile_predall), row.names=FALSE, quote = FALSE, sep = "\t", append=file.exists(outfile_predall))
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
