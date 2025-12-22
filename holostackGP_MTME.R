#!/usr/bin/env Rscript

## ---------------------------
## CONFIGURATION
## ---------------------------
config <- list(
  traits = c("VIR1","VIR2"),
  environments = data.frame(
    EnvFile = c(
      "VIR1_ABI19A_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt",
      "VIR1_ABI19B_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt",
      "VIR1_NGE19B_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt",
      "VIR1_NGE20A_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt",
      "VIR2_ABI19A_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt",
      "VIR2_ABI19B_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt",
      "VIR2_NGE19B_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt",
      "VIR2_NGE20A_metagenome2x_gBLUP_nometag_ALL_PREDICTIONS.txt"
    ),
    Trait = c("VIR1","VIR1","VIR1","VIR1","VIR2","VIR2","VIR2","VIR2"),
    Environment = c("ABI19A","ABI19B","NGE19B","NGE20A","ABI19A","ABI19B","NGE19B","NGE20A"),
    stringsAsFactors = FALSE
  ),
  methods = c("MBLUP","rrBLUP","RKHS","BRR","BayesA","BayesB","BayesC","BL"),
  stacking_methods = c("linear","ridge","lasso","elasticnet","bayesian"),
  prediction_dir = "./",
  cv_rep_col = "rep",
  taxa_col = "Taxa",
  min_common_taxa = 5
)

## ---------------------------
## LIBRARIES
## ---------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(glmnet)
  library(BGLR)
})

## ---------------------------
## UTILITY
## ---------------------------
safe_cor <- function(y, yhat) {
  keep <- complete.cases(y, yhat)
  y <- y[keep]; yhat <- yhat[keep]
  if (length(y) < 3 || sd(y, na.rm=TRUE)==0 || sd(yhat, na.rm=TRUE)==0) return(NA_real_)
  cor(y, yhat)
}

## ---------------------------
## LOAD DATA
## ---------------------------
build_mtme_dataframe <- function(config) {
  message("Building unified MTME dataframe...")
  map_dfr(seq_len(nrow(config$environments)), function(i) {
    env <- config$environments[i, ]
    f <- file.path(config$prediction_dir, env$EnvFile)
    if (!file.exists(f)) stop("Missing file: ", f)
    df <- fread(f)
    if ("Row.names" %in% names(df)) df$Row.names <- NULL

    pheno <- setdiff(names(df), c(config$methods, config$taxa_col, config$cv_rep_col))
    if (length(pheno) != 1) stop("Expected exactly one phenotype in ", env$EnvFile)

    df %>%
      rename(y = all_of(pheno)) %>%
      mutate(Trait = env$Trait,
             Environment = env$Environment,
             ID = .data[[config$taxa_col]])
  })
}

## ---------------------------
## STAGE 1: GP STACKING PER ENVIRONMENT
## ---------------------------
stack_gp_models <- function(df_env, gp_cols, method="ridge") {
  df_env <- as.data.frame(df_env)
  gp_cols_exist <- intersect(gp_cols, names(df_env))
  if (length(gp_cols_exist) < 2) return(list(yhat=rep(NA_real_, nrow(df_env)), ability=NA_real_))

  X <- df_env[, gp_cols_exist, drop=FALSE]
  y <- df_env$y
  if (nrow(X) != length(y)) stop("X and y row counts differ")

  keep <- complete.cases(X, y)
  yhat <- rep(NA_real_, nrow(df_env))

  if (sum(keep) >= 5 && sd(y[keep], na.rm=TRUE) > 0) {
    X_scaled <- scale(X[keep,,drop=FALSE])
    y_sub <- y[keep]

    yhat_sub <- tryCatch({
      if (method=="linear") {
        fit <- lm(y_sub ~ X_scaled)
        as.numeric(cbind(1, X_scaled) %*% coef(fit))
      } else if (method %in% c("ridge","lasso","elasticnet")) {
        alpha <- switch(method, ridge=0, lasso=1, elasticnet=0.5)
        cv <- cv.glmnet(X_scaled, y_sub, alpha=alpha)
        as.numeric(predict(cv, X_scaled, s="lambda.min"))
      } else if (method=="bayesian") {
        fit <- BGLR(y=y_sub, ETA=list(list(X=X_scaled, model="BRR")),
                    nIter=6000, burnIn=1000, verbose=FALSE)
        fit$yHat
      } else stop("Unknown stacking method")
    }, error=function(e) rep(NA_real_, length(y_sub)))

    yhat[keep] <- yhat_sub
  }

  list(yhat=yhat, ability=safe_cor(y[keep], yhat[keep]))
}

## ---------------------------
## STAGE 2: ME STACKING
## ---------------------------
stack_environments <- function(df_rep, stack_method) {
  # select all environment-specific Stage1 columns
  yhat_cols <- grep(paste0("^GPstacked_", stack_method, "_"), names(df_rep), value=TRUE)
  if(length(yhat_cols) < 2) return(NA_real_)

  df_wide <- df_rep %>%
    select(ID, y, all_of(yhat_cols)) %>%
    filter(rowSums(!is.na(select(., all_of(yhat_cols)))) > 0)

  # Compute row-wise mean across available environment predictions
  yhat_matrix <- as.matrix(df_wide[, -c(1,2)])
  yhat_avg <- rowMeans(yhat_matrix, na.rm=TRUE)

  safe_cor(df_wide$y, yhat_avg)
}

## ---------------------------
## STAGE 3: MT STACKING
## ---------------------------
stack_traits <- function(df) {
  df <- df %>% filter(!is.na(Ability))
  if (nrow(df)<2 || length(unique(df$Trait))<2) return(NA_real_)
  df$Trait <- factor(df$Trait)
  fit <- lm(Ability ~ Trait, data=df)
  safe_cor(df$Ability, predict(fit))
}

## ---------------------------
## MAIN MTME PIPELINE
## ---------------------------
run_mtme_pipeline <- function(config) {
  data <- build_mtme_dataframe(config)
  results <- list()
  diagnostics <- list()

  for (trait in unique(data$Trait)) {
    df_trait <- data %>% filter(Trait==trait)
    for (rep_i in unique(df_trait[[config$cv_rep_col]])) {
      df_rep <- df_trait %>% filter(.data[[config$cv_rep_col]]==rep_i)

      for (stack_method in config$stacking_methods) {
        # Stage1: environment-specific predictions
        for (env in unique(df_rep$Environment)) {
          df_env <- df_rep %>% filter(Environment==env)
          s1 <- stack_gp_models(df_env, config$methods, method=stack_method)
          colname <- paste0("GPstacked_", stack_method, "_", env)
          df_rep[df_rep$Environment==env, colname] <- s1$yhat

          results[[length(results)+1]] <- data.frame(
            Trait=trait, Method="GPstacked", StackingMethod=stack_method,
            Rep=rep_i, Stage="Stage1_GP", Environment=env, Ability=s1$ability
          )
        }

        # Stage2: ME stacking
        me_ability <- stack_environments(df_rep, stack_method)
        results[[length(results)+1]] <- data.frame(
          Trait=trait, Method="GPstacked", StackingMethod=stack_method,
          Rep=rep_i, Stage="Stage2_ME", Ability=me_ability
        )
      }

      # Diagnostics
      diagnostics[[length(diagnostics)+1]] <- data.frame(
        Trait=trait, Rep=rep_i,
        N_env=n_distinct(df_rep$Environment),
        N_overlap=length(Reduce(intersect, split(df_rep$ID, df_rep$Environment)))
      )
    }
  }

  results <- bind_rows(results)

  # Stage3: MT stacking
  mt_results <- results %>%
    filter(Stage=="Stage2_ME", !is.na(Ability)) %>%
    group_by(Method, StackingMethod, Rep) %>%
    summarise(
      Ability = if(n()>=2) mean(Ability) else NA_real_,
      .groups="drop"
    ) %>%
    mutate(Trait="MT", Stage="Stage3_MT")

  list(
    results = bind_rows(results, mt_results),
    diagnostics = bind_rows(diagnostics)
  )
}

## ---------------------------
## SUMMARISE RESULTS
## ---------------------------
summarise_mtme <- function(fit) {
  fit$results %>%
    group_by(Trait, Method, StackingMethod, Stage) %>%
    summarise(
      MedianAbility = median(Ability, na.rm=TRUE),
      MeanAbility   = mean(Ability, na.rm=TRUE),
      SDAbility     = sd(Ability, na.rm=TRUE),
      .groups="drop"
    )
}

## ---------------------------
## RUN PIPELINE
## ---------------------------
mtme_results <- run_mtme_pipeline(config)
summary <- summarise_mtme(mtme_results)
write.csv(summary, "mtme_summary.csv", row.names=FALSE)
