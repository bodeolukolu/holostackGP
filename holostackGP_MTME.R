#!/usr/bin/env Rscript

## ===========================
## CONFIGURATION
## ===========================
config <- list(
  traits = c("VIR1", "VIR2"),
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
    Trait = rep(c("VIR1", "VIR2"), each = 4),
    Environment = rep(c("ABI19A","ABI19B","NGE19B","NGE20A"), 2),
    stringsAsFactors = FALSE
  ),
  methods = c("MBLUP","rrBLUP","RKHS","BRR","BayesA","BayesB","BayesC","BL"),
  stacking_methods = c("linear","ridge","lasso","elasticnet","bayesian"),
  prediction_dir = "./",
  taxa_col = "Taxa",
  cv_rep_col = "rep"
)

## ===========================
## LIBRARIES
## ===========================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glmnet)
  library(BGLR)
})

MTME_type <- "specific"

if (MTME_type == "overall"){
  ## ===========================
  ## UTILITIES
  ## ===========================
  safe_cor <- function(y, yhat) {
    keep <- complete.cases(y, yhat)
    if (sum(keep) < 3) return(NA_real_)
    if (sd(y[keep]) == 0 || sd(yhat[keep]) == 0) return(NA_real_)
    cor(y[keep], yhat[keep])
  }

  ## ===========================
  ## LOAD DATA
  ## ===========================
  build_mtme_dataframe <- function(config) {
    message("Building unified MTME dataframe...")

    map_dfr(seq_len(nrow(config$environments)), function(i) {
      env <- config$environments[i, ]
      f <- file.path(config$prediction_dir, env$EnvFile)
      if (!file.exists(f)) stop("Missing file: ", f)

      df <- fread(f) %>% as.data.frame()

      if ("Row.names" %in% colnames(df))
        df$Row.names <- NULL

      trait <- env$Trait
      environment <- env$Environment

      # Identify phenotype column
      pheno_pattern <- paste0("^", trait, "_", environment, "$")
      pheno_col <- grep(pheno_pattern, colnames(df), value = TRUE)

      if (length(pheno_col) != 1) {
        stop(
          "Expected exactly one phenotype column matching ",
          trait, "_", environment,
          " in ", env$EnvFile,
          "\nFound: ", paste(pheno_col, collapse = ", ")
        )
      }

      # Check GP columns
      missing_gp <- setdiff(config$methods, colnames(df))
      if (length(missing_gp) > 0) {
        stop(
          "Missing GP columns in ", env$EnvFile, ": ",
          paste(missing_gp, collapse = ", ")
        )
      }

      df %>%
        transmute(
          ID          = .data[[config$taxa_col]],
          y           = as.numeric(.data[[pheno_col]]),
          Trait       = trait,
          Environment = environment,
          rep         = .data[[config$cv_rep_col]],
          across(all_of(config$methods), as.numeric)
        )
    })
  }

  ## ===========================
  ## STAGE 1: GP STACKING
  ## ===========================
  stack_gp_models <- function(df, gp_cols, method) {
    X <- as.matrix(df[, gp_cols, drop = FALSE])
    y <- df$y
    keep <- complete.cases(X, y)
    if (sum(keep) < 5 || ncol(X) < 2 || sd(y[keep]) == 0)
      return(list(yhat = rep(NA_real_, nrow(df)), ability = NA_real_))

    Xk <- scale(X[keep, , drop = FALSE])
    yk <- y[keep]

    fit_yhat <- tryCatch({
      if (method == "linear") {
        coef <- coef(lm(yk ~ Xk))
        as.numeric(cbind(1, Xk) %*% coef)
      } else if (method %in% c("ridge","lasso","elasticnet")) {
        alpha <- switch(method, ridge=0, lasso=1, elasticnet=0.5)
        fit <- cv.glmnet(Xk, yk, alpha = alpha)
        as.numeric(predict(fit, Xk, s = "lambda.min"))
      } else if (method == "bayesian") {
        fit <- BGLR(y = yk, ETA = list(list(X = Xk, model = "BRR")),
                    nIter = 6000, burnIn = 1000, verbose = FALSE)
        fit$yHat
      }
    }, error = function(e) NULL)

    yhat <- rep(NA_real_, nrow(df))
    if (!is.null(fit_yhat))
      yhat[keep] <- fit_yhat

    list(
      yhat = yhat,
      ability = safe_cor(y, yhat)
    )
  }

  ## ===========================
  ## STAGE 2: ME STACKING
  ## ===========================
  stack_me <- function(df) {

    # Aggregate across reps FIRST
    df_avg <- df %>%
      group_by(ID, Environment) %>%
      summarise(
        y = mean(y, na.rm = TRUE),
        GPstacked = mean(GPstacked, na.rm = TRUE),
        .groups = "drop"
      )

    # Pivot GP predictions
    X_df <- df_avg %>%
      select(ID, Environment, GPstacked) %>%
      pivot_wider(names_from = Environment, values_from = GPstacked)

    # Build y FROM THE SAME IDS AND ORDER
    y_df <- df_avg %>%
      group_by(ID) %>%
      summarise(y = mean(y, na.rm = TRUE), .groups = "drop")

    # Join explicitly to guarantee alignment
    model_df <- inner_join(y_df, X_df, by = "ID")

    y <- model_df$y
    X <- as.matrix(model_df[, setdiff(colnames(model_df), c("ID","y"))])

    keep <- complete.cases(X, y)
    if (sum(keep) < 5 || ncol(X) < 2 || sd(y[keep]) == 0)
      return(rep(NA_real_, nrow(model_df)))

    fit <- cv.glmnet(scale(X[keep, ]), y[keep], alpha = 0)

    yhat <- rep(NA_real_, nrow(model_df))
    yhat[keep] <- as.numeric(
      predict(fit, scale(X[keep, ]), s = "lambda.min")
    )

    yhat
  }

  ## ===========================
  ## STAGE 3: MT STACKING
  ## ===========================
  stack_mt <- function(df) {

    # Aggregate across reps
    df_avg <- df %>%
      group_by(ID, Trait, Environment) %>%
      summarise(
        y = mean(y, na.rm = TRUE),
        GPstacked = mean(GPstacked, na.rm = TRUE),
        .groups = "drop"
      )

    # Pivot predictions
    X_df <- df_avg %>%
      unite("TE", Trait, Environment, remove = FALSE) %>%
      select(ID, TE, GPstacked) %>%
      pivot_wider(names_from = TE, values_from = GPstacked)

    # Build y consistently
    y_df <- df_avg %>%
      group_by(ID) %>%
      summarise(y = mean(y, na.rm = TRUE), .groups = "drop")

    model_df <- inner_join(y_df, X_df, by = "ID")

    y <- model_df$y
    X <- as.matrix(model_df[, setdiff(colnames(model_df), c("ID","y"))])

    keep <- complete.cases(X, y)
    if (sum(keep) < 5 || ncol(X) < 2 || sd(y[keep]) == 0)
      return(rep(NA_real_, nrow(model_df)))

    fit <- cv.glmnet(scale(X[keep, ]), y[keep], alpha = 0)

    yhat <- rep(NA_real_, nrow(model_df))
    yhat[keep] <- as.numeric(
      predict(fit, scale(X[keep, ]), s = "lambda.min")
    )

    yhat
  }

  ## ===========================
  ## MAIN PIPELINE
  ## ===========================
  run_mtme_pipeline <- function(config) {

    data <- build_mtme_dataframe(config)

    results <- list()
    stage1_preds <- list()

    ## ---- Stage 1 ----
    for (trait in unique(data$Trait)) {
      for (env in unique(data$Environment[data$Trait == trait])) {
        df_env <- data %>% filter(Trait == trait, Environment == env)
        for (m in config$stacking_methods) {
          s1 <- stack_gp_models(df_env, config$methods, m)
          stage1_preds[[length(stage1_preds) + 1]] <- data.frame(
            ID = df_env$ID,
            Trait = trait,
            Environment = env,
            rep = df_env$rep,
            y = df_env$y,
            StackingMethod = m,
            GPstacked = s1$yhat
          )
          results[[length(results) + 1]] <- data.frame(
            Trait = trait,
            Environment = env,
            Method = "GPstacked",
            StackingMethod = m,
            Stage = "Stage1_GP",
            Ability = s1$ability
          )
        }
      }
    }

    stage1_preds <- bind_rows(stage1_preds)

    ## ---- Stage 2 ----
    for (trait in unique(stage1_preds$Trait)) {
      for (m in config$stacking_methods) {
        df_tm <- stage1_preds %>%
          filter(Trait == trait, StackingMethod == m)
        if (n_distinct(df_tm$Environment) < 2) next

        yhat2 <- stack_me(df_tm)
        envs <- unique(df_tm$Environment)

        for (env in envs) {
          idx <- df_tm$Environment == env
          results[[length(results) + 1]] <- data.frame(
            Trait = trait,
            Environment = env,
            Method = "GPstacked",
            StackingMethod = m,
            Stage = "Stage2_ME",
            Ability = safe_cor(df_tm$y[idx], yhat2[idx])
          )
        }
      }
    }

    ## ---- Stage 3 ----
    for (m in config$stacking_methods) {
      df_m <- stage1_preds %>% filter(StackingMethod == m)
      yhat3 <- stack_mt(df_m)

      for (trait in unique(df_m$Trait)) {
        for (env in unique(df_m$Environment)) {
          idx <- df_m$Trait == trait & df_m$Environment == env
          if (sum(idx) < 3) next
          results[[length(results) + 1]] <- data.frame(
            Trait = trait,
            Environment = env,
            Method = "GPstacked",
            StackingMethod = m,
            Stage = "Stage3_MT",
            Ability = safe_cor(df_m$y[idx], yhat3[idx])
          )
        }
      }
    }

    list(results = bind_rows(results))
  }

  ## ===========================
  ## SUMMARY
  ## ===========================
  summarise_mtme <- function(fit) {
    fit$results %>%
      group_by(Trait, Environment, Method, StackingMethod, Stage) %>%
      summarise(
        MedianAbility = median(Ability, na.rm = TRUE),
        MeanAbility   = mean(Ability, na.rm = TRUE),
        SDAbility     = sd(Ability, na.rm = TRUE),
        .groups = "drop"
      )
  }

  ## ===========================
  ## RUN
  ## ===========================
  mtme_results <- run_mtme_pipeline(config)
  summary <- summarise_mtme(mtme_results)
}

if (MTME_type == "specific") {
  ## ---------------------------
  ## UTILITIES
  ## ---------------------------
  safe_cor <- function(y, yhat) {
    keep <- complete.cases(y, yhat)
    if (sum(keep) < 5) return(NA_real_)
    if (sd(y[keep]) == 0 || sd(yhat[keep]) == 0) return(NA_real_)
    cor(y[keep], yhat[keep])
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

      pheno <- setdiff(
        names(df),
        c(config$methods, config$taxa_col, config$cv_rep_col)
      )

      if (length(pheno) != 1)
        stop("Expected exactly one phenotype in ", env$EnvFile)

      df %>%
        rename(y = all_of(pheno)) %>%
        mutate(
          Trait = env$Trait,
          Environment = env$Environment,
          ID = .data[[config$taxa_col]]
        )
    })
  }

  ## ---------------------------
  ## STAGE 1: GP STACKING (PER ENV × REP)
  ## ---------------------------
  stack_gp_models <- function(df_env, gp_cols, method) {

    gp_cols <- intersect(gp_cols, names(df_env))
    if (length(gp_cols) < 2)
      return(list(yhat = rep(NA_real_, nrow(df_env)), ability = NA_real_))

    X <- as.matrix(df_env %>% dplyr::select(all_of(gp_cols)))
    y <- df_env$y

    keep <- complete.cases(X, y)
    yhat <- rep(NA_real_, length(y))

    if (sum(keep) >= 5 && sd(y[keep]) > 0) {

      Xs <- scale(X[keep, , drop = FALSE])
      ys <- y[keep]

      preds <- tryCatch({
        if (method == "linear") {
          fit <- lm(ys ~ Xs)
          as.numeric(cbind(1, Xs) %*% coef(fit))
        } else if (method %in% c("ridge","lasso","elasticnet")) {
          alpha <- c(ridge = 0, lasso = 1, elasticnet = 0.5)[method]
          cv <- cv.glmnet(Xs, ys, alpha = alpha)
          as.numeric(predict(cv, Xs, s = "lambda.min"))
        } else if (method == "bayesian") {
          fit <- BGLR(
            y = ys,
            ETA = list(list(X = Xs, model = "BRR")),
            nIter = 6000,
            burnIn = 1000,
            verbose = FALSE
          )
          fit$yHat
        } else stop("Unknown stacking method")
      }, error = function(e) rep(NA_real_, length(ys)))

      yhat[keep] <- preds
    }

    list(
      yhat = yhat,
      ability = safe_cor(y, yhat)
    )
  }


  ## ---------------------------
  ## STAGE 2: MULTI-ENV STACKING (WITHIN REP)
  ## ---------------------------
  stack_environments <- function(df_rep, stack_method) {

    yhat_cols <- grep(
      paste0("^GPstacked_", stack_method, "_"),
      names(df_rep),
      value = TRUE
    )

    if (length(yhat_cols) < 2)
      return(NA_real_)

    df_use <- df_rep %>%
      dplyr::select(ID, y, all_of(yhat_cols)) %>%
      as.data.frame()

    if (nrow(df_use) < 5)
      return(NA_real_)

    X <- as.matrix(df_use[, yhat_cols, drop = FALSE])
    yhat_avg <- rowMeans(X, na.rm = TRUE)

    safe_cor(df_use$y, yhat_avg)
  }


  ## ---------------------------
  ## MAIN PIPELINE
  ## ---------------------------
  run_mtme_pipeline <- function(config) {

    data <- build_mtme_dataframe(config)
    results <- list()

    for (trait in unique(data$Trait)) {

      df_trait <- data %>% filter(Trait == trait)

      for (rep_i in sort(unique(df_trait[[config$cv_rep_col]]))) {

        df_rep <- df_trait %>%
          filter(.data[[config$cv_rep_col]] == rep_i)

        for (stack_method in config$stacking_methods) {

          # ----- Stage 1 -----
          for (env in unique(df_rep$Environment)) {

            df_env <- df_rep %>% filter(Environment == env)

            s1 <- stack_gp_models(df_env, config$methods, stack_method)

            colname <- paste0("GPstacked_", stack_method, "_", env)
            df_rep[df_rep$Environment == env, colname] <- s1$yhat

            results[[length(results) + 1]] <- data.frame(
              Trait = trait,
              Environment = env,
              Method = "GPstacked",
              StackingMethod = stack_method,
              Rep = rep_i,
              Stage = "Stage1_GP",
              Ability = s1$ability
            )
          }

          # ----- Stage 2 -----
          me_ability <- stack_environments(df_rep, stack_method)

          results[[length(results) + 1]] <- data.frame(
            Trait = trait,
            Environment = "ALL",
            Method = "GPstacked",
            StackingMethod = stack_method,
            Rep = rep_i,
            Stage = "Stage2_ME",
            Ability = me_ability
          )
        }
      }
    }

    results <- bind_rows(results)

    # ----- Stage 3: MULTI-TRAIT -----
    mt_results <- results %>%
      filter(Stage == "Stage2_ME", !is.na(Ability)) %>%
      group_by(StackingMethod, Rep) %>%
      summarise(
        Ability = mean(Ability),
        .groups = "drop"
      ) %>%
      mutate(
        Trait = "MT",
        Environment = "ALL",
        Method = "GPstacked",
        Stage = "Stage3_MT"
      )

    bind_rows(results, mt_results)
  }

  ## ---------------------------
  ## SUMMARY
  ## ---------------------------
  summarise_mtme <- function(results) {
    results %>%
      group_by(Trait, Environment, Method, StackingMethod, Stage) %>%
      summarise(
        MedianAbility = median(Ability, na.rm = TRUE),
        MeanAbility   = mean(Ability, na.rm = TRUE),
        SDAbility     = sd(Ability, na.rm = TRUE),
        .groups = "drop"
      )
  }

  ## ---------------------------
  ## RUN
  ## ---------------------------
  mtme_results <- run_mtme_pipeline(config)
  summary <- summarise_mtme(mtme_results)
}

write.csv(summary, "mtme_summary.csv", row.names = FALSE)
