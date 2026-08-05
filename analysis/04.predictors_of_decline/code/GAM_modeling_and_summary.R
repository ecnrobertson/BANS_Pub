fit_esugams_and_summarize <- function(
    data_path,
    out_dir,
    run_name,
    response = "abd_trend",
    esu_col = "ESU",
    focal_vars,
    var_map,
    k_global,
    k_fs,
    drop_esu_levels = NULL,
    na_action = "drop",
    save_models = TRUE,
    overwrite = FALSE
) {
  stopifnot(file.exists(data_path))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  run_dir <- file.path(out_dir, run_name)
  if (dir.exists(run_dir) && !overwrite) {
    stop("run_dir already exists and overwrite=FALSE: ", run_dir)
  }
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  
  suppressPackageStartupMessages({
    library(mgcv)
  })
  
  df <- read.csv(data_path)
  
  # Basic checks
  needed_cols <- unique(c(response, esu_col, unname(var_map[focal_vars])))
  missing_cols <- setdiff(needed_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in data: ", paste(missing_cols, collapse=", "))
  }
  
  # Handle NA strategy
  if (na_action == "drop") {
    df <- df[complete.cases(df[, needed_cols]), , drop = FALSE]
  } else if (na_action == "fail") {
    if (any(!complete.cases(df[, needed_cols]))) {
      stop("NA values present in required columns; set na_action='drop' to remove.")
    }
  } else {
    stop("na_action must be 'drop' or 'fail'")
  }
  
  # Drop ESU levels if requested
  df[[esu_col]] <- as.factor(df[[esu_col]])
  if (!is.null(drop_esu_levels)) {
    df <- df[!(df[[esu_col]] %in% drop_esu_levels), , drop = FALSE]
    df[[esu_col]] <- droplevels(df[[esu_col]])
  }
  
  # Helper: extract R2
  get_r2 <- function(m) {
    s <- summary(m)
    c(R2_adj = s$r.sq, Dev_expl = s$dev.expl)
  }
  
  # Helper: anova row
  get_anova_row <- function(m0, m1, label) {
    a <- anova(m0, m1, test = "Chisq")
    data.frame(
      Comparison = label,
      Delta_EDF  = a$Df[2],
      Delta_Dev  = a$Deviance[2],
      p_value    = a$`Pr(>Chi)`[2],
      stringsAsFactors = FALSE
    )
  }
  
  # Build base RHS (global smooths for all focal vars)
  # Example: s(ppt_june_delta,k=15) + s(soil_june_delta,k=20) + ...
  global_terms <- vapply(
    focal_vars,
    function(v) sprintf("s(%s, k=%d, bs='tp')", var_map[[v]], k_global[[v]]),
    character(1)
  )
  rhs_base <- paste(global_terms, collapse = " + ")
  
  # Fit base models (same for all predictors)
  f_null <- as.formula(sprintf("%s ~ %s", response, rhs_base))
  gam_null2 <- gam(f_null, data = df, method = "REML")
  
  f_modG <- as.formula(sprintf("%s ~ %s + s(%s, k=2, bs='re')", response, rhs_base, esu_col))
  modG <- gam(f_modG, data = df, method = "REML")
  
  if (save_models) {
    saveRDS(gam_null2, file.path(run_dir, "gam_null2.rds"))
    saveRDS(modG,      file.path(run_dir, "modG.rds"))
  }
  
  # Loop over predictors, fit GS and S variants, and write per-predictor CSVs
  all_rows <- list()
  model_paths <- list()
  
  for (v in focal_vars) {
    x <- var_map[[v]]
    
    # modGS: global smooth + factor-smooth deviation
    f_gs <- as.formula(sprintf(
      "%s ~ %s + s(%s, %s, k=%d, bs='fs')",
      response, rhs_base, x, esu_col, k_fs
    ))
    modGS <- gam(f_gs, data = df, method = "REML")
    
    # modS: factor-smooth only for focal var (no global smooth for x)
    rhs_no_x <- paste(global_terms[focal_vars != v], collapse = " + ")
    f_s <- as.formula(sprintf(
      "%s ~ %s + s(%s, %s, k=%d, bs='fs')",
      response, rhs_no_x, x, esu_col, k_fs
    ))
    modS <- gam(f_s, data = df, method = "REML")
    
    if (save_models) {
      p1 <- file.path(run_dir, paste0("modGS_", v, ".rds"))
      p2 <- file.path(run_dir, paste0("modS_",  v, ".rds"))
      saveRDS(modGS, p1)
      saveRDS(modS,  p2)
      model_paths[[paste0("modGS_", v)]] <- p1
      model_paths[[paste0("modS_",  v)]] <- p2
    }
    
    # R2 values
    r2_null <- get_r2(gam_null2)
    r2_G    <- get_r2(modG)
    r2_GS   <- get_r2(modGS)
    r2_S    <- get_r2(modS)
    
    # comparisons: Null -> modG -> modGS -> modS
    res <- rbind(
      data.frame(
        Predictor = v,
        Model = "Null",
        Comparison = "Null (no ESU)",
        Delta_EDF = NA, Delta_Dev = NA, p_value = NA,
        R2_adj = r2_null["R2_adj"], Dev_expl = r2_null["Dev_expl"],
        stringsAsFactors = FALSE
      ),
      cbind(
        data.frame(Predictor=v, Model="modG", stringsAsFactors=FALSE),
        get_anova_row(gam_null2, modG, "modG vs Null"),
        t(r2_G)
      ),
      cbind(
        data.frame(Predictor=v, Model="modGS", stringsAsFactors=FALSE),
        get_anova_row(modG, modGS, paste0("modGS_", v, " vs modG")),
        t(r2_GS)
      ),
      cbind(
        data.frame(Predictor=v, Model="modS", stringsAsFactors=FALSE),
        get_anova_row(modGS, modS, paste0("modS_", v, " vs modGS")),
        t(r2_S)
      )
    )
    
    # p-value formatting (optional)
    res$p_value <- ifelse(
      is.na(res$p_value), NA,
      ifelse(res$p_value < 0.001, "<0.001", as.character(round(res$p_value, 3)))
    )
    
    # write per-predictor csv
    out_csv <- file.path(run_dir, paste0("GAM_model_comparison_", v, ".csv"))
    write.csv(res, out_csv, row.names = FALSE)
    
    all_rows[[v]] <- res
  }
  
  # write master csv
  master <- do.call(rbind, all_rows)
  write.csv(master, file.path(run_dir, "GAM_model_comparison_ALL.csv"), row.names = FALSE)
  
  invisible(list(
    run_dir = run_dir,
    n = nrow(df),
    esu_levels = levels(df[[esu_col]]),
    model_paths = model_paths
  ))
}
