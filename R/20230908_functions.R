# Data wrangling (basics)

mean_cl_and_bees <- function(df, varx, vary, test = NA, comparisons = NULL, cex = 2){
  p <- ggplot(df, aes_string(x = varx, y = vary, colour = varx)) +
    geom_quasirandom(alpha = .7, cex = cex) + 
    stat_summary(fun.data = "mean_cl_normal", geom = "point", size = 2, colour = "black") +
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = .2, size = 1, colour = "black") +
    ylab(gsub("_", " ", vary)) +
    theme_classic(base_size = 10) + # was 16
    theme(legend.position = "none", ) +
    ylab(vary) +
    xlab(str_replace_all(varx, "_", " ")) +
    scale_colour_manual(values = group_cols[[varx]])
  if(!is.na(test)){
    p <- p + stat_pwc(method = test, label = "p = {p.format}", ref.group = comparisons) #, bracket.size = .4)
  }
  return(p)}

clin_setup <- function(files){
  df1 <- read.csv(files[1], stringsAsFactors = TRUE)
  df2 <- read.csv(files[2], stringsAsFactors = TRUE)
  genetics_ids = read_excel(files[3])
  genetics_dat = read_excel(files[4])
  
  # Merge the dfs
  df <- merge(df1, df2)
  
  # Filter by age
  df = df %>% 
    filter(Group %in% c("HCNC", "HCMCI") & Age_at_assessment_years >= 50) %>% 
    left_join(genetics_ids %>% 
                left_join(genetics_dat) %>% 
                dplyr::select(-Gender, -Constipation, -Group, -Years_of_education))
  
  df_fin = df %>% 
    filter(Celiac_disease == "NO",
           Chronic_inflammatory_bowel_disease == "NO") %>% 
    drop_na(Group, 
            Age_at_assessment_years, 
            First_language, 
            Maritial_status, 
            Gender, 
            BDI_I, 
            BMI, 
            Years_of_education, 
            APOE_genotype, 
            ATB_in_last_6_months) %>% 
    mutate(MCI = factor(sub("HC", "", Group), levels = c("NC", "MCI")),
           MCI_bin = factor(ifelse(MCI == "MCI", 1, 0)),
           BDI_I_mild = factor(ifelse(BDI_I > 9, "Yes", "No"), levels = c("Yes", "No")),
           `Years_of_Education` = factor(case_when(
             Years_of_education <= 10 ~ "0-10",
             Years_of_education > 10 & Years_of_education <= 16 ~ "11-16",
             Years_of_education > 16 ~ "16+"),
             levels = c("0-10", "11-16", "16+")),
           `Living_With_Partner` = factor(case_when(
             Maritial_status %in% c("Never married", "Widowed", "Separated", "Divorced") ~ "No",
             Maritial_status %in% c("Married", "Domestic partnership") ~ "Yes"),
             levels = c("No", "Yes")),
           First_Language = factor(ifelse(grepl("French|Luxembourgish|German", First_language) == T, 
                                          "French/Luxembourgish/German", "Other"),
                                   levels = c("Other", "French/Luxembourgish/German")),
           Age = Age_at_assessment_years,
           Age_categorical = factor(case_when(
             Age < 65 ~ "<65",
             Age >= 65 & Age < 70 ~ "65-69",
             Age >= 70 ~ "70+")),
           Gender = factor(Gender, levels = c("Female", "Male")),
           APOE4_count = ifelse(grepl("E4/E4", APOE_genotype), 2,
                          ifelse(grepl("E4", APOE_genotype), 1, 0)),
           APOE4 = factor(ifelse(APOE4_count == 0, "None", "At Least 1"),
                          levels = c("At Least 1", "None")),
           ATB_in_last_6_months = factor(ATB_in_last_6_months, levels = c("YES", "NO"))) %>% 
    select(StoolKitID,
           MCI, MCI_bin, 
           Age, Age_categorical, Gender, First_Language, First_language, Living_With_Partner, 
           BDI_I_mild, BMI, 
           Years_of_Education, 
           APOE4, 
           ATB_in_last_6_months)
  
  # Set rownames
  rownames(df_fin) <- gsub("-", ".", df_fin$StoolKitID)
  
  return(df_fin)
}

clin_summary <- function(df, grpvar, clinvars){
  
  mean_sd_string <- function(x, d = 2){
    if(!is.na(mean(x, na.rm = TRUE)) & !is.na(sd(x, na.rm = TRUE))) {
      paste(round(mean(x, na.rm = TRUE), digits = d), "±",
            round(sd(x, na.rm = TRUE), digits = d))
    } else {
      return(NA)}
  }
  
  n_freqs <- function(x){
    if(sum(is.na(x)) == length(x)){
      return(NA)
    } else {
      paste(table(x, useNA = "ifany"), collapse = "/")}
  }
  
  # Continuous variables
  num_basics <- sapply(
    clinvars$num_vars,
    function(x) rbind(tapply(df[,x], df[, grpvar], mean_sd_string)))
  rownames(num_basics) <- levels(df[,grpvar])
  num_basics <- as.data.frame(t(num_basics))
  
  if(length(levels(df[,grpvar])) > 2){
    num_basics$p <- sapply(rownames(num_basics), function(x)
      summary(
        aov(as.formula(paste(x, "~", grpvar)),
            data = df))[[1]][["Pr(>F)"]][1])
    num_basics$p_test <- "anova"
  } else {
    num_basics$p <- 
      sapply(rownames(num_basics), function(x)
        tryCatch(
          t.test(as.formula(paste(x, "~", grpvar)),
                 df)$p.value, error=function(err) NA))
    num_basics[!is.na(num_basics$p), "p_test"] <- "t"
  }
  
  # Categorical variables
  cat_basics <- sapply(
    clinvars$cat_vars,
    function(x) rbind(tapply(df[,x], df[,grpvar], n_freqs)))
  rownames(cat_basics) <- levels(df[,grpvar])
  cat_basics <- as.data.frame(t(cat_basics))
  cat_basics$p <- sapply(rownames(cat_basics), function(x)
    fisher.test(table(df[,c(x, grpvar)]))$p.value)
  cat_basics$p_test <- "fisher"
  rownames(cat_basics) <- paste0(rownames(cat_basics), " (",
                                 sapply(df[,rownames(cat_basics)],
                                        function(x) paste(sort(unique(x),
                                                               na.last = TRUE),
                                                          collapse = "/")), ")")
  
  # Put together
  df_out <- rbind(n = c(table(df[,grpvar]), NA, NA),
                  cat_basics,
                  num_basics)
  
  colnames(df_out) <- gsub("_", " ", colnames(df_out))
  rownames(df_out) <- gsub("_", " ", rownames(df_out))
  
  return(df_out)
}

# Data wrangling (phyloseq)

phyloseq_setup <- function(file, df){
  full_df <- read.table(file, sep = "\t", header = TRUE, row.names = 1)
  
  # Set up count table
  counts_mat <- as.matrix(full_df[, grep("mothur",
                                          colnames(full_df),
                                          invert = TRUE)])
  colnames(counts_mat) <- sub("^NCER", "",
                              sub("\\.STLOMN.*", "", colnames(counts_mat)))
    
  # Set up taxonomy
  tax_mat <- full_df[,grep("mothur", colnames(full_df))]
  tax_mat <- tax_mat[,2:ncol(tax_mat)]
  colnames(tax_mat) <- sub("\\.mothur", "", colnames(tax_mat))
  tax_mat <- as.matrix(tax_mat)
  
  # Fill empty spaces with "unclassified"
  tax_mat[tax_mat == ""] <- "unclassified"
  
  tax_mat[grep("metagenome", tax_mat)] <- "unclassified"
  
  tax_mat[grep("unidentified", tax_mat)] <- "unclassified"
  
  # Fix some problematic taxon names
  tax_mat <- gsub("-", "_", tax_mat)
  
  # Additional cleanup
  
  # Make sure sample names match
  # (there is one duplicate that needs to be kicked out of the ASV data)
  counts_mat <- counts_mat[,intersect(colnames(counts_mat),
                                     rownames(df))]
  df <- df[intersect(colnames(counts_mat),
                     rownames(df)),]
  
  # Make phyloseq
  final_phy <- phyloseq(otu_table(counts_mat, taxa_are_rows = TRUE),
                       tax_table(tax_mat),
                       sample_data(df))
}

phyloseq_trim <- function(phylo_obj){
  # Trim zero or singleton taxa
  phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 1, phylo_obj)
  
  # Trim out taxa present in less than 1/10th of samples
  phylo_obj <- filter_taxa(phylo_obj, function(x)
    sum(x > 0) > round(nsamples(phylo_obj)/10), prune = TRUE)

  # Trim any remaining contaminants/misclassifications that
  # are labeled as Mitochondria or Chloroplasts
  phylo_obj <- subset_taxa(phylo_obj, Order != "Chloroplast" &
                                 Family != "Mitochondria")
  
  # Trim samples with < 10000 reads
  phylo_obj <- prune_samples(sample_sums(phylo_obj) > 10000, phylo_obj)
}

phyloseq_summarize <- function(phylo_obj, levels){

  # Function for collapsing to a specific taxonomic level
  collapseTaxLevel <- function(phylo_obj, level){
    
    levelNum <- grep(level, colnames(tax_table(phylo_obj)))
    
    otu <- cbind(tax_table(phylo_obj)[, levelNum],
                 as.data.frame(as(otu_table(phylo_obj), "matrix")))
    otu <- melt(otu, id = level)
    otu <- acast(otu, as.formula(paste(level, "~variable", sep = "")),  sum)
    
    tax <- unique(as(tax_table(phylo_obj), "matrix")[,1:levelNum])
    rownames(tax) <- unname(tax[,levelNum])
    
    # Combine unclassifieds at the chosen level into
    # one "unclassified" bin (if there are any)
    if(length(grep("unclassified", tax[,level]) > 0)){
      tax <- tax[-grep("unclassified", tax[,level]),]
      tax <- rbind(tax, unclassified = rep("unclassified", ncol(tax)))
    }
    
    tax <- tax[order(rownames(tax)),]
    
    new_phylo_obj <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                              sample_data(phylo_obj),
                              tax_table(tax))
    return(new_phylo_obj)
  }
  
  # Phyloseq objects for each level:
  phyloseqs_list <- c(
    lapply(levels, function(x)
      collapseTaxLevel(phylo_obj, x)),
    phylo_obj)
  names(phyloseqs_list) <- c(levels, "ASV")
  
  return(phyloseqs_list)
}

# Alpha diversity

adiv_setup <- function(phylo_obj){
  # Subsample to even depth
  phylo_obj_r <- rarefy_even_depth(phylo_obj, rngseed = 8299347)
  
  # Calculate alpha diversity index values and
  # combine with metadata
  adivs <- data.frame(
    estimate_richness(phylo_obj_r,
                      measures = c("Chao1", "Shannon", "InvSimpson")),
    as(sample_data(phylo_obj_r), "data.frame"))
  
  adivs = adivs %>% 
    mutate(zChao1 = as.vector(scale(Chao1)),
           zShannon = as.vector(scale(Shannon)),
           zInvSimpson = as.vector(scale(InvSimpson)))
  
  return(adivs)
}

# Beta dispersion

bdiv_ord <- function(phylo_obj, age_filter){
  # Subsample to even depth
  if(age_filter) {phylo_obj = subset_samples(phylo_obj, Age >= 65)}
  
  phylo_obj_r <- rarefy_even_depth(phylo_obj, rngseed = 3541451)
  
  # Calculate ordination
  ord <- ordinate(phylo_obj_r, method = "NMDS", distance = "bray",
                       trymax = 500, trace = FALSE)
  
  # Calculate BC distances
  set.seed(123)
  cat(" seed = 123 initialized")
  dist <- vegdist(
    as.data.frame(t(as(otu_table(phylo_obj_r), "matrix"))),
    method = "bray")
  
  ord_obj <- list(phyloseq = phylo_obj_r, ord = ord, dist = dist)
}

bdiv_ad <- function(bdiv_list, var){
  phylo_obj <- bdiv_list[["phyloseq"]]
  dist <- bdiv_list[["dist"]]
  
  res <- adonis2(as.formula(paste0("dist ~ ", 
                                   var, 
                                   "+ Gender + Age + ATB_in_last_6_months + BDI_I_mild + First_Language + Living_With_Partner + APOE4")),
          data = as(sample_data(phylo_obj), "data.frame"),
          perm = 999, by = "margin")
}

# Differential abundance

getTaxRA <- function(taxon, level, phylo_obj){
  log10(
    t(
      prop.table(
        otu_table(phylo_obj),
        2)
      )[,taxon] + 0.001
    )
}

ds2_main_var <- function(phylo_list, var, adjust){
  
  prevalenceTrim <- function(phylo_obj, acutoff = 0, ncutoff = 4){
    phylo_obj <- filter_taxa(phylo_obj,
                             function(x) sum(x > acutoff) > 
                               nsamples(phylo_obj)/ncutoff,
                             prune = TRUE)
    return(phylo_obj)
  }
  
  ds2conf <- function(phylo_obj){
    
    phylo_obj <- prevalenceTrim(phylo_obj)
    phylo_obj <- prune_samples(
      !is.na(unlist(sample_data(phylo_obj)[,var])), phylo_obj)
    
    sample_data(phylo_obj)$Years_of_Education = factor(case_when(
      sample_data(phylo_obj)$Years_of_Education == "0-10" ~ "low",
      sample_data(phylo_obj)$Years_of_Education == "11-16" ~ "int",
      sample_data(phylo_obj)$Years_of_Education == "16+" ~ "high"), levels = c("low", "int", "high"))
    sample_data(phylo_obj)$First_Language = factor(case_when(
      sample_data(phylo_obj)$First_Language == "French/Luxembourgish/German" ~ "native",
      sample_data(phylo_obj)$First_Language == "Other" ~ "other"), levels = c("other", "native"))
    sample_data(phylo_obj)$APOE4 = factor(case_when(
      sample_data(phylo_obj)$APOE4 == "None" ~ "None",
      sample_data(phylo_obj)$APOE4 == "At Least 1" ~ "At.least.1"), levels = c("At.least.1", "None"))
    
    ds2 <- phyloseq_to_deseq2(
              phylo_obj, 
              as.formula(
                paste0("~ ", 
                       paste(adjust, collapse = " + "), " + ",
                       var))) # needs to be last

    ds2 <- DESeq(ds2, fitType = "local", sfType = "poscounts")
    
    full_res <- do.call(
      "rbind", 
      lapply(
        data.frame(combn(levels(as_factor(unlist(sample_data(phylo_obj)[,var]))), 2)),
        function(x) data.frame(results(ds2,
                                       contrast = c(var, x[1], x[2]),
                                       alpha = 0.05),
                               Taxon = rownames(ds2),
                               Contrast = paste(x[1], x[2], sep = "."))))
    rownames(full_res) <- NULL
    
    return(full_res)
  }
  
  ds2_res <- do.call("rbind",
                     lapply(names(phylo_list), function(x)
                       data.frame(ds2conf(phylo_list[[x]]),
                                  Level = x)))
}

ds2_restable = function(res, wide = FALSE){
  
  res = subset(res, padj < 0.05)[,c("Contrast", "Taxon", "Level", "baseMean", "log2FoldChange", "lfcSE", "padj")] %>% 
    arrange(Contrast, Level)
  
  if(nrow(res) > 0){
    if(wide){
      res %>% 
        pivot_wider(id_cols = c(Level, Taxon),
                    names_from = Contrast,
                    values_from = padj) %>% 
        kable(booktabs = TRUE, 
              digits = c(NA, NA, 3, 3),
              row.names = FALSE) %>% 
        add_footnote("Shown results are restricted to p < .05.", notation = "alphabet") %>% 
        kable_classic(full_width = FALSE, html_font = "arial") %>% 
        kable_styling(font_size = 16, position = "l") %>% 
        row_spec(0, bold = TRUE)
    }
    
    else
      res %>% kable(booktabs = TRUE, 
                    digits = c(NA, NA, NA, 1, 2, 2, 3),
                    row.names = FALSE) %>% 
        add_footnote("Shown results are restricted to p < .05.", notation = "alphabet") %>% 
        kable_classic(full_width = FALSE, html_font = "arial") %>% 
        kable_styling(font_size = 16, position = "l") %>% 
        row_spec(0, bold = TRUE)
  }
  else {cat("No significant hits.")}
}

abc_main_var <- function(phylo_list, var, adjust){
  
  abcGroup <- function(phylo_obj){
    
    sample_data(phylo_obj)$Years_of_Education = factor(case_when(
      sample_data(phylo_obj)$Years_of_Education == "0-10" ~ "low",
      sample_data(phylo_obj)$Years_of_Education == "11-16" ~ "int",
      sample_data(phylo_obj)$Years_of_Education == "16+" ~ "high"), levels = c("low", "int", "high"))
    
    # differential abundance comparison
    out <- ancombc(phyloseq = phylo_obj,
                   formula = paste0(paste(adjust, collapse = " + "), " + ", var),
                   group = var,
                   p_adj_method = "fdr",
                   prv_cut = .25,
                   struc_zero = TRUE, neg_lb = TRUE,
                   global = T, tax_level = NULL)
    
    # may be included as function argument to alter applicability to other variables
    group_contr <- c("Years_of_Educationint", "Years_of_Educationhigh")
    
    res <- do.call("rbind",
                   lapply(group_contr,
                          function(x) sapply(out$res, function(y) y[, x])))
    res <- as.data.frame(res)
    res$Taxon <- rep(rownames(out$feature_table), length(group_contr))
    res$Contrast <- rep(group_contr,
                        each = nrow(out$feature_table))
    
    
    if(!is.null(out$res_global)){
      res <- rbind(res,
                   data.frame(lfc = NA,
                              se = NA,
                              W = out$res_global$W,
                              p_val = out$res_global$p_val,
                              q_val = out$res_global$q_val,
                              diff_abn = out$res_global$diff_abn,
                              Taxon = out$res_global$taxon,
                              Contrast = "global",
                              row.names = NULL))
    }
    
    return(res)
  }
  
  abc_res <- do.call("rbind",
                     lapply(names(phylo_list), function(x)
                       data.frame(abcGroup(phylo_list[[x]]),
                                  Level = x)))
}

abc_restable = function(res, wide = FALSE){
  
  res = subset(res, q_val < 0.05)[,c("Contrast", "Level", "Taxon", "lfc", "q_val")]
  
  if(nrow(res) > 0){
    if(wide){
      res %>%
        pivot_wider(id_cols = c(Level, Taxon),
                    names_from = Contrast,
                    values_from = q_val) %>% 
        kable(booktabs = TRUE, 
              digits = c(NA, NA, 3, 3), 
              row.names = FALSE) %>% 
        add_footnote("Shown results are restricted to p < .05.", notation = "alphabet") %>% 
        kable_classic(full_width = FALSE, html_font = "arial") %>% 
        kable_styling(font_size = 16, position = "l") %>% 
        row_spec(0, bold = TRUE)
    }
    
    else
      res %>% 
        kable(booktabs = TRUE, 
              digits = c(NA, NA, NA, 2, 3), 
              row.names = FALSE) %>% 
        add_footnote("Shown results are restricted to p < .05.", notation = "alphabet") %>% 
        kable_classic(full_width = FALSE, html_font = "arial") %>% 
        kable_styling(font_size = 16, position = "l") %>% 
        row_spec(0, bold = TRUE)
  }
  
  else cat("No significant hits.")
}

plot_ra_hits = function(i){
  p = ggplot(tax_ra[!duplicated(names(tax_ra))], aes_string(x = "Years_of_Education", y = top_hits[i], colour = "Years_of_Education")) +
    geom_quasirandom(alpha = .7, cex = 1.5) + # was 3
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      geom = "linerange", size = .7, colour = "black") + # was 1.3
    stat_summary(
      fun = median,
      geom = "point", size = 5, colour = "black") +
    theme_classic(base_size = 10) + # was 16
    theme(legend.position = "none", 
          axis.ticks = element_blank(), 
          panel.grid = element_blank(), 
          strip.text = element_text(face = "bold")) +
    ylab("Log10 (RA + 0.001)") + xlab("Years of Education") +
    scale_colour_manual(values = group_cols[["Years_of_Education"]]) +
    ggtitle(str_replace_all(str_replace_all(top_hits[i], "OTU", "ASV"), "_", " ")) +
    scale_fill_viridis_d(option = "D")
  return(p) 
}

# Mediation

do_cmest = function(data, mediator, a, astar = "0-10", mval, yreg = "logistic", int = FALSE, print = F){
  set.seed(123)
  med_analysis_1 <- cmest(data = data %>% mutate(MCI_bin = as.numeric(as.character(MCI_bin))), 
                          model = "rb", 
                          outcome = "MCI_bin", 
                          exposure = "Years_of_Education", 
                          mediator = c(mediator), 
                          basec = c("Age", 
                                    "Gender",
                                    "ATB_in_last_6_months", 
                                    #"BMI",
                                    "BDI_I_mild", 
                                    "First_Language",
                                    "Living_With_Partner", 
                                    "APOE4"),
                          mreg = list("linear"), 
                          yreg = yreg, 
                          EMint = int,
                          astar = astar, 
                          a = a,
                          mval = mval,
                          yval = 1,  
                          estimation = "imputation",  
                          inference = "bootstrap", 
                          nboot = 5000) 
  if(print) {med_analysis_1 %>% summary() %>% print()}
  return(med_analysis_1)
}

neat_cmest_res = function(cmest_res, cap){
  res = data.frame(1:length(cmest_res$effect.pe))
  res %>%
    transmute("Effect" = cmest_res$effect.ci.high %>% names,
           "Estimate [95% CI]" = paste0(cmest_res$effect.pe %>% round(digits = 2) %>% format(nsmall = 2), " [", 
                                        cmest_res$effect.ci.low %>% round(digits = 2) %>% format(nsmall = 2), " to ", 
                                        cmest_res$effect.ci.high %>% round(digits = 2) %>% format(nsmall = 2), "]"),
           "P Value" = paste(cmest_res$effect.pval %>% Hmisc::format.pval(digits = 3, eps = 0.001, nsmall = 3))) %>% 
    mutate(Sig = case_when(`P Value` >= .05 ~ "",
                            `P Value` < .05 & `P Value` >= .01 ~ "'*'",
                            `P Value` < .01 & `P Value` >= .001 ~ "'**'", 
                            `P Value` < .001 ~ "'***'")) %>% 
    kable(caption = cap, "html") %>% 
    kable_classic(full_width = FALSE, html_font = "arial") %>% 
    kable_styling(font_size = 16, position = "l") %>% 
    row_spec(0, bold = TRUE)
}

cmest_coef_ci = function(cmest_res){
  
  yres = data.frame(
    cmest_res$reg.output$yreg$coefficients,
    summary(cmest_res$reg.output$yreg)$coefficients[,4],
    confint(cmest_res$reg.output$yreg))
  
  names(yres) = c("coef", "p.value", "ci.low", "ci.high")
  yres = yres %>% 
    transmute("Estimate [95% CI]" = paste0(coef %>% round(digits = 2) %>% format(nsmall = 2), " [", 
                                           ci.low %>% round(digits = 2) %>% format(nsmall = 2), " to ", 
                                           ci.high %>% round(digits = 2) %>% format(nsmall = 2), "]"),
              "P Value" = paste(p.value %>% Hmisc::format.pval(digits = 3, eps = 0.001, nsmall = 3))) %>% 
    mutate(Sig = case_when(`P Value` >= .05 ~ "",
                            `P Value` < .05 & `P Value` >= .01 ~ "'*'",
                            `P Value` < .01 & `P Value` >= .001 ~ "'**'", 
                            `P Value` < .001 ~ "'***'"))
  
  mres = data.frame(
    cmest_res$reg.output$mreg[[1]]$coefficients,
    summary(cmest_res$reg.output$mreg[[1]])$coefficients[,4],
    confint(cmest_res$reg.output$mreg[[1]]))
  
  names(mres) = c("coef", "p.value", "ci.low", "ci.high")
  mres = mres %>% 
    transmute("Estimate [95% CI]" = paste0(coef %>% round(digits = 2) %>% format(nsmall = 2), " [", 
                                           ci.low %>% round(digits = 2) %>% format(nsmall = 2), " to ", 
                                           ci.high %>% round(digits = 2) %>% format(nsmall = 2), "]"),
              "P Value" = paste(p.value %>% Hmisc::format.pval(digits = 3, eps = 0.001, nsmall = 3))) %>% 
    mutate(Sig = case_when(`P Value` >= .05 ~ "",
                            `P Value` < .05 & `P Value` >= .01 ~ "'*'",
                            `P Value` < .01 & `P Value` >= .001 ~ "'**'", 
                            `P Value` < .001 ~ "'***'"))
  
  kables(list(
    yres %>% 
      kable(caption = "Outcome Model", "html") %>% 
      kable_classic(full_width = FALSE, html_font = "arial") %>% 
      kable_styling(font_size = 16, position = "left") %>% 
      row_spec(0, bold = TRUE),
    mres %>% 
      kable(caption = "Mediator Model", "html") %>% 
      kable_classic(full_width = FALSE, html_font = "arial") %>% 
      kable_styling(font_size = 16, position = "left") %>% 
      row_spec(0, bold = TRUE))) 
}


