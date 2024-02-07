library("targets")
library("tarchetypes")

# Define custom functions and other global objects.
sapply(list.files("R/", "20230908_functions", full.names = TRUE), source)

# Set generally needed packages.
tar_option_set(packages = c("tidyverse", 
                            "readxl",
                            "phyloseq", 
                            "Hmisc", 
                            "reshape2", 
                            "vegan", 
                            "broom", 
                            "DESeq2",
                            "ANCOMBC",
                            "CMAverse",
                            "LDM"))

# List of target objects
list(
  
  # Set up data -----
  tar_file(asv_file, "./data/Amplicons/mci.qc.geno.tax.tsv"),
  tar_target(ncer_phyloseq, phyloseq_setup(asv_file, meta_data)),
  tar_target(ncer_phy_trim, phyloseq_trim(ncer_phyloseq)), 
  tar_target(tax_levels, c("Phylum", "Class", "Order", "Family", "Genus")),
  tar_target(ncer_phyloseqs, phyloseq_summarize(ncer_phy_trim, tax_levels)),
  
  tar_files(meta_files,
            c("./data/Sample_data/2021-09-17_MCI_BIOME_clin.csv",
              "./data/Sample_data/2021-09-17_MCI_BIOME_clin2.csv",
              "./data/LuxPark_MCI_PD_data_codes.xlsx",
              "./data/LuxPark_MCI_PD_data.xlsx")),
  
  tar_target(meta_data, clin_setup(meta_files)),
  tar_target(clin_vars,
             list(num_vars = c("Age", 
                               "BMI", 
                               "Chao1", 
                               "Shannon", 
                               "InvSimpson"),
                  cat_vars = c("Gender",
                               "Years_of_Education", 
                               "First_Language",
                               "Living_With_Partner",
                               "BDI_I_mild", 
                               "APOE4",
                               "ATB_in_last_6_months"))),

  # Alpha diversity (species level) -----
  tar_target(ncer_adiv, adiv_setup(ncer_phy_trim)),
  
  # Beta diversity -----
  tar_target(ncer_bdiv, bdiv_ord(ncer_phy_trim, age_filter = FALSE)),
  tar_target(ncer_bdiv65, bdiv_ord(ncer_phy_trim, age_filter = TRUE)),
  tar_target(ad2_test_edu, bdiv_ad(ncer_bdiv, var = "Years_of_Education")),
  tar_target(ad2_test_edu65, bdiv_ad(ncer_bdiv65, var = "Years_of_Education")),
  
  # Differential abundance -----
  tar_target(basic_adjustments, c("Age", "Gender", "BMI", "ATB_in_last_6_months")),
  tar_target(further_adjustments, c("Living_With_Partner", "First_Language", "BDI_I_mild")),
  tar_target(APOE4, c("APOE4")),
  
  tar_target(dseq2_edu, ds2_main_var(ncer_phyloseqs, "Years_of_Education", adjust = c(basic_adjustments))),
  tar_target(dseq2_edu_further, ds2_main_var(ncer_phyloseqs, "Years_of_Education", adjust = c(basic_adjustments, further_adjustments))),
  tar_target(dseq2_edu_further_apoe, ds2_main_var(ncer_phyloseqs, "Years_of_Education", adjust = c(basic_adjustments, further_adjustments, APOE4))),
  
  tar_target(ancombc_edu, abc_main_var(phylo_list = ncer_phyloseqs, var = "Years_of_Education", adjust = c(basic_adjustments))),
  tar_target(ancombc_edu_further, abc_main_var(phylo_list = ncer_phyloseqs, var = "Years_of_Education", adjust = c(basic_adjustments, further_adjustments))),
  tar_target(ancombc_edu_further_apoe, abc_main_var(phylo_list = ncer_phyloseqs, var = "Years_of_Education", adjust = c(basic_adjustments, further_adjustments, "APOE4"))),
  
  # Mediation CMAverse -----
  tar_target(res11_Chao1, do_cmest(ncer_adiv, mediator = "zChao1", a = "11-16", mval = list(0))),
  tar_target(res16_Chao1, do_cmest(ncer_adiv, mediator = "zChao1", a = "16+", mval = list(0))),
  tar_target(res11_Chao1_multi, do_cmest(ncer_adiv, mediator = "zChao1", a = "11-16", mval = list(0), yreg = "multinomial")),
  tar_target(res16_Chao1_multi, do_cmest(ncer_adiv, mediator = "zChao1", a = "16+", mval = list(0), yreg = "multinomial")),
  tar_target(res11_Chao1_int, do_cmest(ncer_adiv, mediator = "zChao1", a = "11-16", mval = list(0), int = T)),
  tar_target(res16_Chao1_int, do_cmest(ncer_adiv, mediator = "zChao1", a = "16+", mval = list(0), int = T)),
  tar_target(res11_Chao1_int_multi, do_cmest(ncer_adiv, mediator = "zChao1", a = "11-16", mval = list(0), yreg = "multinomial", int = T)),
  tar_target(res16_Chao1_int_multi, do_cmest(ncer_adiv, mediator = "zChao1", a = "16+", mval = list(0), yreg = "multinomial", int = T)),
  
  tar_target(res11_Shannon, do_cmest(ncer_adiv, mediator = "zShannon", a = "11-16", mval = list(0))),
  tar_target(res16_Shannon, do_cmest(ncer_adiv, mediator = "zShannon", a = "16+", mval = list(0))),
  tar_target(res11_Shannon_int, do_cmest(ncer_adiv, mediator = "zShannon", a = "11-16", mval = list(0), int = T)),
  tar_target(res16_Shannon_int, do_cmest(ncer_adiv, mediator = "zShannon", a = "16+", mval = list(0), int = T)),
  
  tar_target(res11_InvSimpson, do_cmest(ncer_adiv, mediator = "zInvSimpson", a = "11-16", mval = list(0))),
  tar_target(res16_InvSimpson, do_cmest(ncer_adiv, mediator = "zInvSimpson", a = "16+", mval = list(0))),
  tar_target(res11_InvSimpson_int, do_cmest(ncer_adiv, mediator = "zInvSimpson", a = "11-16", mval = list(0), int = T)),
  tar_target(res16_InvSimpson_int, do_cmest(ncer_adiv, mediator = "zInvSimpson", a = "16+", mval = list(0), int = T))
)




