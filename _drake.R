options(future.globals.maxSize = +Inf)
source("R/packages.R")

project_folder = "source_data/"
snp_data_file = "GSA_chr7.vcf.gz"
#snp_data_file = "first_1000_chr7.vcf"
flow_data_file = "Flu_Project_Subset_for_IKZF1_2_quant_trait_analysis_wdedit.xlsx"
lower_bound = 50204068
upper_bound = 50505101

BPPARAM = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)
future::plan(future::multisession)

source("R/plan.R")
drake_config(plan = analysis_plan,
             verbose = 2,
             parallelism = "future",
             jobs = parallel::detectCores(),
             lock_envir = TRUE)
