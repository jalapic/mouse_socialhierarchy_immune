###############################################################
# Data and R code for Lee et al.
###############################################################


source('Won_shortcut.R')

theme_set(theme_bw(base_size = 10))

# figure save template
# unhash and change figure names  
# png(filename = glue("results_figures/xx_{my_tissue}.png"),
#     width = 18, height = 12, units = "cm", res = 600)
# invisible(dev.off())


###############################################################
# Install packages that you don't have on your computer 
###############################################################
# source('code_carpentry/00_install_packages.R') # unhash when you are running this first time

###############################################################
# Load libraries  
###############################################################
source('code_carpentry/00_load_libraries.R')

###############################################################
# Preset global ggplot styles 
###############################################################
backup_options <- options() # in case I want to come back to default options 
# options(backup_options) # unhash and use it when you want default graph settings
source('code_carpentry/00_preset_ggplot_options.R')


###############################################################
# Load custom functions 
###############################################################000
source("code_functions/functions_behavior.R")
source("code_functions/functions_plot.R")
source("code_functions/geom_boxjitter_ggparl_erocoar_github_notmine.R")
source("code_functions/functions_network.R")
source("code_functions/functions_rnaseq.R")


###############################################################
# 
# Analysis of Social Behavior
#
###############################################################

###############################################################
# 1. Clean up raw behavioral observation data
###############################################################
source("code_carpentry/01_import_data.R") # raw behavioral data and biological data 
source("code_carpentry/02_fix_errors.R") # fix observation error as noted on the comment column
source("code_carpentry/03_list_dataframes.R") # put behavior datafrmaes into a list 
source("code_carpentry/04_behavior_dataframe.R") # tidy up dfs for downstream analysis

###############################################################
# 2. Behavior: Social hierarchy characteristics  
###############################################################
source("code_analysis/Behav01_observation_descriptives.R")
source("code_analysis/Behav02_hierarchy_characteristics.R")
source("code_analysis/Behav03_ranking_ds.R") # social rank by David's score
source("code_analysis/Behav04_isi_matrixplot.R") # I&SI ranking 
source("code_analysis/Behav05_glicko.R") # Glicko ranking 
source("code_analysis/Behav06_network_measures.R") # aggressive behavior social network 

###############################################################
# 3. Behavior: Behavioral phenotype in response to aggressive encounters
###############################################################
source("code_carpentry/05_get_subordinate_behavior.R") # clean up df more to extract subordinate behavior pattern
source("code_analysis/Behav07_response_cluster_analysis.R") # cluster individuals by behavior observation data

source("code_carpentry/06_gather_all_behavior_measures.R") # merge 

source("code_analysis/Behav08_subordinate_behaviors.R") # multinomial regression: response phenotype ~ social status
source("code_analysis/Behav09_kmeans_cluster_behavior_explore.R") # which behavior drives the cluster assignment
source("code_analysis/Behav10_kmeans_cluster_behavior_stat.R") 


######## SHORTCUT to all the analysis above: ##################
# resultsdf <- readRDS("results_statRDS/resultsdf.RDS")
# all_behavior <- readRDS("results_statRDS/all_behavior.RDS")
###############################################################


###############################################################
# 
# Analysis of Biological measurments
# a) Immunophenotyping - Flow cytometry
# b) fkbp5 DNA methylation level 
# c) corticosterone assay (before and after social housing)
# d) RNA sequencing data (Liver and Spleen)
#
###############################################################

###############################################################
# 4. Clean up raw biological measurememt data, then gather all in one dataframe
###############################################################
source("code_carpentry/07_get_flow_data.R")
source("code_carpentry/08_get_cort_data.R")
source("code_carpentry/09_gather_alldata.R")

######## SHORTCUT to all behavioral and biological data above #
# flow_all <- readRDS('results_statRDS/flow_all.RDS')
# alldata <- readRDS('results_statRDS/alldata.RDS') # except RNAseq data
###############################################################


###############################################################
# 5. Anogenital distance (AGD) and spleen weight (GD14)
###############################################################
source("code_analysis/Bio01_spleen_weight.R") # spleen weight ~ social status 
source("code_analysis/Bio02_anogenital_distance.R") #does AGD predict pair or group status?

###############################################################
# 6.Fkbp5 DNA methylation data 
###############################################################
source("code_analysis/Bio03_fkbp5_methylation_data.R")

source("code_analysis/Bio_Figure_spleen_fkbp5.R") # figure for manuscript

###############################################################
# 7. Flow cytometry data - Immunophenotyping
###############################################################
source("code_analysis/Bio04_flow_data_PD09_pairstatus.R")
source("code_analysis/Bio05_flow_data_PD09_predict_finalrank.R")

source("code_analysis/Bio06_flow_data_GD14_davids_score.R")
source("code_analysis/Bio07_flow_data_GD14_dominance_certainty.R")
source("code_analysis/Bio08_flow_data_GD14_kmeans_cluster.R")
source("code_analysis/Bio09_flow_data_plasticity_PD09_to_GD14_davids_score.R")

source("code_analysis/Bio10_flow_data_GD14_social_network_measures.R")
source("code_analysis/Bio_Figure_PD09_GD14_brms_estimates.R") # figure for manuscript
source("code_analysis/Bio_Figure_flow_data_plasticity.R") # figure for manuscript

###############################################################
# 8. Corticosterone assay
###############################################################
source("code_analysis/Bio11_corticosterone_status.R")
source("code_analysis/Bio12_corticosterone_against_other_bio_measures.R")

###############################################################
# 9. Looking at the big picture: across all biological data (PCA)
###############################################################
source("code_analysis/Bio13_extract_principal_components_from_all_bio_measures.R")
source("code_analysis/Bio14_all_biolmeasure_alpha_despotism.R")

###############################################################
# 10. RNAseq data 
# Sample quality check
# Differentially expressed genes (DEG)
# Rank Rank Hypergeometric Overlap analysis (RRHO)
# Weighted Gene Co-expression Network Analysis (WGCNA)
###############################################################
source("code_carpentry/10_get_RNAseq_data.R")
source("code_analysis/RNAseq01_explore_sample_qualities.R")

source("code_analysis/RNAseq02_DESeq2_DEG.R")
source("code_analysis/RNAseq03_DESeq2_DEG_GOanalysis.R")
source("code_analysis/RNAseq04_limma_DEG_status_cort_interaction.R")
source("code_analysis/RNAseq05_limma_DEG_status_cort_interaction_GOanalysis.R")

# source("code_analysis/RNAseq06_RRHO.R")

# construct WGCNA for each liver and spleen tissue type 
source("code_analysis/RNAseq07_WGCNA_step1.R")
source("code_analysis/RNAseq08_WGCNA_step2.R")
source("code_analysis/RNAseq09_WGCNA_module_GO.R")
source("code_analysis/RNAseq10_WGCNA_step3.R")

# plot module eigengene against other behavioral and biological traits 
source("code_analysis/RNAseq11_WGCNA_ME_against_bioPCs.R")
source("code_analysis/RNAseq12_WGCNA_MEs_both_tissue_corrplot.R")


# Hypothesis driven gene set exploration 


# Differential coexpression analysis - CoDiNA 
source("code_analysis/RNAseq13_differential_coexpression_analysis.R")


# Genewalk 
source("code_carpentry/11_gather_GeneWalk_results.R")
source("code_analysis/RNAseq26_GeneWalk_kIN_overlap.R")



###############################################################
# 
# Figure codes 
#
###############################################################



