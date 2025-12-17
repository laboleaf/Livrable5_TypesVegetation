rm(list = ls())

# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(VegetationNunavik.data)
library(crew)
library(tidyverse)
library(vegan)

# Set target options:
tar_option_set(
  packages = c("VegetationNunavik.data", "tidyverse", "vegan",
               "mgcv", "Rtsne", "ggpubr", "ggrepel", "adespatial"), # Packages that your targets need for their tasks.
  controller = crew_controller_local(workers = 1),
  format = "qs"
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.


# Extraire les noms de formes de croissance pour faire les graphiq --------
formes_croissance <- Inventaires_Abondances_Tad |>
  retrait_sp_non_ident() |>
  filter(forme %in% c("Arbres", "Arbustes", "Herbes")) |>  ## TEMPORAIRE
  pull(forme.croissance) |>
  unique()

sous.regne.choix <-  c("Tout")

## Valeur de hauteur d'arbre pour la coupe
k_vec <- c(2, 10, 20)
groupes <- c(1,2)


# Replace the target list below with your own:
targets <- list(
    tar_target(Commu, import_commu(seuil = 2)),
    tar_target(Commu_forme, join_forme(Commu)),
    tar_target(NMDS, fit_NMDS(Commu = Commu_forme, k = 4)),
    tar_target(Commu_scores, extr_scores_NMDS(Commu = Commu_forme, NMDS = NMDS)),
    tar_target(plot_sp, graph_sp(Commu_scores)),
    tar_target(PAM, classif_pam(Commu_scores, kmin = 4, kmax = 30)),
    tar_target(Commu_pam, extr_cluster_pam(Commu_scores, PAM)),
    tar_target(Beta_flex, classif_beta_flex(Commu_pam)),
    tar_target(Commu_bflex, extr_cluster_bflex(Commu_pam, Beta_flex, kmin = 2, kmax = 30)),
    # tar_target(Sp_indic, calcul_sp_indic(Commu_clust, max.order = 1))
    tar_target(Commu_div, calcul_diversite(Commu_bflex))
)

list(targets)

