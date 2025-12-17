# Vider l'environnement
rm(list = ls())

## Importer les librairies
library(targets)
library(tidyverse)
library(VegetationNunavik.data)

# Importer les fonctions --------------------------------------------------
tar_source()

# Importation des données -------------------------------------------------
Commu <- import_commu(seuil = 2)

# Jonction des abondances par forme de croissance -------------------------
Commu_forme <- join_forme(tar_read(Commu))

# Estimation de la NMDS ---------------------------------------------------
NMDS <- fit_NMDS(tar_read(Commu), k = 4)
ordiplot(NMDS)

# Extraction des scores --------------------------------------------------
Commu_scores <- extr_scores_NMDS(tar_read(Commu_forme))

plot(Commu_scores$Codes_sp$NMDS1, Commu_scores$Codes_sp$NMDS2)

# Graphique des scores des espèces ----------------------------------------
plot_sp <- graph_sp(Commu = tar_read(Commu_scores))

plot_sp[[1]]

# Classification à l'aide de PAM ------------------------------------------
PAM <- classif_pam(Commu = tar_read(Commu_scores), kmin = 4, kmax = 6)

# Exraction de l'identité des clusters ------------------------------------
Commu_pam <- extr_cluster_pam(tar_read(Commu_scores), tar_read(PAM))

# Classification par béta flexible ----------------------------------------
Beta_flex <- classif_beta_flex(tar_read(Commu_pam))

# Exraction de l'identité des clusters ------------------------------------
Commu_bflex <- extr_cluster_bflex(tar_read(Commu_pam), tar_read(Beta_flex), kmin = 2, kmax = 30)

# Calcul des espèces indicatrices -----------------------------------------
Sp_indic <- calcul_sp_indic(tar_read(Commu_clust), max.order = 2)

# Calcul des indices de diversité -----------------------------------------
Diversite <- calcul_diversite(tar_read(Commu))

# Graphique des scores ----------------------------------------------------
plot_PCoA <- graph_scores(tar_read(Commu_scores))
plot_PCoA

# Stocker les noms d'espèces ----------------------------------------------
formes_croissance <- Inventaires_Abondances_Tad |>
  retrait_sp_non_ident() |>
  filter(forme %in% c("Arbres", "Arbustes", "Herbes")) |>  ## TEMPORAIRE
  pull(forme.croissance) |>
  unique()

noms_sp <- stock_nom_sp(Inventaires = Inventaires_Abondances_Tad,
                        forme_croissance = formes_croissance[2], nbre = 10)


# Calcul des moyennes pondérées des espèces -------------------------------
Sp_scores <- moy_pond_sp(Commu_Tout = tar_read(Commu), Commu_scores = tar_read(Commu_scores), noms_sp = tar_read(noms_sp))
summary(Sp_scores)
plot(Sp_scores$Dim1, Sp_scores$Dim2)




# Séparer la matrice de communauté ----------------------------------------
sep_matrice(tar_read(Commu_scores), groupe = 1)

# Jonction de la richesse en espèce --------------------------------------
Commu_rich <- joindre_richesse(tar_read(Commu_scores))


# Tracer un graphique de la richesse en espèces ---------------------------
graph_rich(Commu_rich, plot_tsne)

# Jonction des variables environnementales --------------------------------
Commu_env <- joindre_enviro(tar_read(Commu_rich))


# Tracer un graphique avec les variables environnementales ----------------
graph_env(tar_read(Commu_env), tar_read(plot_tsne))


# Estimations des gams ----------------------------------------------------
GAM <- estim_gam_sp(Commu_scores = Commu_scores, sp = "Vacc.uligi")
GAM


# Graphique des gams ------------------------------------------------------
gam_sp <- tar_read(gam_sp_Arbre)[[1]]
graph_gam_sp(gam_sp)

