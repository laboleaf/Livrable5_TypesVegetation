
# Importer la matrice de communautés --------------------------------------
import_commu <- function(seuil = 2){
  require(VegetationNunavik.data)
  require(tidyverse)
  require(sf)

  ## Importer les données d'abondance
  D_long <- Inventaires_Abondances_Tad |>
    filter(epithete != "sp." &
             complete.cases(famille) &
             famille != "NA") |>  # Retirer les espèces non identifiées
    group_by(parcelle) |>
    mutate(richesse = length(unique(sp.nom))) |>
    ungroup() |>
    filter(richesse > 1)
  ## Stocker les noms d'espèces et recoder les formes de croissance pour les graphiques futures
  Codes_sp <- D_long |>
    group_by(sp.nom, sp.code, famille, sous.regne, forme, forme.croissance, feuillage.vasc, pos.sporo.bryo, photobionte) |>
    reframe()

  ## Pivoter la matrice de communauté
  D_large <- D_long |>
    pivot_matrice_communaute(seuil = seuil, # Retirer les espèces dont occurence = 1
                             seuil_type = "occurence",
                             matrice = "Abondance", valeur = "Pourc") |>
    left_join(Inventaires_Enviro_Tad, by = c("campagne", "parcelle", "lat.dd", "lon.dd")) |>
    left_join(as_tibble(Inventaires_Spatial) |> select(-geometry),
              by = c("campagne", "parcelle", "type.releve"))

  D_list <- list(D_large = D_large, Codes_sp = Codes_sp)

  return(D_list)
}

# Matrice des forme de croissance -----------------------------------------
join_forme <- function(Commu){
  ## Synthétiser par forme de croissance
  Forme <- Inventaires_Abondances_Tad |>
    group_by(campagne, parcelle, lat.dd, lon.dd, type.releve, forme.croissance) |>
    summarise(couvert.pourc = sum(couvert.pourc)) |>
    pivot_wider(values_from = couvert.pourc, names_from = forme.croissance, values_fill = 0)

  Commu$D_large <- left_join(Commu$D_large, Forme)
  return(Commu)
}

# Estimer la NMDS ---------------------------------------------------------
fit_NMDS <- function(Commu, k){
  require(tidyverse)
  require(VegetationNunavik.data)
  require(vegan)

  # Extraire la matrice de communauté
  X <- Commu$D_large |>
    select(any_of(Commu$Codes_sp$sp.code))

  # Estimer la NMDS
  set.seed(999)
  NMDS <- metaMDS(comm = X,
                  distance = "hellinger",
                  k = k,
                  model = "global",
                  stress = 1,
                  weakties = F,
                  autotransform = F,
                  noshare = F,
                  try = 100, trymax = 500,
                  wascores = T, expand = T,
                  scaling = T, pc = T,
                  parallel = 4,
                  smin = 1e-4, sfgrmin = 1e-7, sratmax=0.999999)

  return(NMDS)
}

# Extraire les scores de la NMDS ------------------------------------------
extr_scores_NMDS <- function(Commu, NMDS){
  Scores <- scores(NMDS, tidy = T)

  Scores_sites <- Scores |> filter(score == "sites") |>
    select(-score, - label)

  Scores_sp <- Scores |> filter(score == "species") |>
    select(-score) |>
    rename(sp.code = label)

  Commu$D_large <- bind_cols(Commu$D_large, Scores_sites)
  Commu$Codes_sp <- left_join(Commu$Codes_sp, Scores_sp)

  return(Commu)
}

# Afficher les moyennes pondérées des espèces sélectionnées ---------------
graph_sp <- function(Commu){
  require(tidyverse)
  require(ggrepel)

  Commu$Codes_sp <- Commu$Codes_sp |>
    mutate(graph = case_when(
      forme.croissance == "Arbre" ~ "Arbre",
      forme.croissance == "Arbuste erige" ~ "Arbustes érigés, nains et prostrés",
      forme.croissance == "Arbuste nain" ~ "Arbustes érigés, nains et prostrés",
      forme.croissance == "Arbuste prostre" ~ "Arbustes érigés, nains et prostrés",
      famille == "Cyperaceae" ~ "Graminoïdes - Cyperaceae",
      forme.croissance == "Graminoide" & famille != "Cyperaceae" ~ "Graminoïdes - Autres",
      famille %in% c("Asteraceae", "Caryophyllaceae" , "Orobanchaceae", "Ranunculaceae", "Fabaceae") ~
        "Phorbes à fleurs - familles les plus diversifiées",
      forme.croissance == "Phorbe a fleurs" &
        famille %nin% c("Asteraceae", "Caryophyllaceae" , "Orobanchaceae", "Ranunculaceae", "Fabaceae") ~
        "Phorbes à fleurs - familles moins diversifiées",
      forme.croissance == "Phorbe a spores" ~ "Phorbes à spores",
      famille %in% c("Polytrichaceae", "Dicranaceae", "Bryaceae", "Mniaceae") ~ "Mousses acrocarpes - familles les plus diversifiées",
      famille %nin% c("Polytrichaceae", "Dicranaceae", "Bryaceae", "Mniaceae") & pos.sporo.bryo == "Acrocarpe" ~ "Mousses acrocarpes - familles les moins diversifiées",
      pos.sporo.bryo == "Pleurocarpe" ~ "Mousses pleurocarpes",
      forme.croissance == "Hepatique a feuilles" ~ "Hépatiques à feuilles",
      forme.croissance == "Hepatique a thalle" ~ "Hépatiques à thalle",
      .default = forme.croissance
    ))

  p_sites <- Commu$D_large |>
    ggplot(aes(x = NMDS1, y = NMDS2)) +
    geom_point(alpha = 0.4) +
    labs(x = "NMDS1", y = "NMDS2") +
    theme_minimal()

  ## Tracer le graphique
  p_sp <- unique(Commu$Codes_sp$graph) |>
    set_names() |>
    map(
      function(x){
        Res <- p_sites +
          geom_point(data = filter(Commu$Codes_sp, graph == x),
                     aes(NMDS1, NMDS2, color = forme.croissance), shape = "cross") +
          geom_text_repel(data = filter(Commu$Codes_sp, graph == x),
                          aes(NMDS1, NMDS2, label = sp.code, color = forme.croissance),
                          vjust = 0) +
          ggtitle(paste(x))
      }
    )

  return(p_sp)
}


# Fonction pour calculer la contribution locale à la diversité béta --------------------------------------
calcul_lcbd <- function(Commu){
  require(tidyverse)
  require(adespatial)

  X <- Commu$D_large |>
  select(any_of(Commu$Codes_sp$sp.code))

  Beta <- beta.div(X, method = "hellinger", sqrt.D = F)

  Commu$D_large <-  Commu$D_large |>
    mutate(lcbd = Beta$LCBD)

  Commu$Codes_sp <- Commu$Codes_sp |>
    left_join(
      tibble(sp.code = names(Beta$SCBD), scbd = Beta$SCBD)
    )
}

# Classifier avec PAM -----------------------------------------------------
classif_pam <- function(Commu, kmin, kmax){
  require(tidyverse)
  require(cluster)
  require(purrr)

  # Extraire les valeurs d'axes de nmds
  Scores <- Commu$D_large |>
    select(starts_with("NMDS"))
  # Calculer la matrice de distance
  D_scores <- dist(Scores)

  # Partitionnement
  PAM <- kmin:kmax |>
    set_names() |>
    map(
      \(x)     pam(D_scores, k = x, diss = T,
                   metric = "euclidean",
                   medoids = "random", nstart = 100,
      )
    )


  return(PAM)
}

# Extraire les clusters de PAM ---------------------------------------------------
extr_cluster_pam <- function(Commu, PAM){
  require(tidyverse)
  # Nombre de partitionnements
  n_pam <- length(PAM)

  Clusters_pam <- 1:n_pam |>
    map(
      \(x){
        df <- data.frame(A = PAM[[x]]$clustering)
        colnames(df) <- paste0("cluster.pam.", names(PAM)[x])
        return(df)
      }
    ) |>
    list_cbind()

  Commu$D_large <- bind_cols(
    Commu$D_large, Clusters_pam
  )

  return(Commu)
}


# Classification avec la méthode Béta flexible ----------------------------
classif_beta_flex <- function(Commu){
  require(tidyverse)
  require(vegan)
  require(cluster)

  X <- Commu$D_large |>
    select(any_of(Commu$Codes_sp$sp.code))

  Hell <- vegdist(X, method = "hellinger")

  beta_flex <- agnes(x = Hell, diss = T, method = 'flexible', par.method = 0.625)

  return(beta_flex)
}

# Extraire les clusters de béta flexible ---------------------------------------------------
extr_cluster_bflex <- function(Commu, Beta_flex, kmin, kmax){
  require(tidyverse)
  require(cluster)

  Clusters_bflex <- kmin:kmax |>
    map(
      \(x){
        df <- data.frame(A = cutree(Beta_flex, x))
        colnames(df) <- paste0("cluster.bflex.", x)
        return(df)
      }
    ) |>
    list_cbind()

  Commu$D_large <- bind_cols(
    Commu$D_large, Clusters_bflex
  )

  return(Commu)
}

# Calculer les espèces indicatrices ---------------------------------------
calcul_sp_indic <- function(Commu, max.order){
  require(tidyverse)
  require(indicspecies)
  # Extraire la matrice de communauté
  X <- Commu$D_large |>
    select(any_of(unique(Commu$Codes_sp$sp.code)))
  # Extraire les clusters
  clusters <- Commu$D_large$cluster_pam
  # Calculer et tester les valeurs indicatrices
  sp_indic <- multipatt(X, clusters,
                        func = "IndVal.g",
                        duleg = F,
                        restcomb = NULL,
                        min.order = 1, max.order = max.order,
                        control = how(nperm=999))

  return(sp_indic)
}


# Créer un tableau des formes de croissance -------------------------------
table_formes_croissance <- function(Commu, var_groupe){
  require(tidyverse)
  Table_formes_croissance <- Commu$D_large |>
    # Sélectionnner les colonnes et passer au format long
    select(any_of(var_groupe), parcelle, `Arbuste erige`:`Hepatique a thalle`) |>
    pivot_longer(-parcelle:-any_of(var_groupe),
                 values_to = "couvert",
                 names_to = "forme.croissance") |>
    # Calculer le nombre de parcelles par cluster
    group_by(across(any_of(var_groupe))) |>
    mutate(n.parc = length(unique(parcelle))) |>
    # Calculer les occurrences et abondances moyennes
    group_by(across(any_of(c(var_groupe, "forme.croissance")))) |>
    reframe(occ = sum(ifelse(couvert > 0, 1, 0)),
            n.parc = unique(n.parc),
            occ.pourc = occ/n.parc*100,
            couvert.med = median(couvert),
            couvert.min = min(couvert),
            couvert.max = max(couvert)) |>
    # Créer les classes d'occurence
    mutate(occ.class = case_when(
      occ.pourc < 20 ~ "I",
      occ.pourc >= 20 & occ.pourc < 40 ~ "II",
      occ.pourc >= 40 & occ.pourc < 60 ~ "III",
      occ.pourc >= 60 & occ.pourc < 80 ~ "IV",
      occ.pourc >= 80 ~ "V"
    ),
    couvert = paste0(couvert.med, " (", couvert.min, "-",
                     couvert.max, ")")) |>
    # Réordonner les colonnes
    relocate(forme.croissance, occ.class, couvert, any_of(var_groupe)) |>
    # Réordonner les noms de formes de croissances
    mutate(forme.croissance = factor(
      forme.croissance,
      levels = c("Arbre", "Arbuste erige", "Arbuste nain", "Arbuste prostre",
                 "Graminoide", "Phorbe a fleurs", "Phorbe a spores",
                 "Bryophyte s.s.", "Sphaigne",
                 "Hepatique a feuilles", "Hepatique a thalle",
                 "Lichen fruticuleux", "Lichen foliace",
                 "Lichen crustace", "Lichen squamuleux")
    )) |>
    ungroup()

  return(Table_formes_croissance)
}


# Créer un tableau résumé des espèces -------------------------------------
table_sp <- function(Commu, var_groupe){
  Table_sp <- Commu$D_large |>
    # Sélectionner les colonnes et passer au format long
    select(parcelle, any_of(var_groupe),
           any_of(Commu$Codes_sp$sp.code)) |>
    pivot_longer(-any_of(var_groupe):-parcelle,
                 names_to = "sp.code", values_to = "couvert") |>
    # Calculer le nombre de parcelles par cluster
    group_by(across(any_of(var_groupe))) |>
    mutate(n.parc = length(unique(parcelle))) |>
    # Calculer les occurrences et abondances moyennes
    filter(couvert > 0) |>
    group_by(across(any_of(c(var_groupe, "sp.code")))) |>
    summarise(occ = length(unique(parcelle)),
              n.parc = unique(n.parc),
              occ.pourc = occ/n.parc*100,
              couvert.med = median(couvert),
              couvert.min = min(couvert),
              couvert.max = max(couvert)) |>
    # Créer les classes d'occurence
    mutate(occ.class = case_when(
      occ.pourc < 20 ~ "I",
      occ.pourc >= 20 & occ.pourc < 40 ~ "II",
      occ.pourc >= 40 & occ.pourc < 60 ~ "III",
      occ.pourc >= 60 & occ.pourc < 80 ~ "IV",
      occ.pourc >= 80 ~ "V"
    ),
    couvert = paste0(couvert.med, " (", couvert.min, "-",
                     couvert.max, ")")) |>
    left_join(Commu$Codes_sp |>
                select(sp.code, sp.nom, forme.croissance, sous.regne)) |>
    # Réordonner les colonnes
    relocate(sp.nom, occ.class, couvert, any_of(var_groupe)) |>
    # Réordonner les noms de sous règnes
    mutate(sous.regne = factor(
      sous.regne,
      levels = c("Vasculaires", "Bryophytes", "Lichens")
    )) |>
    ungroup()

  return(Table_sp)
}


# Tableau présentant l'environnement --------------------------------------
table_enviro <- function(Commu, var_groupe){
  Table_enviro <- Commu$D_large |>
    group_by(across(any_of(var_groupe))) |>
    add_count() |>
    add_count(drain.class.1, name = "drain.class.1.count") |>
    add_count(mat.org.type, name = "mat.org.type.count") |>
    add_count(depot.surf.text.1.nom, name = "depot.surf.text.1.nom.count") |>
    add_count(situ.topo.nom, name = "situ.topo.nom.count")  |>
    add_count(depot.surf.melccfp.group.nom, name = "depot.surf.melccfp.group.nom.count") |>
    add_count(expo.class, name = "expo.class.count") |>
    mutate(across(ends_with("count"), .fns = ~round(.x/n, 2)*100)) |>
    mutate(
      drainage = paste0(drain.class.1, ": ",
                        drain.class.1.count),
      matiere.organique = paste0(mat.org.type, ": ",
                                 mat.org.type.count),
      texture = paste0(depot.surf.text.1.nom, ": ",
                       depot.surf.text.1.nom.count),
      topographie = paste0(situ.topo.nom, ": ",
                           situ.topo.nom.count),
      depots = paste0(depot.surf.melccfp.group.nom, ": ",
                      depot.surf.melccfp.group.nom.count),
      exposition = paste0(expo.class, ": ",
                      expo.class.count)
    ) |>
    arrange(desc(drain.class.1.count)) |>
    mutate(drainage = paste0(unique(drainage), collapse = "; ")) |>
    arrange(desc(mat.org.type.count)) |>
    mutate(matiere.organique = paste0(unique(matiere.organique), collapse = "; ")) |>
    arrange(desc(depot.surf.text.1.nom.count)) |>
    mutate(texture = paste0(unique(texture), collapse = "; ")) |>
    arrange(desc(situ.topo.nom.count)) |>
    mutate(topographie = paste0(unique(topographie), collapse = "; ")) |>
    arrange(desc(depot.surf.melccfp.group.nom.count)) |>
    mutate(depots = paste0(unique(depots), collapse = "; ")) |>
    arrange(desc(expo.class.count)) |>
    mutate(exposition = paste0(unique(exposition), collapse = "; ")) |>
    summarise(across(c(drainage, matiere.organique, texture, topographie, depots, exposition), unique)) |>
    pivot_longer(drainage:exposition, names_to = "Variable", values_to = "Comptes")

  return(Table_enviro)
}



# Calculer les indices de diversité ------------------------------------------
calcul_diversite <- function(Commu){
  require(VegetationNunavik.data)
  require(tidyverse)
  require(sf)
  require(vegan)

  ## Importer les données d'abondance
  D_long <- Inventaires_Abondances_Tad |>
    filter(epithete != "sp." &
             complete.cases(famille) &
             famille != "NA") |>  # Retirer les espèces non identifiées
    group_by(parcelle) |>
    mutate(richesse = length(unique(sp.nom))) |>
    ungroup() |>
    filter(richesse > 1)
  ## Stocker les noms d'espèces et recoder les formes de croissance pour les graphiques futures
  Codes_sp <- D_long |>
    group_by(sp.nom, sp.code, famille, sous.regne, forme, forme.croissance, feuillage.vasc, pos.sporo.bryo, photobionte) |>
    reframe()

  ## Pivoter la matrice de communauté
  D_large <- D_long |>
    pivot_matrice_communaute(seuil = 0, # Retirer les espèces dont occurence = 1
                             seuil_type = "occurence",
                             matrice = "Abondance", valeur = "Pourc")
  # Richesse totale
  Diversite <- data.frame(richesse = D_large |>
    select(any_of(D_long$sp.code)) |>
    specnumber()
  )
  # Richesse des plantes vasculaires
  Diversite$richesse.vasc <- D_large |>
    select(
      any_of(filter(D_long, sous.regne == "Vasculaires") |> pull(sp.code))
           ) |>
    specnumber()
  # Richesse des bryophytes
  Diversite$richesse.bryo <- D_large |>
    select(
      any_of(filter(D_long, sous.regne == "Bryophytes") |> pull(sp.code))
    ) |>
    specnumber()
  # Richesse des lichens
  Diversite$richesse.lichens <- D_large |>
    select(
      any_of(filter(D_long, sous.regne == "Lichens") |> pull(sp.code))
    ) |>
    specnumber()
  # Indice de Shannon total
  Diversite$shannon <- D_large |>
    select(any_of(D_long$sp.code)) |>
    diversity(index = "shannon")
  # Indice de Shannon des plantes vasculaires
  Diversite$shannon.vasc <- D_large |>
    select(filter(D_long, sous.regne == "Vasculaires") |> pull(sp.code)) |>
    diversity(index = "shannon")
  # Indice de Shannon des bryophytes
  Diversite$shannon.bryo <- D_large |>
    select(filter(D_long, sous.regne == "Bryophytes") |> pull(sp.code)) |>
    diversity(index = "shannon")
  # Indice de Shannon des lichens
  Diversite$shannon.lichens <- D_large |>
    select(filter(D_long, sous.regne == "Lichens") |> pull(sp.code)) |>
    diversity(index = "shannon")
  # Calculer l'équitabilité de Pielou
  Diversite <- Diversite |>
    mutate(pielou = shannon/log(richesse),
           pielou.vasc = shannon.vasc/log(richesse.vasc),
           pielou.bryo = shannon.bryo/log(richesse.bryo),
           pielou.lichens = shannon.lichens/log(richesse.lichens),
           campagne = D_large$campagne,
           parcelle = D_large$parcelle)

  Commu$D_large <- left_join(Commu$D_large, Diversite)

  return(Commu)
}

# # Tracer un graphique avec la richesse en espèces -------------------------
# graph_rich <- function(Commu_rich, plot_tsne){
#   plot_tsne +
#     geom_point(data = Commu_rich, aes(size = richesse, color = richesse)) +
#     scale_colour_viridis_c() +
#     scale_size_continuous(range = c(0.5, 6))
#
# }

# # Joindre les variables environnementales ------------------------------------------
# joindre_enviro <- function(Commu_rich){
#   require(tidyverse)
#   require(VegetationNunavik.data)
#
#   Env <- Inventaires_Enviro_Tad
#
#   Res <- left_join(Commu_rich, Env)
#
#   return(Res)
# }
#

#
#
# # Estimer les gams des abondances d'espèces -------------------------------
# estim_gam_sp <- function(Commu_scores, sp){
#   require(mgcv)
#
#   ## Créer la formule
#   form <- as.formula(
#     paste0(sp, " ~ tsne1 * tsne2")
#   )
#
#   ## Estimer le modèle
#   Res <- gam(formula = form, data = Commu_scores, family = binomial)
#
#   return(Res)
# }
#
# # Ajout les surface des espèces -------------------------------------------
# graph_gam_sp <- function(gam_sp){
#   require(mgcv)
#   require(tidyverse)
#   require(modelr)
#   require(metR)
#
#   ## Extraire les données
#   D <- gam_sp$model
#
#   ## Créer la surface de prédictions
#   Pred <- D |>
#     data_grid(tsne1 = seq_range(tsne1, 50), tsne2 = seq_range(tsne2, 50)) |>
#     add_predictions(gam_sp, type = "response")
#
#   ## Créer le graphique
#   D |> ggplot(aes(tsne1, tsne2)) +
#     geom_tile(data = Pred, aes(fill = pred)) +
#     geom_contour(data = Pred, aes(z = pred), color = "white") +
#     geom_text_contour(data = Pred, aes(z = pred), color = "white") +
#     geom_point() +
#     ggtitle(colnames(D)[1]) +
#     scale_fill_distiller(palette = "YlGnBu") +
#     theme_minimal()
# }
#
#
# # Combiner les graphiques d'espèces ---------------------------------------
# combiner_graph_gam <- function(plot_gam){
#   require(ggpubr)
#   Res <- ggarrange(plotlist = plot_gam,
#                    ncol = 2, nrow = 2)
# }




