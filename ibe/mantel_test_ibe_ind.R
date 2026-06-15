library(vegan)
genetic <- read.table("genetic_distances_square_ind_pi_comp.tsv",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)
env <- read.table("canberra_distances_square_ind_comp.tsv",
                  header = TRUE,
                  row.names = 1,
                  sep = "\t",
                  check.names = FALSE)

geo <- read.table("harvesine_distances_square_ind_comp.tsv",
                  header = TRUE,
                  row.names = 1,
                  sep = "\t",
                  check.names = FALSE)
genetic_dist <- as.dist(genetic)
env_dist     <- as.dist(env)
geo_dist     <- as.dist(geo)
mantel_geo <- mantel(genetic_dist, geo_dist,
                     method = "pearson",
                     permutations = 9999)

mantel_geo

mantel_env <- mantel(genetic_dist, env_dist,
                     method = "pearson",
                     permutations = 9999)
mantel_env

partial_mantel <- mantel.partial(genetic_dist,
                                 env_dist,
                                 geo_dist,
                                 method = "pearson",
                                 permutations = 9999)

partial_mantel

partial_mantel_geo <- mantel.partial(genetic_dist,
                                     geo_dist,
                                     env_dist,
                                     method = "pearson",
                                     permutations = 9999)

partial_mantel_geo
