#' """Code to run mixing models and trophic niche
#'    @author: Ryan James
#'    Date: 12/17/24"""

library(MixSIAR)
library(tidyverse)
library(furrr)
options(max.print = 6000000)

#### mix by time ####
#### predators fin clips ----
# load mixture files
mix = load_mix_data(file("data/pred_iso.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c('collection', 'Species'),
                    fac_random=c(F, F),
                    fac_nested=c(F, F),
                    cont_effects=NULL)

# load source files
source = load_source_data(file("data/source4.csv"),
                          source_factors=NULL,
                          conc_dep=TRUE,
                          data_type="means",
                          mix)

# load TEF
discr = load_discr_data(file("data/prTEF_pred4_2.csv"), mix)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags_pred = run_model(run="long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)

saveRDS(jags_pred, 'data/pred_mm4.rds')
#jags_pred = readRDS('data/pred_mm4.rds')
# Process JAGS output
output = list(summary_save = TRUE,
              summary_name = "data/predFin4_ss",
              sup_post = FALSE,
              plot_post_save_pdf = FALSE,
              plot_post_name = "lower_posterior_density",
              sup_pairs = FALSE,
              plot_pairs_save_pdf = FALSE,
              plot_pairs_name = "lower_pairs_plot",
              sup_xy = TRUE,
              plot_xy_save_pdf = FALSE,
              plot_xy_name = "lower_xy_plot",
              gelman = TRUE,
              heidel = FALSE,
              geweke = TRUE,
              diag_save = TRUE,
              diag_name = "data/predFin4_diag",
              indiv_effect = FALSE,
              plot_post_save_png = F,
              plot_pairs_save_png = FALSE,
              plot_xy_save_png = FALSE)

output_JAGS(jags_pred, mix, source, output)

# jags = readRDS('data/preds_mm.rds')

# extract values ----
jags = readRDS('data/pred_mm4.rds')

# MixSIAR function
mcmc.chains = jags$BUGSoutput$n.chains
N = mix$N
n.re = mix$n.re
n.effects = mix$n.effects
if (n.re == 1) {
  random_effects = ifelse(mix$FAC[[1]]$re, mix$FAC[[1]]$name, 
                           mix$FAC[[2]]$name)
}
if (n.re == 2) {
  random_effects = mix$factors
}
n.sources = source$n.sources
source_names = source$source_names
jags.mcmc = coda::as.mcmc(jags)
n.draws = length(jags$BUGSoutput$sims.list$p.global[,1])
fac2_lookup = list()
for (f1 in 1:mix$FAC[[1]]$levels) {
  fac2_lookup[[f1]] = unique(mix$FAC[[2]]$values[which(mix$FAC[[1]]$values == 
                                                          f1)])
}
ilr.both = array(NA, dim = c(n.draws, mix$FAC[[1]]$levels, 
                             mix$FAC[[2]]$levels, n.sources - 1))
p.both = array(NA, dim = c(n.draws, mix$FAC[[1]]$levels, 
                           mix$FAC[[2]]$levels, n.sources))
cross.both = array(data = NA, dim = c(n.draws, mix$FAC[[1]]$levels, 
                                      mix$FAC[[2]]$levels, n.sources, n.sources - 1))
e = matrix(rep(0, n.sources * (n.sources - 1)), nrow = n.sources, 
           ncol = (n.sources - 1))
for (i in 1:(n.sources - 1)) {
  e[, i] = exp(c(rep(sqrt(1/(i * (i + 1))), i), -sqrt(i/(i + 1)), rep(0, n.sources - i - 1)))
  e[, i] = e[, i]/sum(e[, i])
}
for (i in 1:n.draws) {
  for (f1 in 1:mix$FAC[[1]]$levels) {
    for (f2 in fac2_lookup[[f1]]) {
      for (src in 1:(n.sources - 1)) {
        ilr.both[i, f1, f2, src] = jags$BUGSoutput$sims.list$ilr.global[i,src] + 
          jags$BUGSoutput$sims.list$ilr.fac1[i,f1, src] + 
          jags$BUGSoutput$sims.list$ilr.fac2[i, f2, src]
        cross.both[i, f1, f2, , src] = (e[, src]^ilr.both[i, f1, f2, src])/sum(e[, src]^ilr.both[i,f1, f2, src])
      }
      for (src in 1:n.sources) {
        p.both[i, f1, f2, src] = prod(cross.both[i,f1, f2, src, ])
      }
      p.both[i, f1, f2, ] = p.both[i, f1, f2, ]/sum(p.both[i,f1, f2, ])
    }
  }
}

df_pred = as.data.frame(as.table(p.both)) |> 
  as_tibble() |> 
  mutate(Source = source$source_names[as.numeric(Var4)],
         Collection = mix$FAC[[1]]$labels[as.numeric(Var2)],
         Species = mix$FAC[[2]]$labels[as.numeric(Var3)],
         group = 'Predator',
         i = as.numeric(Var1)) |> 
  dplyr::select(Species, group, Collection, Source, Prop = Freq, i) |> 
  drop_na()


#### omnivores mixing model----
# load mix file
mix = load_mix_data(file("data/omn_iso.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c('collection', 'Species'),
                    fac_random=c(F, F),
                    fac_nested=c(F, F),
                    cont_effects=NULL)

# load source file
source = load_source_data(file("data/source4.csv"),
                          source_factors=NULL,
                          conc_dep=TRUE,
                          data_type="means",
                          mix)

#load TEF file
discr = load_discr_data(file("data/prTEF_omn4.csv"), mix)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags_omn = run_model(run="long", mix, source, discr, model_filename,
                     alpha.prior = 1, resid_err, process_err)

# saveRDS(jags_omn, 'data/omn_mm4.rds')
# jags_omn = readRDS('data/omn_mm4.rds')
# Process JAGS output
output = list(summary_save = TRUE,
              summary_name = "data/omn4_ss",
              sup_post = FALSE,
              plot_post_save_pdf = FALSE,
              plot_post_name = "lower_posterior_density",
              sup_pairs = FALSE,
              plot_pairs_save_pdf = FALSE,
              plot_pairs_name = "lower_pairs_plot",
              sup_xy = TRUE,
              plot_xy_save_pdf = FALSE,
              plot_xy_name = "lower_xy_plot",
              gelman = TRUE,
              heidel = FALSE,
              geweke = TRUE,
              diag_save = TRUE,
              diag_name = "data/omn4_diag",
              indiv_effect = FALSE,
              plot_post_save_png = F,
              plot_pairs_save_png = FALSE,
              plot_xy_save_png = FALSE)

output_JAGS(jags_omn, mix, source, output)

# jags = readRDS('data/omn_mm.rds')
# extract values ----
jags = readRDS('data/omn_mm4.rds')

# MixSIAR function
mcmc.chains = jags$BUGSoutput$n.chains
N = mix$N
n.re = mix$n.re
n.effects = mix$n.effects
if (n.re == 1) {
  random_effects = ifelse(mix$FAC[[1]]$re, mix$FAC[[1]]$name, 
                           mix$FAC[[2]]$name)
}
if (n.re == 2) {
  random_effects = mix$factors
}
n.sources = source$n.sources
source_names = source$source_names
jags.mcmc = coda::as.mcmc(jags)
n.draws = length(jags$BUGSoutput$sims.list$p.global[,1])
fac2_lookup = list()
for (f1 in 1:mix$FAC[[1]]$levels) {
  fac2_lookup[[f1]] = unique(mix$FAC[[2]]$values[which(mix$FAC[[1]]$values == 
                                                          f1)])
}
ilr.both = array(NA, dim = c(n.draws, mix$FAC[[1]]$levels, 
                             mix$FAC[[2]]$levels, n.sources - 1))
p.both = array(NA, dim = c(n.draws, mix$FAC[[1]]$levels, 
                           mix$FAC[[2]]$levels, n.sources))
cross.both = array(data = NA, dim = c(n.draws, mix$FAC[[1]]$levels, 
                                      mix$FAC[[2]]$levels, n.sources, n.sources - 1))
e = matrix(rep(0, n.sources * (n.sources - 1)), nrow = n.sources, 
           ncol = (n.sources - 1))
for (i in 1:(n.sources - 1)) {
  e[, i] = exp(c(rep(sqrt(1/(i * (i + 1))), i), -sqrt(i/(i + 1)), rep(0, n.sources - i - 1)))
  e[, i] = e[, i]/sum(e[, i])
}
for (i in 1:n.draws) {
  for (f1 in 1:mix$FAC[[1]]$levels) {
    for (f2 in fac2_lookup[[f1]]) {
      for (src in 1:(n.sources - 1)) {
        ilr.both[i, f1, f2, src] = jags$BUGSoutput$sims.list$ilr.global[i,src] + 
          jags$BUGSoutput$sims.list$ilr.fac1[i,f1, src] + 
          jags$BUGSoutput$sims.list$ilr.fac2[i, f2, src]
        cross.both[i, f1, f2, , src] = (e[, src]^ilr.both[i, f1, f2, src])/sum(e[, src]^ilr.both[i,f1, f2, src])
      }
      for (src in 1:n.sources) {
        p.both[i, f1, f2, src] = prod(cross.both[i,f1, f2, src, ])
      }
      p.both[i, f1, f2, ] = p.both[i, f1, f2, ]/sum(p.both[i,f1, f2, ])
    }
  }
}

df_omn = as.data.frame(as.table(p.both)) |> 
  as_tibble() |> 
  mutate(Source = source$source_names[as.numeric(Var4)],
         Collection = mix$FAC[[1]]$labels[as.numeric(Var2)],
         Species = mix$FAC[[2]]$labels[as.numeric(Var3)],
         group = 'Omnivore',
         i = as.numeric(Var1)) |> 
  dplyr::select(Species, group, Collection, Source, Prop = Freq, i) |> 
  drop_na()

# combine all data mixing model data----
df_all = bind_rows(df_pred, df_omn) |> 
  mutate(Collection = factor(Collection, levels = c('Dec-22', 'Aug-23', 'May-24')))

# hypervolume niche analysis -----
library(tidyverse)
library(hypervolume)

# subset mixing model data of species caught in all 3 collections ----
df_mm = df_all |> 
  filter(Species %in% c('Graysby', 'Red Hind',
                        'Yellowtail Snapper', 'Longspine Squirrelfish')) |> 
  pivot_wider(names_from = Source, values_from = Prop) |> 
  select(-i, -group) |> 
  group_by(Species, Collection) |> 
  nest()

# create replicates and z-score data across reps
reps = 100

set.seed(14)
df_z = df_mm |> 
  slice(rep(1:n(), each=reps))|> 
  group_by(Species, Collection) |> 
  mutate(i = row_number(),
         samp = map(data, \(data) slice_sample(data, n = 100))) |> 
  ungroup() |> 
  select(-data) |> 
  unnest(samp) |> 
  group_by(i) |> 
  mutate(across(Algae:POM, scale))


# generate hypervolumes
df = df_z |> 
  group_by(Species, Collection, i) |> 
  nest() |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = paste(Species, Collection, i,sep = '_'),
                                                     samples.per.point = 500,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)),
         centroid = map(hv, \(hv) get_centroid(hv)))

# saveRDS(df,'data/PRfish_hvs.rds')
# df = readRDS('data/PRfish_hvs.rds')
# occupancy
df_occ = df |> 
  ungroup() |> 
  select(Species, Collection, i, hv) |> 
  nest(.by = i) |> 
  mutate(hv_list= map(data, \(data) hypervolume_join(data$hv)),
         hv_occ = map2(data, hv_list, \(data, hv_list)
                       hypervolume_n_occupancy(hv_list, 
                                               method = "box", 
                                               classification = data$Collection, 
                                               box_density = 5000, 
                                               FUN = mean,
                                               verbose = F)),
         size = map(hv_occ, \(hv) get_volume(hv)),
         cent = map(hv_occ, \(hv) get_centroid_weighted(hv)),
         hv_df = map(hv_occ, \(data) hypervolume_to_data_frame(data)))

# write_rds(df_occ, 'data/hvOCC.rds')
#df_occ = read_rds('data/hvOCC.rds')

n_spp = length(unique(df$Species))

# calculate proportion of occupancy > 1 sp
oc_p = df_occ  |> 
  select(i,hv_df) |> 
  unnest(hv_df) |> 
  rename(Collection = Name, occupancy = ValueAtRandomPoints) |> 
  mutate(cat = if_else(occupancy > 1/n_spp, 'Multi species', 'Single species')) |> 
  group_by(i, Collection, cat) |> 
  summarize(n = n(), .groups = 'drop') |> 
  group_by(i, Collection) |> 
  mutate(prop = n/sum(n))

# calculate the size of each category
oc_s = df_occ |> 
  select(i, size) |> 
  unnest_longer(size, indices_to = 'Collection') |> 
  left_join(oc_p) |> 
  mutate(size_occ = prop*size)

# write_csv(oc_s, 'data/occ_size.csv')