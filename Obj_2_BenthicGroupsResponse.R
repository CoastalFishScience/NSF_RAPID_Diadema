
library(tidyverse)
# library(hypervolume)
library(viridis)
library(ggpubr)


#Models libraries
library(mgcv)
library(glmmTMB)
# library(gam)
# library(gamm4)
# library(nlme)

#Model check libraries
library(visreg)
library(DHARMa)
library(performance)
#library(lsmeans) - predecessor for emmeans
library(emmeans)
library(MuMIn)
library(ggeffects)

library(multcomp)
library(multcompView)

#Community Analysis
library(vegan)
library(RVAideMemoire)

#trajectory analysis
library(ecotraj)

# Load clean datasets -----------------------------------------------------

pcover = read_csv("./data/pcover.csv") |> 
  mutate(period = factor(period, levels =c("baseline", "dieoff", "after 6mo", "after 12mo")),location=as.factor(location), site=as.numeric(site), plot=as.numeric(plot))

urchins = read_csv("./data/urchins_clean.csv") |> 
  drop_na(location)|> 
  mutate(period = factor(period, levels = c("baseline", "dieoff", "after 6mo", "after 12mo")))

urchins_counts = read_csv("./data/urchin_counts.csv") |> 
  drop_na(location) |> 
  mutate(period = factor(period, levels = c("baseline", "dieoff", "after 6mo", "after 12mo")))

community_matrix<-read_csv("C:/Users/nicor/OneDrive - Florida International University/Documents/Git/RAPID_benCom/Data/community.csv") |> 
  mutate(site = as.numeric(site),plot = as.numeric(plot))

nmds <- readRDS("C:/Users/nicor/OneDrive - Florida International University/Documents/Git/RAPID_benCom/rds/nmds_results.rds")

nmds_coords <- as.data.frame(nmds$points) 

colnames(nmds_coords) <- c("NMDS1", "NMDS2")  

nmds_coords$period <- community_matrix$period   
nmds_coords$location <- community_matrix$location  

stress_value <- nmds$stress   
nmds_coords$period <- factor(nmds_coords$period, levels = c("baseline", "dieoff", "after 6mo", "after 12mo"))

nmdscom = ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = period)) +  
  geom_point(size = 3, alpha = 0.7) +  
  labs(title = "NMDS of Community Data", x = "NMDS 1", y = "NMDS 2", color = 
         "Periods") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = .90,
           label = paste("Stress =", round(stress_value, 4)), 
           size = 5, fontface = "bold") + 
  theme_minimal() +   
  theme(axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), # Increase legend title size
        legend.text = element_text(size = 12)) +
  scale_color_viridis_d(option = "turbo", 
                        labels = c(
                          "baseline" = "Baseline",
                          "dieoff" = "Die-off",
                          "after 6mo" = "After 4 Months",
                          "after 12mo" = "After 12 Months"
                        )
  )


nmdscom


nmds_coordsM = nmds_coords |> 
  filter(location == "Maguey")

nmds_coordsT = nmds_coords |>
  filter(location == "Tampico")

nmdslocM = ggplot(nmds_coordsM, aes(x = NMDS1, y = NMDS2, color = period)) + 
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Maguey",x = "NMDS 1", y = "NMDS 2", color = "Periods") +
  theme_minimal() +   
  theme(axis.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  scale_color_viridis_d(option = "turbo", 
                        labels = c(
                          "baseline" = "Baseline",
                          "dieoff" = "Die-off",
                          "after 6mo" = "After 4 Months",
                          "after 12mo" = "After 12 Months"
                        )
  )

nmdslocM

nmdslocT = ggplot(nmds_coordsT, aes(x = NMDS1, y = NMDS2, color = period)) + 
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Tampico",x = "NMDS 1", y = "NMDS 2", color = "Periods") +
  theme_minimal() +   
  theme(axis.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"), # Increase legend title size
        legend.text = element_text(size = 12),
        legend.position = "none") + 
  scale_color_viridis_d(option = "turbo", 
                        labels = c(
                          "baseline" = "Baseline",
                          "dieoff" = "Die-off",
                          "after 6mo" = "After 4 Months",
                          "after 12mo" = "After 12 Months"
                        )
  )

nmdslocT

community_matrix <- pcover |>
  mutate(plot = as.numeric(plot), site = as.numeric(site),location = as.factor(location)) |>
  filter(category %in% c("coral", "gorgonian", "sponge", "algae","seagrass")) |>
  dplyr::select(location, site, plot, period, category, percentcover) |>
  pivot_wider(
    names_from = category,
    values_from = percentcover,
    values_fill = 0
  )


distance_matrix <- vegdist(community_matrix[5:9], method = "bray") 

permanova_results <- adonis2(distance_matrix ~ location * period, 
                             data = community_matrix, 
                             permutations = 999,
                             strata = community_matrix$site, 
                             by ="terms"
) 

permanova_results


# Perform pairwise comparisons for 'location'
pairwise_results_location <- pairwise.adonis2(distance_matrix ~ location, data = community_matrix, nperm = 999)

# Perform pairwise comparisons for 'period'
pairwise_results_period <- pairwise.adonis2(distance_matrix ~ period, data = community_matrix, nperm = 999, p.method = "bonferroni")

# Combine results into a list or data frame
combined_pairwise_results <- list(Location = pairwise_results_location, Period = pairwise_results_period)
combined_pairwise_results





# Univariate Analysis -----------------------------------------------------



# Benthic Cover Analysis --------------------------------------------------

coral = pcover |> 
  filter(category %in% c("coral"))

algae = pcover |> 
  filter(category %in% c("algae")) 

sponge = pcover |> 
  filter(category %in% c("sponge")) 

## Exploratory analysis ####

ggplot(coral, aes(period, percentcover, fill=period)) +
  geom_boxplot()+
  facet_wrap(~location)+
  scale_fill_viridis_d()

ggplot(algae, aes(period, percentcover, fill=period)) +
  geom_violin()+
  facet_wrap(~location)+
  scale_fill_viridis_d()

ggplot(sponge, aes(period, percentcover, fill=period)) +
  geom_violin()+
  facet_wrap(~location)+
  scale_fill_viridis_d()



##Beta models####

### Corals ####

coralbeta = glmmTMB(percentcover ~ location*period + (1|site), 
                    family = beta_family(), data = coral)

check_model(coralbeta)
summary(coralbeta)
performance(coralbeta)


car::Anova(coralbeta)

#No interaction effect was seen between period*location so we will change to additive

coralbeta2 = glmmTMB(percentcover ~ location + period + (1|site), 
                     family = beta_family(), data = coral)

check_model(coralbeta2)
summary(coralbeta2)
performance(coralbeta2)
car::Anova(coralbeta2)

compare_performance(coralbeta, coralbeta2)

#coralbeta2 is the better model. Tampico is sig different and so are perioddieoff and periodafter 12 mo


gg.coralperiod <-ggpredict(coralbeta2, terms = c("period")) |> 
  rename(period = x, cover = predicted)


model_means_coralperiodbeta2 = emmeans(coralbeta2, specs = ~ period)

model_means_cld_coralperiodbeta2 = cld(object = model_means_coralperiodbeta2,
                                       Letters = letters,
                                       alpha = 0.05)

gg.coralperiod2 <- merge(gg.coralperiod, model_means_cld_coralperiodbeta2, by = c("period"))

cor = ggplot(gg.coralperiod2, aes(period, cover*100, colour = period))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position = position_dodge(width = 0.5))+
  geom_text(aes(label = .group,y=1.3 # Set y based on location
  ), show.legend = F, size = 6,
  position = position_dodge(width = 0.5))+
  labs(x = NULL, y = "% Coral Cover",color= "period")+
  scale_x_discrete(
    labels = c(
      "baseline" = "Baseline",
      "dieoff" = "Die-off",
      "after 6mo" = "After 4 Months",
      "after 12mo" = "After 12 Months"
    )
  ) +
  scale_color_viridis_d(option = 'turbo', 
                        labels = c(
                          "baseline" = "Baseline",
                          "dieoff" = "Die-off",
                          "after 6mo" = "After 4 Months",
                          "after 12mo" = "After 12 Months"
                        ))+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),legend.position = "none")

gg.coralloc <- ggpredict(coralbeta2, terms = c("location")) |> 
  rename(location = x, cover = predicted)

model_means_corallocbeta2 = emmeans(coralbeta2, specs = ~ location)

model_means_cld_corallocbeta2 = cld(object = model_means_corallocbeta2,
                                    Letters = letters,
                                    alpha = 0.05)

gg.coralloc2 <- merge(gg.coralloc, model_means_cld_corallocbeta2, by = "location")


cor2 = ggplot(gg.coralloc2, aes(location, cover*100, colour = location))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position = position_dodge(width = 0.5))+
  geom_text(aes(label = .group, y = ifelse(location == "Tampico", 3, 3) # Set y based on location
  ), show.legend = F, size = 6,
  position = position_dodge(width = 0.5))+
  labs(x = NULL, y = "% Coral Cover",color= "location")+
  scale_x_discrete(
    labels = c(
      "baseline" = "Baseline",
      "dieoff" = "Die-off",
      "after 6mo" = "After 4 Months",
      "after 12mo" = "After 12 Months"
    )
  ) +
  scale_color_viridis_d(option = 'turbo', 
                        labels = c(
                          "baseline" = "Baseline",
                          "dieoff" = "Die-off",
                          "after 6mo" = "After 4 Months",
                          "after 12mo" = "After 12 Months"
                        ))+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),legend.position = "none")


#Average coral covers

avg_coral = coral |> 
  #filter(period == "baseline") |>
  group_by(location, period) |> 
  summarize(mean_cover = mean(percentcover)*100, sd_cover = sd(percentcover)*100)

#Shapiro test for nornamlity
coral |> 
  filter(period == "baseline") |> 
  group_by(location) |> 
  summarise(shapiro_p = shapiro.test(percentcover)$p.value)

#Data is not normal so using wilcox rank test 

coral |> 
  filter(period == "baseline") |> 
  wilcox.test(percentcover ~ location, data = _)

#Percent increase
percent_increase <- avg_coral |> 
  dplyr::select(1:3) |> 
  filter(period %in% c("baseline", "after 12mo")) |> 
  pivot_wider(names_from = period, values_from = mean_cover, names_prefix = "cover_") |> 
  mutate(percent_increase = ((`cover_after 12mo` - cover_baseline) / cover_baseline) * 100) |> 
  summarise(location = location, average_increase = mean(percent_increase, na.rm = TRUE))

# Calculate the overall average increase across both locations
mean(percent_increase$average_increase, na.rm = TRUE)


###Algae####

algaebeta = glmmTMB(percentcover ~ location*period + (1|site), 
                    family = beta_family(), data = algae)

check_model(algaebeta)
summary(algaebeta)
performance(algaebeta)

car::Anova(algaebeta)

#Anova shows no interaction effect, so we will remove it

algaebeta2 = glmmTMB(percentcover ~ location + period + (1|site), 
                     family = beta_family(), data = algae)

check_model(algaebeta2)
summary(algaebeta2)
performance(algaebeta2)
car::Anova(algaebeta2)

compare_performance(algaebeta, algaebeta2)

#algaebeta2 is the better model. TEverything is showing significance

gg.algaeperiod = ggpredict(algaebeta2, terms = c("period")) |> 
  rename(period = x, cover = predicted)

gg.algaeloc = ggpredict(algaebeta2, terms = c("location")) |> 
  rename(location = x, cover = predicted)

model_means_algaelocbeta2 = emmeans(algaebeta2, specs = ~ location)

model_means_algaeperiodbeta2 = emmeans(algaebeta2, specs = ~ period)

model_means_cld_algaelocbeta2 = cld(object = model_means_algaelocbeta2,
                                    Letters = letters,
                                    alpha = 0.05)

model_means_cld_algaeperiodbeta2 = cld(object = model_means_algaeperiodbeta2,
                                       Letters = letters,
                                       alpha = 0.05)

gg.algaeperiod2 <- merge(gg.algaeperiod, model_means_cld_algaeperiodbeta2, by = "period")

gg.algaeloc2 <- merge(gg.algaeloc, model_means_cld_algaelocbeta2, by = "location")



alg <- ggplot(gg.algaeperiod2, aes(period, cover*100, colour = period))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position =position_dodge(width = 0.5))+
  labs(x = NULL, y = "% Algae Cover")+
  scale_x_discrete(
    labels = c(
      "baseline" = "Baseline",
      "dieoff" = "Die-off",
      "after 6mo" = "After 4 Months",
      "after 12mo" = "After 12 Months"
    )
  ) +
  scale_color_viridis_d(option = 'turbo', 
                        labels = c(
                          "baseline" = "Baseline",
                          "dieoff" = "Die-off",
                          "after 6mo" = "After 4 Months",
                          "after 12mo" = "After 12 Months"
                        ))+
  geom_text(aes(label=.group,  y=30,
  ), size=6,position = position_dodge(width = 0.5))+
  scale_color_viridis(option = 'turbo', discrete = T)+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12))





alg2 <- ggplot(gg.algaeloc2, aes(location, cover*100, colour = location))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position =position_dodge(width = 0.5))+
  labs(x = NULL, y = "% Algae Cover")+
  geom_text(aes(label=.group,  y=42,
  ), size=6,position = position_dodge(width = 0.5))+
  scale_color_viridis(option = 'turbo', discrete = T)+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12))


#Average algae covers

avg_algae = algae |> 
  #filter(period == "baseline") |>
  group_by(location, period) |> 
  summarize(mean_cover = mean(percentcover)*100, sd_cover = sd(percentcover)*100)

#Shapiro test for nornamlity
algae |> 
  filter(period == "baseline") |> 
  group_by(location) |> 
  summarise(shapiro_p = shapiro.test(percentcover)$p.value)

#Maguey is not normal so using wilcox rank test 

algae |> 
  filter(period == "baseline") |> 
  wilcox.test(percentcover ~ location, data = _)

#Percent increase
percent_increase_algae <- avg_algae |> 
  dplyr::select(1:3) |> 
  filter(period %in% c("baseline", "after 12mo")) |> 
  pivot_wider(names_from = period, values_from = mean_cover, names_prefix = "cover_") |> 
  mutate(percent_increase = ((`cover_after 12mo` - cover_baseline) / cover_baseline) * 100) |> 
  summarise(location = location, average_increase = mean(percent_increase, na.rm = TRUE))

# Calculate the overall average increase across both locations
mean(percent_increase_algae$average_increase, na.rm = TRUE)


###Sponges####

spongebeta = glmmTMB(percentcover ~ location*period + (1|site), 
                     family = beta_family(), data = sponge)

check_model(spongebeta)
summary(spongebeta)
performance(spongebeta)

car::Anova(spongebeta)

#Anova shows an interaction effect

spongebeta2 = glmmTMB(percentcover ~ location + period + (1|site), 
                      family = beta_family(), data = sponge)

check_model(spongebeta2)
summary(spongebeta2)
performance(spongebeta2)

compare_performance(spongebeta, spongebeta2)

#spongebeta is the better model. TEverything is showing significance

sponge_periods = ggpredict(spongebeta, terms = c("period")) 

plot(sponge_periods)

sponge_locations = ggpredict(spongebeta, terms = c("location"))

plot(sponge_locations)


#Final model for sponge


gg.spongebeta = ggpredict(spongebeta, terms = c("period", "location")) |> 
  rename(period = x, location = group, cover = predicted)

model_means_spongebeta = emmeans(spongebeta, specs = ~ period | location)


model_means_cld_spongebeta = cld(object = model_means_spongebeta,
                                 Letters = letters,
                                 alpha = 0.05)

gg.spongebeta2 <- merge(gg.spongebeta, model_means_cld_spongebeta, by = c("period", "location"))

spo = ggplot(gg.spongebeta2, aes(period, cover*100, colour = location))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position =position_dodge(width = 0.5))+
  geom_text(aes(label=.group, y = ifelse(period == "baseline", cover*100 + 2,4) 
  ),size=6,position = position_dodge(width = 0.5))+
  labs(x = NULL, y = "% Sponge Cover")+
  scale_x_discrete(
    labels = c(
      "baseline" = "Baseline",
      "dieoff" = "Die-off",
      "after 6mo" = "After 4 Months",
      "after 12mo" = "After 12 Months"
    )
  ) +
  scale_color_viridis_d(option = 'turbo', 
                        labels = c(
                          "baseline" = "Baseline",
                          "dieoff" = "Die-off",
                          "after 6mo" = "After 4 Months",
                          "after 12mo" = "After 12 Months"
                        ))+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12))




# Correlation Plots -------------------------------------------------------


cover2<-pcover |> 
  dplyr::select( period, category, percentcover, site, plot) |>
  pivot_wider(names_from = category, values_from = percentcover)


correlation_results <- cover2 %>%
  mutate(
    Sponge_Coral_corr = cor(sponge, coral, use = "complete.obs"),
    Sponge_Coral_p = cor.test(sponge, coral)$p.value,
    
    Algae_Coral_corr = cor(algae, coral, use = "complete.obs"),
    Algae_Coral_p = cor.test(algae, coral)$p.value,
    
    Algae_Sponge_corr = cor(algae, sponge, use = "complete.obs"),
    Algae_Sponge_p = cor.test(algae, sponge)$p.value,
    
    Algae_Gorgonian_corr = cor(algae, gorgonian, use = "complete.obs"),
    Algae_Gorgonian_p = cor.test(algae, gorgonian)$p.value
  )




#Algae VS Sponge

algae_sponge_results <- correlation_results %>%
  dplyr::select(Algae_Sponge_corr, Algae_Sponge_p)



AS = ggplot(correlation_results, aes(algae,sponge))+
  geom_point(aes(color = period))+
  geom_smooth(method="lm")+
  scale_color_viridis(option = 'turbo', discrete = T, 
                      labels = c(
                        "baseline" = "Baseline",
                        "dieoff" = "Die-off",
                        "after 6mo" = "After 4 Months",
                        "after 12mo" = "After 12 Months"
                      ))+
  #facet_wrap(~location)+
  theme_classic()+
  labs(color = "Periods", y = "Sponge", x= "Aglae")+
  theme(legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 14, face = "bold")
  ) + 
  geom_text(
    aes(
      x = Inf, y = Inf,  # Position in the top right corner of each facet
      label = paste("r =", round(Algae_Sponge_corr, 2), "\np =", format.pval(Algae_Sponge_p, digits = 3, eps = 0.001))
    ),
    hjust = 1.2, vjust = 1.2, size = 6, color = "black",
    data = algae_sponge_results
  )

#Algae VS Coral

algae_coral_results <- correlation_results %>%
  dplyr::select(Algae_Coral_corr, Algae_Coral_p)



AC = ggplot(correlation_results, aes(algae,coral))+
  geom_point(aes(color = period))+
  scale_color_viridis(option = 'turbo', discrete = T, 
                      labels = c(
                        "baseline" = "Baseline",
                        "dieoff" = "Die-off",
                        "after 6mo" = "After 4 Months",
                        "after 12mo" = "After 12 Months"
                      ))+
  geom_smooth(method="lm")+
  #facet_wrap(~location)+
  theme_classic()+
  labs(color = "Periods", y = "Coral", x= "Algae")+
  theme(legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 14, face = "bold")
  )+
  geom_text(
    aes(
      x = Inf, y = Inf,  # Position in the top right corner of each facet
      label = paste("r =", round(Algae_Coral_corr, 2), "\np =", format.pval(Algae_Coral_p, digits = 3, eps = 0.001))
    ),
    hjust = 1.2, vjust = 1.2, size = 6, color = "black",
    data = algae_coral_results
  )


urch2<-urchins_counts |> 
  dplyr::select(plot, period, counts, mean_size) |> 
  mutate(plot = as.numeric(plot))


urch_join <- cover2 %>%
  left_join(urch2, by = c("plot", "period"))




correlation_results2 <- urch_join %>%
  summarize(
    algae_count_corr = cor(algae, counts, use = "complete.obs"),
    algae_count_p = cor.test(algae, counts)$p.value,
    
    algae_size_corr = cor(algae, mean_size, use = "complete.obs"),
    algae_size_p = cor.test(algae, mean_size)$p.value,
  )





AUC = ggplot(urch_join, aes(counts,algae))+
  geom_point(aes(color = period))+
  scale_color_viridis(option = 'turbo', discrete = T, 
                      labels = c(
                        "baseline" = "Baseline",
                        "dieoff" = "Die-off",
                        "after 6mo" = "After 4 Months",
                        "after 12mo" = "After 12 Months"
                      ))+
  geom_smooth(method="lm")+
  #facet_wrap(~location)+
  theme_classic()+
  labs(
    x = "Urchin Counts",  # Update x-axis label
    y = "Algae",        # Optional: Update y-axis label for clarity
    color = "Periods"           # Legend title
  )+
  theme(legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 14, face = "bold")
  ) + 
  geom_text(
    aes(
      x = Inf, y = Inf,  # Position in the top right corner of each facet
      label = paste("r =", round(algae_count_corr, 2), "\np =", format.pval(algae_count_p, digits = 3, eps = 0.001))
    ),
    hjust = 1.2, vjust = 1.2, size = 6, color = "black",
    data = correlation_results2
  )





AUS = ggplot(urch_join, aes(mean_size,algae))+
  geom_point(aes(color = period))+
  scale_color_viridis(option = 'turbo', discrete = T, 
                      labels = c(
                        "baseline" = "Baseline",
                        "dieoff" = "Die-off",
                        "after 6mo" = "After 4 Months",
                        "after 12mo" = "After 12 Months"
                      ))+
  geom_smooth(method="lm")+
  #facet_wrap(~location)+
  theme_classic()+
  labs(color = "Periods")+
  theme(legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  # Legend title size
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 14, face = "bold")
  ) + 
  labs(
    x = "Urchin Cluster Size",  # Update x-axis label
    y = "Algae",        # Optional: Update y-axis label for clarity
    color = "Periods"           # Legend title
  )+
  geom_text(
    aes(
      x = Inf, y = Inf,  # Position in the top right corner of each facet
      label = paste("r =", round(algae_size_corr, 2), "\np =", format.pval(algae_size_p, digits = 3, eps = 0.001))
    ),
    hjust = 1.2, vjust = 1.2, size = 6, color = "black",
    data = correlation_results2)




# Macroalgae expansion ----------------------------------------------------





jpgarea = read_csv("C:/Users/nicor/OneDrive - Florida International University/Documents/Git/RAPID_benCom/Data/jpg.csv") |> 
  dplyr::select(cells:year) |> 
  mutate(
    site = as.integer(site),period = case_when(
      year == 2021 ~ "baseline",
      year == 202206 ~ "dieoff",
      year == 202210 ~ "after 6mo",
      year == 202305 ~ "after 12mo",
      TRUE ~ as.character(year) # Default to the original value if none match
    )
  )

lsmetrics = read_csv("C:/Users/nicor/OneDrive - Florida International University/Documents/Git/RAPID_benCom/Data/lsm.csv")


lambda_df = read_csv("C:/Users/nicor/OneDrive - Florida International University/Documents/Git/RAPID_benCom/Data/lambda_df.csv") |> 
  mutate(period = if_else(timepoint == 1, "baseline", "dieoff")) 

df = lsmetrics |> 
  left_join(jpgarea) |> 
  mutate(ta = cells * 1/(scale^2),
         pland = ca/ta *100, pd = np/ta) |> 
  #filter( period != "after 6mo") |> 
  left_join(lambda_df) |> 
  mutate( period = factor (period, levels = c("baseline", "dieoff", "after 6mo","after 12mo")))

df2 = lsmetrics |> 
  left_join(jpgarea) |> 
  mutate(ta = cells * 1/(scale^2),
         pland = ca/ta *100, pd = np/ta) |>
  pivot_longer(cols = c(ca:np, pland, pd), names_to = "metric", values_to = "value") |>
  nest(.by=metric) |>
  mutate(model = map(data, \(x) aov(value ~ period, data = x)),
         summary = map(model, summary),
         p = map_dbl(summary, \(x) x[[1]]$`Pr(>F)`[1])) |> 
  filter(metric %in% c("pland", "pd", "clumpy", "gyrate_mn")) |> 
  unnest(data) |> 
  mutate(period = case_when(
    period == "baseline" ~ "Baseline",
    period == "dieoff" ~ "Die-off",
    period == "after 6mo" ~ "After 4 Months",
    period == "after 12mo" ~ "After 12 Months",
    TRUE ~ as.character(period)),
    # Explicitly set levels in the desired order
    period = factor(period, levels = c("Baseline", "Die-off", "After 4 Months", "After 12 Months")))




macro = ggplot(df2, aes(x = period, y = value, color = period)) +
  geom_boxplot()+
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(~metric, scales = "free_y",labeller = labeller(metric = c(
    "clumpy" = "Clumpy",
    "gyrate_mn" = "Gyration Mean",
    "pd" = "Patch Density ",
    "pland" = "Percent Cover"
  )))+
  labs(
    x = "Period",       # Capitalize the x-axis title
    y = "Value" ,
    color = "Period",
    caption = "Figure 9: Multiple `seascape` metrics used to asses the expansion of macroalgae across periods. "# Capitalize the y-axis title
  )+
  scale_color_viridis_d(option = "turbo",
                        limits = c("Baseline", "Die-off", "After 4 Months","After 12 Months" ) # Adjust order
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    plot.caption = element_text(size = 12, hjust = .5) # Align caption
  )

macro



# Lambda Modeling ---------------------------------------------------------


df = df |> 
  drop_na(lambda)
df
model = glmmTMB(lambda ~ pland + pd + gyrate_mn + clumpy + (1|site), family = gaussian, data = df, na.action = na.fail)
check_model(model)

dmodel = dredge(model) |> 
  filter(delta < 4)

top = which.min(dmodel$df)

top_model = get.models(dmodel, subset = top)[[1]]

summary(top_model)

