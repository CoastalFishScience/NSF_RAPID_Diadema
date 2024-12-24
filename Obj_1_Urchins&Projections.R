

# Load Libraries ----------------------------------------------------------


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

urchins = read_csv("./data/urchins_clean.csv") |> 
  drop_na(location)|> 
  mutate(period = factor(period, levels = c("baseline", "dieoff", "after 6mo", "after 12mo")))

urchins_counts = read_csv("./data/urchin_counts.csv") |> 
  drop_na(location) |> 
  mutate(period = factor(period, levels = c("baseline", "dieoff", "after 6mo", "after 12mo")))



# Urchins Analysis --------------------------------------------------------


##Counts analysis####

ggplot(urchins_counts, aes(period, counts, fill = location))+
  #geom_boxplot(outliers = F)+
  geom_point(aes(color = location), position = position_dodge(width = .75))+
  geom_boxplot(outliers = F, alpha = 0.5)

urchins_counts_summary = urchins_counts |> 
  group_by(location, period) |> 
  summarise(mean_c = mean(counts), sd_c = sd(counts))

#Calculate rellative differences from baseline
relative_differences <- urchins_counts_summary |>
  pivot_wider(names_from = period, values_from = c(mean_c, sd_c), 
              names_glue = "{.value}_{period}") |>
  mutate(
    relative_difference_6mo = (`mean_c_after 6mo` - mean_c_baseline) / mean_c_baseline * 100,
    relative_difference_12mo = (`mean_c_after 12mo` - mean_c_baseline) / mean_c_baseline * 100
  ) |>
  dplyr::select(location, relative_difference_6mo, `sd_c_after 6mo`, 
                relative_difference_12mo, `sd_c_after 12mo`)

# View the results
print(relative_differences)


# Relative difference between 6 and 12 months
relative_diff_6_to_12mo <- urchins_counts_summary %>%
  pivot_wider(names_from = period, values_from = c(mean_c, sd_c), 
              names_glue = "{.value}_{period}") %>%
  mutate(
    relative_difference_6_to_12mo = (`mean_c_after 12mo` - `mean_c_after 6mo`) / `mean_c_after 6mo` * 100
  ) %>%
  dplyr::select(location, relative_difference_6_to_12mo, `sd_c_after 6mo`, `sd_c_after 12mo`)

# View the results
print(relative_diff_6_to_12mo)




###
#Models
###


###Negative binomial####

#Final model for counts
counts_periods = glmmTMB(counts ~ period + (1|site), 
                         family = nbinom2(link = "log"), data = urchins_counts)
car::Anova(counts_periods)
summary(counts_periods)
#capture.output(summary(counts_periods), file="./tables/summary.counts.period.txt")

gg.counts_period = ggpredict(counts_periods, terms = c("period")) |> 
  rename(period = x, counts = predicted)

summary(gg.counts_period)

#check region comparisons
model_means_nbcounts = emmeans(object = counts_periods,
                               specs = "period")

# pairs(model_means_nbcounts, adjust = 'bonf')
pairs(model_means_nbcounts)

contrast(model_means_nbcounts)

#capture.output(pairs(model_means_nbcounts))
# add letters to each mean
model_means_cld_nbcounts = cld(object = model_means_nbcounts,
                               Letters = letters,
                               alpha = 0.05)



#plot of fit for location
num_clusters = ggplot(gg.counts_period, aes(period, counts, colour = period))+
  geom_pointrange(aes(ymax = conf.high, ymin = conf.low), size = 1, linewidth = 2)+
  labs(x = "Periods relative to urchin dieoff", y = "Counts of urchin clusters")+
  scale_color_viridis(option = 'turbo', discrete = TRUE)+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
num_clusters



## Cluster size analysis ####

ggplot(urchins_counts, aes(period, mean_size, fill = location))+
  #geom_boxplot(outliers = F)+
  geom_point(aes(color = location), position = position_dodge(width = .75))+
  geom_boxplot(outliers = F, alpha = 0.5)


###
#Models
###


###Gaussian####

#For the models, we considered only period since we only testing for cluster size
#when clusters are present and because there are combination with location and period with
#with zero clusters

cluster_size = urchins_counts |> 
  dplyr::select(period, site, mean_size) |> 
  filter(mean_size > 0)

#Final model for cluster size
size_gaussian = glmmTMB(mean_size ~ period + (1|site), 
                        family = gaussian(link = "log"), data = cluster_size)

car::Anova(size_gaussian)
#capture.output(summary(size_gaussian), file="./tables/summary.size.period.txt")

gg.size_period = ggpredict(size_gaussian, terms = c("period")) |> 
  rename(period = x, size = predicted)

#check region comparisons
model_means_size = emmeans(object = size_gaussian,
                           specs = "period")

# pairs(model_means_nbcounts, adjust = 'bonf')
pairs(model_means_size)

contrast(model_means_size)

#capture.output(pairs(model_means_size))
# add letters to each mean
model_means_cld_size = cld(object = model_means_size,
                           Letters = letters,
                           alpha = 0.05)



#plot of fit for location
size_clusters = ggplot(gg.size_period, aes(period, size, colour = period))+
  geom_pointrange(aes(ymax = conf.high, ymin = conf.low), size = 1, linewidth = 2)+
  labs(x = "Periods relative to urchin dieoff", y = "Size of urchin clusters")+
  scale_color_viridis(option = 'turbo', discrete = TRUE)+
  # geom_hline(yintercept = 0, linetype = 2, color = "red") +
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
size_clusters




ggarrange(num_clusters, size_clusters,
          labels = c('a)', 'b)'),
          nrow = 1, legend = 'none')




# Load Data ---------------------------------------------------------------

t1demographydata = read_csv("Data/t1demographyfull.csv") 

t2demographydata = read_csv("Data/t2demographyfull.csv")



# Porites astreoides ------------------------------------------------------

#Matrix 2 2022=20234

##Matrix 2 2022-2023 ####

t2p = t2demographydata |>  
  mutate(class=case_when(class == "empty" ~ "porites astreoides",
                         TRUE ~ class)) |> 
  filter(class == "porites astreoides")

t2precruits = sum(t2p$action == "born")

#70 is the size for fertile adults [citation needed]
adults_year1t2p <- sum(t2p$area1 > 70, na.rm = TRUE) * 0.61


#number of recruits divided by total number of FERTILE adults

#t2p_fecundity_rate = t2precruits / (adults_year1t2pdull)
t2p_fecundity_rate = t2precruits / (adults_year1t2p)


### Create matrix model ####

t2pdf = t2p |> 
  filter(action!="born")

t2pdf = t2pdf |>
  mutate(juvenile = ifelse(area1>70,t2p_fecundity_rate,0),
         coral=row_number()) |> 
  dplyr::select(coral,stage,fate,juvenile)

stages= c("juvenile","intermediate","adult")

t2pdf$stage = ordered(t2pdf$stage, stages)


#popbio requires dataframe not tibble

t2pdf = as.data.frame(t2pdf)

#Creates population matrix

t2pmatrix = projection.matrix(t2pdf, sort=stages,TF=F)


###Population projections####

#initial population

#juveniles, intermediate, adults

n2p=c(6400,901,444)

p2=pop.projection(t2pmatrix,n2p,15)



t2p_stage_vectors = p2$stage.vectors



# Convert the stage vectors to a data frame

t2p_df_stage = as.data.frame(t(t2p_stage_vectors))

t2p_df_stage$year = rownames(t2p_df_stage)

rownames(t2p_df_stage) = NULL


# Reshape the data frame into long format for ggplot

t2p_df_stage_long = gather(t2p_df_stage, class, value,-year)

t2p_df_stage_long$year=as.numeric(t2p_df_stage_long$year)


t2p_sum_values = t2p_df_stage_long |>
  group_by(year) |>
  summarize(sum_value = sum(value))

##Matrix 1 2021-2022 ####


###Calculate  Fecundity####

t1p = t1demographydata |> 
  mutate(class=case_when(class == "empty" ~ "porites astreoides",
                         TRUE ~ class)) |>  
  filter(class == "porites astreoides")



t1precruits = sum(t1p$action == "born") 

#70 is the size for fertile adults [citation needed]

adults_year1t1p <- sum(t1p$area1 > 70, na.rm = TRUE) * 0.61


#number of recruits divided by total number of FERTILE adults 
t1p_fecundity_rate = t1precruits / (adults_year1t1p)


###Create matrix model ####

t1pdf = t1p |> 
  filter(action!="born")

t1pdf = t1pdf |> 
  mutate(juvenile = ifelse(area1>70,t1p_fecundity_rate,0),
         coral=row_number()) |> 
  dplyr::select(coral,stage,fate,juvenile)

stages= c("juvenile","intermediate","adult")

t1pdf$stage = ordered(t1pdf$stage, stages)

#popbio requires dataframe not tibble
t1pdf = as.data.frame(t1pdf)

#Creates population matrix
t1pmatrix = projection.matrix(t1pdf, sort=stages, TF = F)

###Population projections####

#initial population

#juveniles, intermediate, adults


n1p=c(5075,833,487)

p1=pop.projection(t1pmatrix,n1p,15)

t1p_stage_vectors = p1$stage.vectors

t1p_df_stage = as.data.frame(t(t1p_stage_vectors))
t1p_df_stage$year = rownames(t1p_df_stage)
rownames(t1p_df_stage) = NULL

t1p_df_stage_long = gather(t1p_df_stage, class, value,-year)
t1p_df_stage_long$year=as.numeric(t1p_df_stage_long$year)

t1p_sum_values = t1p_df_stage_long |>
  group_by(year) |>
  summarize(sum_value = sum(value))

t1p_df_stage_long <- t1p_df_stage_long %>%
  mutate(class = case_when(
    class == "juvenile" ~ "Juvenile",
    class == "intermediate" ~ "Intermediate",
    class == "adult" ~ "Adult",
    TRUE ~ class  # Keep other values unchanged if any
  ))

t1plot = ggplot(t1p_df_stage_long, aes(x = year, y = value, color = class)) +
  geom_line(linewidth = 1.8) +
  geom_line(data = t1p_sum_values, aes(x = year, y = sum_value, color = "Total"), linetype = "dotted",linewidth = 1.5) +
  labs(x = "Time", y = "Population Size", color = "Class", title = "2021-2022 based projection") +
  theme_classic() +
  geom_text(
    x = max(t2p_df_stage_long$year),   # Set the x-coordinate (e.g., last year)
    y = max(t2p_df_stage_long$value), # Set the y-coordinate (e.g., max value)
    label = expression(lambda == 1.1),  # Use bquote for lambda symbol and value
    hjust = 1.1, vjust = -23.7, size = 5, fontface = "bold", color = "black"
  ) +
  #ylim(0, 1000)+
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.caption = element_text(size = 12, hjust = .5,
                                    vjust = -4, margin = margin(t = 10)),
        plot.margin = margin(t = 40, b=40),
        legend.position = "none"
        
  )








# Add the sum_values as another layer to your ggplot

t2plot = ggplot(t2p_df_stage_long, aes(x = year, y = value, color = class)) +
  geom_line(linewidth = 1.8) +
  geom_line(data = t2p_sum_values, aes(x = year, y = sum_value, color = "Total"), linetype = "dotted",linewidth = 1.5) +
  labs(x = "Time", y = "Population Size", color = "Class", title = "2022-2023 based projection") +
  theme_classic() +
  geom_text(
    x = max(t2p_df_stage_long$year),   # Set the x-coordinate (e.g., last year)
    y = max(t2p_df_stage_long$value), # Set the y-coordinate (e.g., max value)
    label = expression(lambda == 0.988),  # Use bquote for lambda symbol and value
    hjust = 1.1, vjust = -5.5, size = 5, fontface = "bold", color = "black"
  ) +
  #ylim(0, 1000)+
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.caption = element_text(size = 12, hjust = .5,
                                    vjust = -4, margin = margin(t = 10)),
        plot.margin = margin(t = 40, b=40),
        legend.position = "none"
        
  )



#Stochastic model

matrices2 <- list(t1pmatrix, t2pmatrix)

n <- c(100,50,30)

stochastic_projection <- stoch.projection(matrices2, n,prob = c(.1,.9),tmax=1000, nreps=5000, verbose=T)
stochGR = stoch.growth.rate(matrices2,prob = c(.1,.9), verbose=T)

stoch.projection.trajectories <- function(matrices, n0, tmax = 50, nreps = 5000, prob = NULL, 
                                          nmax = NULL, sumweight = rep(1, length(n0)), verbose = FALSE) {
  if (!is.list(matrices)) {
    stop("Please use a list of matrices as input")
  }
  # Initialize a 3D array to store trajectories: [time, stages, replicates]
  trajectories <- array(0, dim = c(tmax, length(n0), nreps))
  for (i in 1:nreps) {
    A <- sample(matrices, tmax, replace = TRUE, prob = prob)
    n <- n0
    if (verbose) {
      if (i == 1 || i%%100 == 0) {
        message("Starting projections for nrep ", i)
      }
    }
    for (j in 1:tmax) {
      n <- A[[j]] %*% n
      if (!is.null(nmax)) {
        if (sum(n * sumweight) > nmax) 
          n <- n * (nmax/sum(n * sumweight))
      }
      trajectories[j, , i] <- n  # Store population size at each time step
    }
  }
  trajectories
}

#Use new function
trajectories <- stoch.projection.trajectories(matrices2, n, prob = c(.9, .1), tmax = 100, nreps = 5000, verbose = TRUE)

# Mean across replicates
average_trajectory <- apply(trajectories, c(1, 2), mean)  
lower_bound <- apply(trajectories, c(1, 2), quantile, probs = 0.05)  # 5th percentile for lower bound
upper_bound <- apply(trajectories, c(1, 2), quantile, probs = 0.95) 
stage_names <- c("Juvenile", "Intermediate", "Adult")


df <- as.data.frame(average_trajectory)
colnames(df) <- stage_names
df <- df %>%
  mutate(Year = 1:nrow(df)) %>%
  pivot_longer(cols = -Year, names_to = "Stage", values_to = "Population")

lower_df <- as.data.frame(lower_bound)
colnames(lower_df) <- stage_names
lower_df <- lower_df %>%
  mutate(Year = 1:nrow(lower_df)) %>%
  pivot_longer(cols = -Year, names_to = "Stage", values_to = "Lower")


upper_df <- as.data.frame(upper_bound)
colnames(upper_df) <- stage_names
upper_df <- upper_df %>%
  mutate(Year = 1:nrow(upper_df)) %>%
  pivot_longer(cols = -Year, names_to = "Stage", values_to = "Upper")

df <- df %>%
  left_join(lower_df, by = c("Year", "Stage")) %>%
  left_join(upper_df, by = c("Year", "Stage"))


stoch1 <- ggplot(df, aes(x = Year, y = Population, color = Stage)) +
  geom_line(size = 1.2) +  # Average trajectories
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Stage), alpha = 0.2, color = NA) +  # Confidence intervals
  labs(
    title = "Average Population Dynamics Over Time (10% Disturbance)",
    x = "Year",
    y = "Population Size",
    color = "Stage"
  ) +
  theme_classic() +
  geom_text(
    x = max(df$Year),  # Set the x-coordinate (e.g., last year)
    y = max(df$Population),  # Set the y-coordinate (e.g., max value)
    label = expression(lambda < 1),  # Use bquote for lambda symbol and value
    hjust = -1, vjust = -1, size = 5, fontface = "bold", color = "black"
  ) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.caption = element_text(size = 12, hjust = .5,
                                    vjust = -4, margin = margin(t = 10)),
        plot.margin = margin(t = 40, b=40),
        legend.position = "none"
        
  )






#Stochastic model2

trajectories2 <- stoch.projection.trajectories(matrices2, n, prob = c(.1, .9), tmax = 100, nreps = 5000, verbose = TRUE)



# Mean across replicates
average_trajectory2 <- apply(trajectories2, c(1, 2), mean)  
lower_bound2 <- apply(trajectories2, c(1, 2), quantile, probs = 0.05)  # 5th percentile for lower bound
upper_bound2 <- apply(trajectories2, c(1, 2), quantile, probs = 0.95) 
stage_names <- c("Juvenile", "Intermediate", "Adult")


df2 <- as.data.frame(average_trajectory2)
colnames(df2) <- stage_names
df2 <- df2 %>%
  mutate(Year = 1:nrow(df2)) %>%
  pivot_longer(cols = -Year, names_to = "Stage", values_to = "Population")

lower_df2 <- as.data.frame(lower_bound2)
colnames(lower_df2) <- stage_names
lower_df2 <- lower_df2 %>%
  mutate(Year = 1:nrow(lower_df2)) %>%
  pivot_longer(cols = -Year, names_to = "Stage", values_to = "Lower")


upper_df2 <- as.data.frame(upper_bound2)
colnames(upper_df2) <- stage_names
upper_df2 <- upper_df2 %>%
  mutate(Year = 1:nrow(upper_df2)) %>%
  pivot_longer(cols = -Year, names_to = "Stage", values_to = "Upper")

df2 <- df2 %>%
  left_join(lower_df2, by = c("Year", "Stage")) %>%
  left_join(upper_df2, by = c("Year", "Stage"))


stoch2 <- ggplot(df2, aes(x = Year, y = Population, color = Stage)) +
  geom_line(size = 1.2) +  # Average trajectories
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Stage), alpha = 0.2, color = NA) +  # Confidence intervals
  labs(
    title = "Average Population Dynamics Over Time (90% Disturbance)",
    x = "Year",
    y = "Population Size",
    color = "Stage"
  ) +
  theme_classic() +
  geom_text(
    x = max(df2$Year),  # Set the x-coordinate (e.g., last year)
    y = max(df2$Population),  # Set the y-coordinate (e.g., max value)
    label = expression(lambda < 1),  # Use bquote for lambda symbol and value
    hjust = -1, vjust = -1, size = 5, fontface = "bold", color = "black"
  ) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.caption = element_text(size = 12, hjust = .5,
                                    vjust = -4, margin = margin(t = 10)),
        plot.margin = margin(t = 40, b=40),
        legend.position = "none"
        
  )


###Sensitivity and Elasticity Analysis -------------------------------------

sens<-stoch.sens(matrices2, tlimit=100)

elasticity_matrix <- sens$elasticities
colnames(elasticity_matrix) <- c("Juvenile", "Intermediate", "Adult")
rownames(elasticity_matrix) <- c("Juvenile", "Intermediate", "Adult")

elas_df <- as.data.frame(elasticity_matrix)

elas_df <- elas_df %>%
  rename(To = Var1, From = Var2, Elasticity = Freq)


#make a heat map
elasticity = ggplot(elas_df, aes(x = From, y = To, fill = Elasticity)) +
  geom_tile() +
  #geom_text(aes(label = round(Elasticity, 3)), color = "white", fontface = "bold") +
  scale_fill_viridis_c() +
  labs(
    title = "Elasticity Matrix",
    x = "From (Contributing Stage)",
    y = "To (Receiving Stage)",
    fill = "Elasticity"
  ) +
  theme_classic()+
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Prepare data for sensitivities
sens_df <- as.data.frame(stoch.sens(matrices2, tlimit=100)$sensitivities)
colnames(sens_df) <- c("Juvenile", "Intermediate", "Adult")
sens_df$Stage <- c("Juvenile", "Intermediate", "Adult")
sens_long <- pivot_longer(sens_df, cols = -Stage, names_to = "To", values_to = "Sensitivity")
sens_long$Stage <- factor(sens_long$Stage, levels = rev(c("Adult", "Intermediate", "Juvenile")))
sens_long$To <- factor(sens_long$To, levels = c("Juvenile", "Intermediate", "Adult"))

# Sensitivity heatmap
sensitivity = ggplot(sens_long, aes(x = Stage, y = To , fill = Sensitivity)) +
  geom_tile() +
  labs(title = "Sensitivity Matrix", x = "From (Contributing Stage)", y = "To (Receiving Stage)", fill = "Sensitivity") +
  scale_fill_viridis_c() +
  theme_classic()+
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


sens_plots <- ggarrange(
  sensitivity, elasticity,
  ncol = 2,
  #align = "v",
  common.legend = F, # To have a common legend for all plots
  legend = "right",
  labels = "AUTO" # Place the legend on the right side if desired
)

sens_plots <- annotate_figure(
  sens_plots, 
  bottom = text_grob(
    "Figure 7. For *Porites astreoides* A) Sensitivity analysis of stochasitc model with equal probabilities of disturbance. B) Elasticity analysis of stochasitc model with equal probabilities of disturbance.",
    size = 14 # Adjust vertical positioning to reduce space
  )
)
sens_plots



# DLAB --------------------------------------------------------------------



## Matrix 1 2021-2022----------------------------------------------------------

###Calculate Fecundity####

t1d = t1demographydata |> 
  filter(class == "dlab")

t1drecruits = sum(t1d$action == "born")

#70 is for porites, find DLAB reproductive size

adults_year1t1d = sum(t1d$area1 > 70, na.rm = TRUE)

#number of recruits divided by total number of FERTILE adults

t1d_fecundity_rate = t1drecruits / (adults_year1t1d)


### Create matrix model ####

t1ddf = t1d |> 
  filter(action!="born")

t1ddf = t1ddf |>
  mutate(juvenile = ifelse(area1>70,t1d_fecundity_rate,0),
         coral=row_number()) |> 
  dplyr::select(coral,stage,fate,juvenile)

stages= c("juvenile","intermediate","adult")

t1ddf$stage = ordered(t1ddf$stage, stages)


#popbio requires dataframe not tibble

t1ddf = as.data.frame(t1ddf)

#Creates population matrix

t1dmatrix = projection.matrix(t1ddf, sort=stages)




##Matrix 2 2022 - 2023 ----------------------------------------------------------


###Calculate Fecundity####

t2d = t2demographydata |> 
  filter(class == "dlab")

t2drecruits = sum(t2d$action == "born")

#70 is the size for fertile adults [citation needed]

adults_year1t2d = sum(t2d$area1 > 70, na.rm = TRUE)

#number of recruits divided by total number of FERTILE adults

t2d_fecundity_rate = t2drecruits / (adults_year1t2d)



### Create matrix model ####

t2ddf = t2d |> 
  filter(action!="born")

t2ddf = t2ddf |>
  mutate(juvenile = ifelse(area1>70,t2d_fecundity_rate,0),
         coral=row_number()) |> 
  dplyr::select(coral,stage,fate,juvenile)


stages= c("juvenile","intermediate","adult")

t2ddf$stage = ordered(t2ddf$stage, stages)


#popbio requires dataframe not tibble

t2ddf = as.data.frame(t2ddf)

#Creates population matrix

t2dmatrix = projection.matrix(t2ddf, sort=stages)

#Stochastic model2

matrices3<- list(t1dmatrix, t2dmatrix)

trajectories3 <- stoch.projection.trajectories(matrices3, n, prob = c(.9, .1), tmax = 100, nreps = 5000, verbose = TRUE)

# Mean across replicates
average_trajectory3 <- apply(trajectories3, c(1, 2), mean)  


stage_names3 <- c("Juvenile", "Intermediate", "Adult")

# Create a data frame for plotting
df3 <- as.data.frame(average_trajectory3)
colnames(df3) <- stage_names3  # Name the stages appropriately
df3 <- df3 %>%
  mutate(Year = 1:nrow(df3)) %>%  # Add time steps
  pivot_longer(cols = -Year, names_to = "Stage", values_to = "Population")  # Convert to long format

# Plot the data
stoch3 = ggplot(df3, aes(x = Year, y = Population, color = Stage)) +
  geom_line(size = 1.2) +
  labs(
    title = "Average Population Dynamics Over Time (10% Disturbance)",
    x = "Year",
    y = "Population Size",
    color = "Stage"
  ) +
  theme_classic()+
  geom_text(
    x = max(df3$year),   # Set the x-coordinate (e.g., last year)
    y = max(df3$value), # Set the y-coordinate (e.g., max value)
    label = expression(lambda < 0.90),  # Use bquote for lambda symbol and value
    hjust = -7.0, vjust = -31.5, size = 5, fontface = "bold", color = "black"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


###Quasi Extinction####
x <- stoch.quasi.ext(matrices3, n, Nx=15, tmax=30, nreps=100,prob = c(.9, .1))

x_df <- as.data.frame(x)  # Convert matrix to data frame
x_df$Years <- 1:nrow(x_df)  # Add a "Years" column
x_long <- pivot_longer(x_df, cols = -Years, names_to = "Scenario", values_to = "Probability")
x_long$Scenario <- as.numeric(as.factor(x_long$Scenario))

ext = ggplot(x_long, aes(x = Years, y = Probability, color = as.factor(Scenario))) +
  geom_line(size = 1.2) +
  labs(
    title = "Quasi-Extinction Probability Over Time",
    x = "Years",
    y = "Quasi-Extinction Probability",
    color = "Iteration"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


sens2<-stoch.sens(matrices3, tlimit=100)

elasticity_matrix2 <- sens2$elasticities
colnames(elasticity_matrix2) <- c("Juvenile", "Intermediate", "Adult")
rownames(elasticity_matrix2) <- c("Juvenile", "Intermediate", "Adult")

elas_df2 <- as.data.frame(elasticity_matrix2)

elas_df2 <- elas_df2 %>%
  rename(To = Var1, From = Var2, Elasticity = Freq)


#make a heat map
elasticity2 = ggplot(elas_df2, aes(x =From , y =To , fill = Elasticity)) +
  geom_tile() +
  #geom_text(aes(label = round(Elasticity, 3)), color = "white", fontface = "bold") +
  scale_fill_viridis_c() +
  labs(
    title = "Elasticity Matrix",
    x = "From (Contributing Stage)",
    y = "To (Receiving Stage)",
    fill = "Elasticity"
  ) +
  theme_classic()+
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Prepare data for sensitivities
sens_df2 <- as.data.frame(stoch.sens(matrices3, tlimit=100)$sensitivities)
colnames(sens_df2) <- c("Juvenile", "Intermediate", "Adult")
sens_df2$Stage <- c("Juvenile", "Intermediate", "Adult")
sens_long2 <- pivot_longer(sens_df2, cols = -Stage, names_to = "To", values_to = "Sensitivity")
sens_long2$Stage <- factor(sens_long2$Stage, levels = rev(c("Adult", "Intermediate", "Juvenile")))
sens_long2$To <- factor(sens_long2$To, levels = c("Juvenile", "Intermediate", "Adult"))

# Sensitivity heatmap
sensitivity2 = ggplot(sens_long2, aes(x =Stage , y =To , fill = Sensitivity)) +
  geom_tile() +
  labs(title = "Sensitivity Matrix", x = "From (Contributing Stage)", y = "To (Receiving Stage)", fill = "Sensitivity") +
  scale_fill_viridis_c() +
  theme_classic()+
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )



dlabplots <- ggarrange(
  stoch3,ext,sensitivity2, elasticity2,
  ncol = 2, nrow=2,
  #align = "v",
  common.legend = F, # To have a common legend for all plots
  legend = "right",
  labels = "AUTO" # Place the legend on the right side if desired
)

dlabplots <- annotate_figure(
  dlabplots, 
  bottom = text_grob(
    "Figure 8. A) Population dynamics of DLAB from stocastic model with a 10% chance of a macroalgal bloom occuring, showing a rapid decline in population size. \n B) Multiple model predictors of time until population collapse (less than 20 individuals). C) Sensitivity analysis of stochastic model with equal chance of \ndisturbance.  D) Elasticity analysis of stochastic model with equal chance of disturbance",
    size = 14 # Adjust vertical positioning to reduce space
  )
)

dlabplots


# Finito ------------------------------------------------------------------




