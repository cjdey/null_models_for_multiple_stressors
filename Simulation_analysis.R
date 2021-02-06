require(tidyverse)
require(distr)
require(cowplot)

# load the simulation results
param <- read_csv("param.csv")

################################################################

# reshape to one row per model per parameter set
results <- param

results <- results %>% rename(
  simple.addition = SEsa,
  effect.addition = SEea,
  dominance = SEd,
  concentration.addition = SEca,
  stressor.addition = SESAM
)

results$simulation <- 1:nrow(results)

statdata <- results %>%
  pivot_longer(
    cols = c(
      "dominance",
      "simple.addition",
      "effect.addition",
      "concentration.addition",
      "stressor.addition"
    ),
    names_to = "model",
    values_to = "null_prediction"
  )

statdata$null_prediction <- ifelse(statdata$null_prediction >= 1, 0.999, statdata$null_prediction)

statdata$model <- fct_relevel(statdata$model, "simple.addition")

#############################
# Model divergence and model deviation
###########################

DIF <- function(a, b) {
  z <- a - b
  return(z)
} # model divergence

SD <- function(a, b, c, d, e) {
  z <- sd(c(a, b, c, d, e))
  return(z)
} # model deviation

param <- param %>%
  mutate(
    MINmod = pmap_dbl(list(SEsa, SEea, SEd, SEca, SESAM), min),
    MAXmod = pmap_dbl(list(SEsa, SEea, SEd, SEca, SESAM), max)
  )

param <- param %>%
  mutate(divergence = map2_dbl(MAXmod, MINmod, DIF))

param <- param %>%
  mutate(deviation = pmap_dbl(list(SEsa, SEea, SEd, SEca, SESAM), SD))

median(param$deviation)

median(param$divergence)


# Table S1
param$fA <- fct_relevel(param$fA, "linear")

param$fB <- fct_relevel(param$fB, "linear")

summary(M1 <- lm(divergence * 100 ~ SIA + SIB + fA + fB, data = param))

summary(M2 <- lm(deviation * 100 ~ SIA + SIB + fA + fB, data = param))


########
# Figure2A
#########

labels <- tibble(
  name = c("Model divergence", "Model deviation"),
  x = c(75, 43),
  y = c(.02, .03)
)

p2a <- ggplot() +
  geom_density(
    data = param,
    aes(deviation * 100, y = ..count.. / 250000),
    fill = "grey15",
    color = "black",
    alpha = 0.75
  ) +
  geom_density(
    data = param,
    aes(divergence * 100, y = ..count.. / 250000),
    fill = " grey85",
    color = "black",
    alpha = 0.75
  ) +
  theme_bw() +
  geom_hline(aes(yintercept = 0),
    color = "white"
  ) +
  scale_y_continuous(name = "Proportion of\nparameter sets") +
  scale_x_continuous(
    name = "Difference between null models\n(estimated mortality)",
    breaks = c(25, 50, 75),
    labels = c("25%", "50%", "75%"),
    limits = c(0, 105)
  ) +
  theme() +
  geom_label(
    data = filter(labels, name == "Model divergence"),
    aes(x = x, y = y, label = name),
    size = 3,
    fill = "grey85",
    alpha = 0.75
  ) +
  geom_label(
    data = filter(labels, name == "Model deviation"),
    aes(x = x, y = y, label = name),
    size = 3,
    fill = "grey15",
    alpha = 0.75
  )

p2a


###############
## Figure 2b
##############
plotdata <- results %>%
  pivot_longer(
    cols = c(
      "dominance",
      "effect.addition",
      "concentration.addition",
      "stressor.addition"
    ),
    names_to = "model",
    values_to = "null_prediction"
  )

plotdata <- plotdata %>%
  mutate(simple.addition = ifelse(simple.addition > 1, 1, simple.addition))
plotdata <- plotdata %>%
  mutate(abs_divergence = (null_prediction - simple.addition))

plotdata <- plotdata %>%
  mutate(
    model = fct_recode(
      model,
      "Dominance" = "dominance",
      "Multiplicative" = "effect.addition",
      "Concentration\naddition" = "concentration.addition",
      "Stressor\naddition" = "stressor.addition"
    )
  )

# this takes a while to plot because the dataset is so large
p2b <- ggplot(data = plotdata) +
  geom_point(aes(y = abs_divergence, x = simple.addition),
    size = 0.1,
    alpha = 0.1
  ) +
  geom_abline(
    aes(slope = 0, intercept = 0),
    size = 0.5,
    alpha = 1,
    color = "grey"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    limits = c(-0.05, 1.05),
    labels = c("0%", "50%", "100%"),
    name = "Mortality estimate by Simple Addition"
  ) +
  scale_y_continuous(
    name = "Difference from\nSimple Addition",
    breaks = c(-0.5, 0, 0.5, 1),
    labels = c("-50%", "0%", "+50%", "+100%")
  ) +
  coord_fixed(ratio = 1) +
  facet_wrap(~model, nrow = 1) +
  theme(strip.background = element_rect(fill = "white", color = NA))

plot_grid(
  p2a,
  p2b,
  ncol = 1,
  align = "v",
  axis = "l",
  rel_heights =  c(1, 1.3)
)
ggsave("Figure2.jpg",
  device = "jpg",
  height = 5,
  width = 7
)

###############
# Table 2 - Explained variance in mortality predictions
##############

mod_null <-
  glm(null_prediction ~ 1, family = binomial, data = statdata)

SIA_R2 <- 1 - logLik(update(mod_null, ~ +SIA)) / logLik(mod_null)
SIB_R2 <- 1 - logLik(update(mod_null, ~ +SIB)) / logLik(mod_null)
fA_R2 <- 1 - logLik(update(mod_null, ~ +fA)) / logLik(mod_null)
fB_R2 <- 1 - logLik(update(mod_null, ~ +fB)) / logLik(mod_null)
model_R2 <- 1 - logLik(update(mod_null, ~ +model)) / logLik(mod_null)

tibble(
  variable = c(
    "Stressor intensity A",
    "Stressor intensity B",
    "Stressor-effect reln A",
    "Stressor-effect reln B",
    "null model"
  ),
  explained_var = c(SIA_R2, SIB_R2, fA_R2, fB_R2, model_R2)
)


#############################################
########## Figure 3 - which null model is most precautionary?
################################

plotdata <- results %>%
  pivot_longer(
    cols = c(
      "dominance",
      "simple.addition",
      "effect.addition",
      "concentration.addition",
      "stressor.addition"
    ),
    names_to = "model",
    values_to = "null_prediction"
  )

plotdata$null_prediction <- ifelse(plotdata$null_prediction > 1, 1, plotdata$null_prediction)

plotdata <- plotdata %>%
  ungroup() %>%
  group_by(simulation) %>%
  mutate(difference_from_max = max(null_prediction) - null_prediction)

plotdata <- plotdata %>%
  mutate(
    model =
      fct_recode(
        model,
        "Simple addition" = "simple.addition",
        "Dominance" = "dominance",
        "Multiplicative" = "effect.addition",
        "Concentration addition" = "concentration.addition",
        "Stressor addition" = "stressor.addition"
      )
  )

prec <- plotdata %>%
  filter(difference_from_max == 0) %>%
  ungroup() %>%
  group_by(model) %>%
  summarise(nsim = dplyr::n())

f3a <- ggplot() +
  geom_bar(
    data = prec,
    aes(x = reorder(model, nsim), y = nsim / 250000),
    stat = "identity"
  ) +
  theme_bw() +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  ylab("Proportion of parameter sets")
f3a

plotdata <- plotdata %>%
  ungroup() %>%
  group_by(SIA, SIB, model) %>%
  summarise(precaution = max(difference_from_max)) %>% # max difference for given SIs
  mutate(precaution_cat = cut(
    precaution,
    breaks = c(0, 0.01, 0.05, 0.1, 1),
    labels = c("< 0.01", "< 0.05", "< 0.1", "> 0.1"),
    include.lowest = TRUE
  ))

plotdata <- plotdata %>%
  filter(precaution_cat != "> 0.1") %>%
  mutate(precaution_cat = fct_drop(precaution_cat)) %>%
  mutate(precaution_cat = fct_relevel(precaution_cat, "< 0.01", "< 0.05", "< 0.1"))

f3b <- ggplot(plotdata) +
  geom_tile(aes(x = SIA, y = SIB, fill = precaution_cat)) +
  theme_bw() +
  scale_x_continuous(name = "Intensity of stressor A", breaks = c(0, 0.5, 1)) +
  scale_y_continuous(name = "Intensity of stressor B", breaks = c(0, 0.5, 1)) +
  facet_wrap(~model, nrow = 2) +
  scale_fill_manual(values = c("black", "grey45", "grey85")) +
  labs(fill = "Precaution level") +
  theme(legend.position = c(0.83, 0.2)) +
  theme(strip.background = element_blank())
f3b

plot_grid(
  f3a,
  f3b,
  ncol = 1,
  align = "v",
  axis = "l",
  rel_heights =  c(1, 1.5)
)

 ggsave("Figure3.jpg", device="jpg", height=6, width=7)


#############
# FIGURE S1
##########

x <- seq(0, 1, 0.002)
stress_effect <- tibble(
  x = x,
  "linear" = qbeta(x, 1, 1),
  "diminishing" = qbeta(x, 3, 1),
  "accelerating" = qbeta(x, 1, 3),
  "steep extremes" = qbeta(x, 3, 3),
  "steep middle" = qbeta(x, 1 / 3, 1 / 3)
)

stress_effect <- stress_effect %>%
  pivot_longer(-x, names_to = "shape", values_to = "prob")


ggplot(stress_effect) +
  geom_line(aes(x = x, y = prob * 100, color = shape), size = 1.2) +
  theme_bw() +
  scale_y_continuous(name = "mortality (%)") +
  scale_x_continuous(name = "stress intensity") +
  theme() +
  labs(color = "")

# ggsave("FigureS1.jpg", device="jpg", height=4, width=5)

################
# Figure S2
##############

plotdata <- results %>%
  pivot_longer(
    cols = c(
      "dominance",
      "simple.addition",
      "effect.addition",
      "concentration.addition",
      "stressor.addition"
    ),
    names_to = "model",
    values_to = "null_prediction"
  )

plotdata$null_prediction <- ifelse(plotdata$null_prediction > 1, 1, plotdata$null_prediction)

# bin indepedent stress effects into 0.025 blocks
roundto <- 0.025
plotdata$SEA <- round(plotdata$SEA / roundto) * roundto
plotdata$SEB <- round(plotdata$SEB / roundto) * roundto

plotdata <- plotdata %>%
  ungroup() %>%
  group_by(SEA, SEB) %>%
  summarize(
    Minimum = min(null_prediction),
    Maximum = max(null_prediction)
  ) %>%
  pivot_longer(
    cols = c("Minimum", "Maximum"),
    names_to = "min_max",
    values_to = "value"
  )

plotdata$min_max <- fct_relevel(plotdata$min_max, "Minimum")

fs2 <- ggplot(data = plotdata) +
  geom_raster(aes(x = SEA, y = SEB, fill = value * 100)) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
    breaks = c(25, 50, 75),
    labels = c("25%", "50%", "75%")
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c(0, 50, 100),
    limits = c(-0.05, 1.1),
    name = "Individual effect of stressor B (% mortality)"
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c(0, 50, 100),
    name = "Individual effect of\nstressor A (% mortality)"
  ) +
  facet_wrap(~min_max, ncol = 2) +
  coord_fixed(ratio = 1) +
  labs(fill = "joint\nmortality") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8)
  ) +
  guides(fill = guide_colorbar(barheight = 4))

fs2
# ggsave("FigureS2.jpg", device="jpg", height=4, width=7)
