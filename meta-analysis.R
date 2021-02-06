require(tidyverse)
require(distr)

##############
# Re-analysis of data from Crain et al. and Darling & Cote meta-analyses
############

meta_crain <- read_csv("Crain et al. survival data.csv")
# note, this data file was provided by Ben Halpern and can be requested by contacting
# them at halpern@nceas.ucsb.edu

meta_crain <- meta_crain %>%
  mutate(
    baseline_mort = 100 - baseline_survival,
    A_mort = 100 - Stressor_effect_A,
    B_mort = 100 - Stressor_effect_B,
    AB_mort = 100 - Stressor_effect_AB
  )

meta_crain <- meta_crain %>%
  mutate(
    SEA = A_mort - baseline_mort,
    SEB = B_mort - baseline_mort,
    SEAB = AB_mort - baseline_mort
  )

meta_crain <- meta_crain %>% select(species, Stressor_A, Stressor_B, baseline_mort, A_mort, B_mort, AB_mort, SEA, SEB, SEAB)


###
meta_darling = read_csv("Darling_data.csv") %>% select(-Study)
# note, this data file was provided by Emily Darling and can be requested by contacting
# them at edarling@wcs.org

meta_darling$Stressor = rep(c("A", "B", "AB"), nrow(meta_darling)/3)

meta_darling = meta_darling %>% pivot_wider(id_cols = c(Expt, Species), 
                             names_from = Stressor, 
                             values_from = c( Stress, Control_mortality, Stress_mortality ) )

meta_darling = meta_darling %>% mutate(baseline_mort = Control_mortality_A)

meta_darling$baseline_mort[which(meta_darling$baseline_mort == "<10%")]<- 10

meta_darling = meta_darling %>% mutate(baseline_mort = as.numeric(baseline_mort))

meta_darling = meta_darling %>% 
  select(-Control_mortality_A, -Control_mortality_B, -Control_mortality_AB, -Stress_AB)

meta_darling <- meta_darling %>%
  mutate(
    SEA = Stress_mortality_A - baseline_mort,
    SEB = Stress_mortality_B - baseline_mort,
    SEAB = Stress_mortality_AB - baseline_mort
  )

meta_darling = meta_darling %>% select(Species, Stress_A, Stress_B, baseline_mort, Stress_mortality_A, Stress_mortality_B, Stress_mortality_AB, SEA, SEB, SEAB)

names(meta_darling) <- names(meta_crain)

meta_data = rbind(meta_crain, meta_darling)

####

meta_data <- meta_data %>%
  filter(SEA >= 0, SEB >= 0)  # remove experiments where 'stressors' increase survival

# uniquely label for plotting
LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))

meta_data$dataset <- LETTERS702[1:nrow(meta_data)]


# define the null models

simple_addition <- function(SEA, SEB, base) {
  pred <- (SEA + SEB + base)
  if (pred > 100) {
    pred <- 100
  }
  return(pred)
}

dominance <- function(SEA, SEB, base) {
  pred <- max(SEA, SEB) + base
  if (pred > 100) {
    pred <- 100
  }
  return(pred)
}

multiplicative <- function(SEA, SEB, base) {
  SEA_proportion <- SEA / 100
  SEB_proportion <- SEB / 100
  pred <- (SEA_proportion + SEB_proportion - (SEA_proportion * SEB_proportion)) * 100 + base
  if (pred > 100) {
    pred <- 100
  }
  return(pred)
}

stressor_addition <- function(SEA, SEB, base) {
  x1 <- qbeta(SEA / 100, 3.2, 3.2)
  x2 <- qbeta(SEB / 100, 3.2, 3.2)
  pred <- pbeta((x1 + x2), 3.2, 3.2) * 100 + base
  if (pred > 100) {
    pred <- 100
  }
  return(pred)
}

# parameters for the stressor-effect relationships of the concentration addition model
se_params <- expand_grid(
  fA = c("linear", "accelerating", "diminishing", "steep_mid", "steep_ext"), # stress-effect shapes
  fB = c("linear", "accelerating", "diminishing", "steep_mid", "steep_ext")
) # stress-effect shapes

se_params$fAa <- se_params$fAb <- se_params$fBa <- se_params$fBb <- NA

se_params$fAa <- ifelse(se_params$fA == "linear", 1,
  ifelse(se_params$fA == "diminishing", 3,
    ifelse(se_params$fA == "accelerating", 1,
      ifelse(se_params$fA == "steep_ext", 3,
        ifelse(se_params$fA == "steep_mid", (1 / 3), NA)
      )
    )
  )
)

se_params$fBa <- ifelse(se_params$fB == "linear", 1,
  ifelse(se_params$fB == "diminishing", 3,
    ifelse(se_params$fB == "accelerating", 1,
      ifelse(se_params$fB == "steep_ext", 3,
        ifelse(se_params$fB == "steep_mid", (1 / 3), NA)
      )
    )
  )
)

se_params$fAb <- ifelse(se_params$fA == "linear", 1,
  ifelse(se_params$fA == "diminishing", 1,
    ifelse(se_params$fA == "accelerating", 3,
      ifelse(se_params$fA == "steep_ext", 3,
        ifelse(se_params$fA == "steep_mid", (1 / 3), NA)
      )
    )
  )
)

se_params$fBb <- ifelse(se_params$fB == "linear", 1,
  ifelse(se_params$fB == "diminishing", 1,
    ifelse(se_params$fB == "accelerating", 3,
      ifelse(se_params$fB == "steep_ext", 3,
        ifelse(se_params$fB == "steep_mid", (1 / 3), NA)
      )
    )
  )
)

concentration_addition <- function(SEA, SEB, base, out) {
  if (!out %in% c("min", "max", "mean", "list")) {
    stop()
  }

  scope <- (100 - base)

  preds <- NA

  # calculate the CA predictions for each combination of stressor-effect relationships
  for (i in 1:nrow(se_params)) {
    B1 <- Beta(shape1 = se_params$fAa[i], shape2 = se_params$fAb[i])
    Z1 <- B1 * scope

    B2 <- Beta(shape1 = se_params$fBa[i], shape2 = se_params$fBb[i])
    Z2 <- B2 * scope

    SIA <- p(Z1)(SEA)
    SIB <- p(Z2)(SEB)

    EC1 <- p(Z1)(scope / 2)
    EC2 <- p(Z2)(scope / 2)

    lambda <- as.numeric(EC1 / EC2)

    corrected_SI <- SIA + (lambda * SIB)
    if (corrected_SI > 1) {
      corrected_SI <- 1
    }

    preds[i] <- (q(Z1)(corrected_SI) + base)
  }

  if (out == "min") {
    return(min(preds))
  }
  if (out == "max") {
    return(max(preds))
  }
  if (out == "mean") {
    return(mean(preds))
  }
  if (out == "list") {
    return(preds)
  }
}


# calculate the null model predictions for each data experiment

meta_data <- meta_data %>% mutate(
  simple.addition = pmap_dbl(list(SEA = SEA, SEB = SEB, base = baseline_mort), simple_addition),
  dominance = pmap_dbl(list(SEA = SEA, SEB = SEB, base = baseline_mort), dominance),
  multiplicative = pmap_dbl(list(SEA = SEA, SEB = SEB, base = baseline_mort), multiplicative),
  stressor.addition = pmap_dbl(list(SEA = SEA, SEB = SEB, base = baseline_mort), stressor_addition),
  concentration.addition.min = pmap_dbl(list(SEA = SEA, SEB = SEB, base = baseline_mort, out = "min"), concentration_addition),
  concentration.addition.max = pmap_dbl(list(SEA = SEA, SEB = SEB, base = baseline_mort, out = "max"), concentration_addition),
  concentration.addition.mean = pmap_dbl(list(SEA = SEA, SEB = SEB, base = baseline_mort, out = "mean"), concentration_addition),
  concentration.addition.list = pmap(list(SEA = SEA, SEB = SEB, base = baseline_mort, out = "list"), concentration_addition)
)


# calculate the minimum and maximum null model predictions (i.e. the interval)
meta_data <- meta_data %>%
  mutate(
    null.model.min = pmap_dbl(list(simple.addition, dominance, multiplicative, stressor.addition, concentration.addition.min), min),
    null.model.max = pmap_dbl(list(simple.addition, dominance, multiplicative, stressor.addition, concentration.addition.max), max)
  )

meta_data = meta_data %>% mutate_at(c("AB_mort", "null.model.min", "null.model.max"),
                                      round, 3)

meta_data = meta_data %>% mutate(outcome = ifelse(AB_mort < null.model.min, "ant",
                                                  ifelse(AB_mort > null.model.max, "syn", "as_pred")))

table(meta_data$outcome)/nrow(meta_data)

meta_data %>% 
  mutate(con_add = ifelse(abs(AB_mort - simple.addition) <= 5, 1,0 )) %>%
  group_by(con_add) %>% summarize(count = dplyr::n())

# re_order the dataset for plotting
meta_data <- meta_data %>% mutate(null.mid = (null.model.max + null.model.min)/2) 

meta_data$dataset <- fct_reorder(meta_data$dataset, meta_data$null.mid)
meta_data$dataset <- fct_reorder(meta_data$dataset, meta_data$outcome)

meta_data = meta_data %>% arrange(dataset)

meta_data%>% select(dataset,  null.mid, outcome)

#Figure 4
roundto <- 1

plotdata = meta_data %>% mutate(null.model.min = round(null.model.min/roundto)*roundto,
                                null.model.max = round(null.model.max/roundto)*roundto,
                                AB_mort = round(AB_mort/roundto)*roundto)

waffle = expand_grid(mortality = (seq(0, 100, roundto)),
                     dataset = plotdata$dataset)

waffle$fill = "white"
waffle$cat = NA
for(i in 1:nrow(waffle)){
  
  m = waffle$mortality[i]
  d = waffle$dataset[i]
  
  nmin = plotdata$null.model.min[which(plotdata$dataset == d)]
  nmax = plotdata$null.model.max[which(plotdata$dataset == d)]
  obs = plotdata$AB_mort[which(plotdata$dataset == d)]
  outcome = plotdata$outcome[which(plotdata$dataset == d)]

  if(m >= nmin & m <= nmax & m != obs){waffle$fill[i]<- "grey"} else
  if( m == obs & outcome == "as_pred"){waffle$fill[i]<- "black"} else
    if( m == obs & outcome == "ant"){waffle$fill[i]<- "royalblue"} else  
      if( m == obs & outcome == "syn"){waffle$fill[i]<- "red"}  
  
  waffle$cat[i] = plotdata$outcome[which(plotdata$dataset == d)]
}

Figure4A = ggplot(waffle)+
  geom_tile(aes(y=dataset, x= mortality, fill = fill),
            height = 0.9)+
  scale_fill_identity(breaks = c("royalblue",
                                 "black",
                                 "red"),
                      labels = c("< null models",
                                 "within null models",
                                 "> null models"),
                      guide = guide_legend(
                    title = "Observed mortality",
                    title.position = "top",
                    title.hjust = 0.5,
                    keywidth = 0.5,
                    keyheight = 0.5,
                    label.vjust = 0.5,
                  title.theme = element_text(face="italic", size=10)
                  
  ))+

  scale_x_continuous(name = "Mortality (%)") +
  coord_fixed(ratio = 1.2) +
  theme_bw()+
  theme(
    axis.text.y= element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1, 0, 1, 0), "cm")
  )



##model accuracy
#Figure 4B

meta_data = meta_data %>% unnest(cols = c(concentration.addition.list))

meta_data = meta_data %>% mutate(
  SA_dif = simple.addition - AB_mort,
  D_dif = dominance - AB_mort,
  M_dif = multiplicative - AB_mort,
  SAM_dif = stressor.addition - AB_mort,
  CA_dif = concentration.addition.list - AB_mort
)

accuracy = meta_data %>% summarise(
  bias_SA = mean( SA_dif), 
  prec_SA =sd(SA_dif),
  bias_D = mean( D_dif), 
  prec_D =sd(D_dif),
  bias_M = mean( M_dif), 
  prec_M =sd(M_dif),
  bias_SAM = mean(SAM_dif), 
  prec_SAM =sd(SAM_dif),
  bias_CA = mean(CA_dif), 
  prec_CA =sd(CA_dif))

accuracy = accuracy %>% pivot_longer(cols=everything())
bias=accuracy %>% filter(str_detect(name, "bias"))%>% mutate(yval = 5:1)
prec=accuracy %>% filter(str_detect(name, "prec")) %>% mutate(yval = 5:1)
bias$prec = prec$value

Figure4B = ggplot()+
  geom_violin(data=meta_data, aes(x=SA_dif, y = 5), fill="grey", adjust=2)+
  geom_violin(data=meta_data, aes(x=D_dif, y = 4), fill="grey", adjust=2)+
  geom_violin(data=meta_data, aes(x=M_dif, y = 3), fill="grey", adjust=2)+
  geom_violin(data=meta_data, aes(x=SAM_dif, y = 2), fill="grey", adjust=2)+
  geom_violin(data=meta_data, aes(x=CA_dif, y = 1), fill="grey", adjust=2)+
  geom_vline(aes(xintercept = 0), color="black", linetype="dashed")+
  geom_errorbarh(data=bias, aes(xmin=value-prec, xmax = value+prec,
                                y = yval, height=0.2), color ="black", size=0.5)+
  geom_point(data=bias, aes(x=value, y = yval), color ="black", size=2.5)+
  
   scale_x_continuous(limits=c(-101,101),
                     breaks = c(-100, -50,0,50,100),
                     labels=c("-100%", "-50%", "0%", "+50%", "+100%"),
                     name= "difference from observed mortality" )+
  scale_y_continuous(breaks=c(1,2,3,4,5), 
                     labels = c("Concentration addition", "Stressor addition",
                                "Multiplicative", "Dominance", "Simple addition"), name="")+
  theme_bw()+
  theme(plot.margin = unit(c(1, 1, 1, 0), "cm"))

require(cowplot)
F4 = plot_grid(Figure4A, Figure4B, labels = "AUTO", rel_widths = c(2, 2))
ggsave(F4, filename = "Figure4AB.jpg", device="jpg", height = 4, width = 7.5, unit="in")

require(gridExtra)
lay=rbind(c(1,1,1,2,2))
grid.arrange(Figure4A, Figure4B, layout_matrix=lay,
             padding=0)
