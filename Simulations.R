require(tidyverse)

# set up parameters for simulation
param <- expand_grid(
  SIA = (seq(0.01, 1, 0.01)),
  # stressor intensity 1
  SIB = (seq(0.01, 1, 0.01)),
  # stressor intensity 2
  fA = c(
    "linear",
    "accelerating",
    "diminishing",
    "steep_mid",
    "steep_ext"
  ),
  # stress-effect shapes
  fB = c(
    "linear",
    "accelerating",
    "diminishing",
    "steep_mid",
    "steep_ext"
  ),
  # stress-effect shapes
  strcap = c(3.2)
)

# create parameters for the stressor-effect curves
# based on a beta distribution

param$fAa <- param$fAb <- param$fBa <- param$fBb <- NA


# parameters for the stress-effect relationships
param$fAa <- ifelse(param$fA == "linear",
  1,
  ifelse(
    param$fA == "diminishing",
    3,
    ifelse(
      param$fA == "accelerating",
      1,
      ifelse(param$fA == "steep_ext", 3,
        ifelse(param$fA == "steep_mid", (1 /
          3), NA)
      )
    )
  )
)

param$fBa <- ifelse(param$fB == "linear",
  1,
  ifelse(
    param$fB == "diminishing",
    3,
    ifelse(
      param$fB == "accelerating",
      1,
      ifelse(param$fB == "steep_ext", 3,
        ifelse(param$fB == "steep_mid", (1 /
          3), NA)
      )
    )
  )
)

param$fAb <- ifelse(param$fA == "linear",
  1,
  ifelse(
    param$fA == "diminishing",
    1,
    ifelse(
      param$fA == "accelerating",
      3,
      ifelse(param$fA == "steep_ext", 3,
        ifelse(param$fA == "steep_mid", (1 /
          3), NA)
      )
    )
  )
)

param$fBb <- ifelse(param$fB == "linear",
  1,
  ifelse(
    param$fB == "diminishing",
    1,
    ifelse(
      param$fB == "accelerating",
      3,
      ifelse(param$fB == "steep_ext", 3,
        ifelse(param$fB == "steep_mid", (1 /
          3), NA)
      )
    )
  )
)


# These are the models that will be compared

param$SEA <- NA # Individual stress effect of stressor A
param$SEB <- NA # Individual stress effect of stressor B
param$SEsa <- NA # Joint Stress Effect - Simple Addition
param$SEea <- NA # Joint Stress Effect - Effect Addition (Multiplicative)
param$SEd <- NA # Joint Stress Effect - Dominance
param$SEca <- NA # Joint Stress Effect - Concentration Addition
param$SESAM <- NA # Joint Stress Effect - Stressor Addition Model

#run the simulations

for (i in 1:nrow(param)) {
  # for each parameter set

  param$SEA[i] <- qbeta(param$SIA[i], param$fAa[i], param$fAb[i])

  param$SEB[i] <- qbeta(param$SIB[i], param$fBa[i], param$fBb[i])

  # simple addition
  param$SEsa[i] <- param$SEA[i] + param$SEB[i]

  # effect addition
  param$SEea[i] <- param$SEA[i] + param$SEB[i] - param$SEA[i] * param$SEB[i]

  # dominance
  param$SEd[i] <- max(param$SEA[i], param$SEB[i])

  # concentration addition
  EC1 <- pbeta(0.5, param$fAa[i], param$fAb[i])
  EC2 <- pbeta(0.5, param$fBa[i], param$fBb[i])
  lambda <- as.numeric(EC1 / EC2)

  corrected_SI <- (param$SIA[i] + (lambda * param$SIB[i]))

  param$SEca[i] <- qbeta(
    ifelse(corrected_SI > 1, 1, corrected_SI),
    param$fAa[i],
    param$fAb[i]
  )

  ## Stressor Addition Model
  x1 <- qbeta(param$SEA[i], param$strcap[i], param$strcap[i])
  x2 <- qbeta(param$SEB[i], param$strcap[i], param$strcap[i])
  param$SESAM[i] <- pbeta((x1 + x2), param$strcap[i], param$strcap[i])

  print(round(100 * i / nrow(param), 1))
}

#####
#write_csv(param, "param.csv")
# write to file to use in the Sim
#####
