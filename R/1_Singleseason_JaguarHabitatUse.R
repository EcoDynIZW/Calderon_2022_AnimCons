######################################################################
#  title: "Jaguar Habitat Use in Central America"                    #
#  authors: "Ana Patricia Calderon & Julie Louvrier"                  #
#  date: "23/11/2021"                                                #
######################################################################

# Load libraries
library(camtrapR)
library(unmarked)
library(tibble)
library(Hmisc)
library(corrplot)
library(ubms)

rm(list = ls())

###############################
#   Load and prepare data     #
###############################

data <- read.csv("./data/jaguar_data.csv", header = TRUE, row.names = 1)

Jaguar2 <- unmarkedFrameOccu(y=data[,1:14],siteCovs=data[,15:24],
                             obsCovs=list(eff=data[,25:38], team=data[,39:52]))

###############################
#     Set detection first     #
###############################

jnull <- occu(~1 ~1, data = Jaguar2)
j46 <- occu(~ eff ~ 1, data=Jaguar2)
j47 <- occu(~ team ~ 1, data=Jaguar2)
j48 <- occu(~ team + eff ~ 1, data=Jaguar2)

modelList.det2 <- fitList(
  jnull= jnull,
  j46 = j46,
  j47 = j47,
  j48 = j48
)

class(modelList.det2)
(modSel.det2 <- modSel(modelList.det2))

# Top detection model j48

###############################
#   Natural habitat models    #
###############################
j49 <- occu(~ team + eff ~ treecov + npp, data=Jaguar2)
j50<- occu(~ team + eff ~ treecov + dist_riv, data=Jaguar2)
j51 <- occu(~ team + eff ~ treecov + agriculture + treecov*agriculture, data=Jaguar2)
j52 <- occu(~ team + eff ~ npp + dist_riv, data=Jaguar2)
j53 <- occu(~ team + eff ~ npp + elev, data=Jaguar2)
j54 <- occu(~ team + eff ~ npp + elev + dist_riv, data=Jaguar2)

###############################
#   Human influence models    #
###############################
j55 <- occu(~ team + eff~ dist_sett, data=Jaguar2)
j56 <- occu(~ team + eff ~ dist_sett + dist_road, data=Jaguar2)
j57 <- occu(~ team + eff ~ dist_sett + dist_road + dist_pa, data=Jaguar2)
j58 <- occu(~ team + eff ~ dist_pa, data=Jaguar2)
j59 <- occu(~ team + eff ~ dist_sett + dist_pa, data=Jaguar2)
j60 <- occu(~ team + eff ~ dist_sett + dist_pa + dist_sett*dist_pa, data=Jaguar2)
j61 <- occu(~ team + eff ~ dist_sett + dist_road + dist_pa + dist_sett*dist_pa, data=Jaguar2)

#######################################################
#     Natural habitat and Human influence models      #
#######################################################
j62 <- occu(~ team + eff ~ treecov + dist_sett + dist_riv, data=Jaguar2)
j63 <- occu(~ team + eff ~ npp + dist_sett + dist_riv, data=Jaguar2)
j69 <- occu(~ team + eff ~ treecov + dist_sett + dist_riv + elev + dist_road + dist_pa, data=Jaguar2)
j71 <- occu(~ team + eff ~ npp + treecov  + dist_sett + dist_riv + agriculture + treecov*agriculture, data=Jaguar2)
j73 <- occu(~ team + eff ~ elev + dist_sett + dist_road + dist_riv + dist_pa + dist_sett*dist_pa, data=Jaguar2)

# Models including dist_road, treecov and agriculture were removed from the analysis because these covariates were uninformative (sensu Arnold 2010) 

modelList.occu1 <- fitList(
  jnull=jnull,
  j48=j48,
  j49=j49,
  j50=j50,
  j51=j51,
  j52=j52,
  j53=j53,
  j54=j54,
  j55=j55,
  j56=j56,
  j57=j57,
  j58=j58,
  j59=j59,
  j60=j60,
  j61=j61,
  j62=j62,
  j63=j63,
  j69=j69,
  j71=j71,
  j73=j73)  

(modSel.occu1 <- modSel(modelList.occu1))

# Top model was j63 with 0.89 AIC cumulative weight
# j63 <- occu(~ team + eff ~ npp + dist_sett + dist_riv, data=Jaguar2)


##############################################################################
##        Running candidate models with year as random effect              ##
##############################################################################

# Evaluate potential effects of year on jaguar occupancy

# Install.dist_packages("devtools")
# devtools::install_github("kenkellner/ubms")

# Install dependencies
# install.dist_packages(c("unmarked", "ggplot2", "gridExtra", "lme4", "loo","Matrix", "Rcpp", "rstan", "rstantools"))

# Download and install ubms
# download.file("https://github.com/kenkellner/ubms/releases/download/v1.0.2/ubms_1.0.2.zip", destfile="ubms_1.0.2.zip")
# install.dist_packages("ubms_1.0.2.zip", repos=NULL)

# Best model: j63 <- occu(~ team + eff ~ npp + dist_sett + dist_riv, data=Jaguar2)
options(mc.cores=3)

###################################
#     Set detection first         #
###################################

# Models with no random effect
jnull <- stan_occu(~1 ~1, data = Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j46 <- stan_occu(~ eff ~ 1, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j47 <- stan_occu(~ team ~ 1, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j48 <- stan_occu(~ team + eff ~ 1, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))

# Models with random effect
jnull_random <- stan_occu(~1 ~1 + (1|year), data = Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j46_random <- stan_occu(~ eff ~ 1 + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j47_random <- stan_occu(~ team ~ 1 + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j48_random<- stan_occu(~ team + eff ~ 1 + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))

# Model selection
modelList.det2 <- fitList(
  jnull= jnull ,
   j46 = j46,
   j47 = j47,
   j48 = j48,
   jnull_random = jnull_random,
   j46_random = j46_random,
   j47_random = j47_random,
   j48_random = j48_random)

class(modelList.det2)
(modSel.det2 <- modSel(modelList.det2))

#                   elpd    ndist_param    elpd_diff   se_diff   weight
# j48          -745.3063 14.703009   0.0000000  0.000000 3.686547e-01
# j48_random   -745.5100 20.166959  -0.2036451  1.888011 4.374888e-01

# elpd is the expected log point predictive density. The "top" ranked model (with highest/least negative elpd) is first in the table, and all other models are being compared to the top model in pairwise comparisons. Models j48 and j48_random are equivalent and best

#################################
#       Set Occupancy           #
#################################

###############################
#   Natural habitat models    #
###############################

# Models with no random effect
j49 <- stan_occu(~ team + eff ~ treecov + npp , data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j50<- stan_occu(~ team + eff ~ treecov + dist_riv , data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j52 <- stan_occu(~ team + eff ~ npp + dist_riv , data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j53 <- stan_occu(~ team + eff ~ npp + elev , data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j54 <- stan_occu(~ team + eff ~ npp + elev + dist_riv , data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))

# Models with random effect
j49_random <- stan_occu(~ team + eff ~ treecov + npp + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j50_random <- stan_occu(~ team + eff ~ treecov + dist_riv + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j52_random <- stan_occu(~ team + eff ~ npp + dist_riv + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j53_random <- stan_occu(~ team + eff ~ npp + elev + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j54_random <- stan_occu(~ team + eff ~ npp + elev + dist_riv + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
?stan_occu

###############################
#   Human influence models    #
###############################

# Models with no random effect
j55 <- stan_occu(~ team + eff~ dist_sett, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j56 <- stan_occu(~ team + eff ~ dist_sett + dist_road, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j57 <- stan_occu(~ team + eff ~ dist_sett + dist_road + dist_pa, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j58 <- stan_occu(~ team + eff ~ dist_pa, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j59 <- stan_occu(~ team + eff ~ dist_sett + dist_pa, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j60 <- stan_occu(~ team + eff ~ dist_sett + dist_pa + dist_sett*dist_pa, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j61 <- stan_occu(~ team + eff ~ dist_sett + dist_road + dist_pa + dist_sett*dist_pa, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))

# Models with random effect
j55_random <- stan_occu(~ team + eff~ dist_sett + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j56_random <- stan_occu(~ team + eff ~ dist_sett + dist_road + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j57_random <- stan_occu(~ team + eff ~ dist_sett + dist_road + dist_pa + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j58_random <- stan_occu(~ team + eff ~ dist_pa + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j59_random <- stan_occu(~ team + eff ~ dist_sett + dist_pa + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j60_random <- stan_occu(~ team + eff ~ dist_sett + dist_pa + dist_sett*dist_pa + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j61_random <- stan_occu(~ team + eff ~ dist_sett + dist_road + dist_pa + dist_sett*dist_pa + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))

#######################################################
#     Natural habitat and Human influence models      #
#######################################################

# Models with no random effect
j51 <- stan_occu(~ team + eff ~ treecov + agriculture + treecov*agriculture , data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j62 <- stan_occu(~ team + eff ~ treecov + dist_sett + dist_riv, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j63 <- stan_occu(~ team + eff ~ npp + dist_sett + dist_riv, data=Jaguar2, chains=3, iter=30000, control = list(adapt_delta = 0.99))
j69 <- stan_occu(~ team + eff ~ treecov + dist_sett + dist_riv + elev + dist_road + dist_pa, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j73 <- stan_occu(~ team + eff ~ elev + dist_sett + dist_road + dist_riv + dist_pa + dist_sett*dist_pa, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j71 <- stan_occu(~ team + eff ~ npp + treecov  + dist_sett + dist_riv + agriculture + treecov*agriculture, data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))

# Models with random effect
j51_random <- stan_occu(~ team + eff ~ treecov + agriculture + treecov*agriculture + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j62_random <- stan_occu(~ team + eff ~ treecov + dist_sett + dist_riv + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j63_random1 <- stan_occu(~ team + eff ~ npp + dist_sett + dist_riv+ (1|year), data=Jaguar2 , chains=3, iter=10000, control = list(adapt_delta = 0.99))
j69_random <- stan_occu(~ team + eff ~ treecov + dist_sett + dist_riv + elev + dist_road + dist_pa + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j71_random <- stan_occu(~ team + eff ~ npp + treecov  + dist_sett + dist_riv + agriculture + treecov*agriculture + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))
j73_random <- stan_occu(~ team + eff ~ elev + dist_sett + dist_road + dist_riv + dist_pa + dist_sett*dist_pa + (1|year), data=Jaguar2, chains=3, iter=10000, control = list(adapt_delta = 0.99))

modelList.occu1 <- fitList(
  jnull=jnull,
  j48=j48,
  j49=j49,
  j50=j50,
  j51=j51,
  j52=j52,
  j53=j53,
  j54=j54,
  j55=j55,
  j56=j56,
  j57=j57,
  j58=j58,
  j59=j59,
  j60=j60,
  j61=j61,
  j62=j62,
  j63=j63,
  j69=j69,
  j71=j71,
  j73=j73,
  jnull_random=jnull_random,
  j48_random=j48_random,
  j49_random=j49_random,
  j50_random=j50_random,
  j51_random=j51_random,
  j52_random=j52_random,
  j53_random=j53_random,
  j54_random=j54_random,
  j55_random=j55_random,
  j56_random=j56_random,
  j57_random=j57_random,
  j58_random=j58_random,
  j59_random=j59_random,
  j60_random=j60_random,
  j61_random=j61_random,
  j62_random=j62_random,
  j63_random1=j63_random1,
  j69_random=j69_random,
  j71_random=j71_random,
  j73_random=j73_random
  )

(modSel.occu1 <- modSel(modelList.occu1))

#                   elpd    ndist_param   elpd_diff    se_diff       weight
# j63          -706.7618 17.146911    0.000000  0.0000000 1.472172e-01
# j63_random   -707.8232 19.772510   -1.061392  0.9674488 4.140711e-02

# The top ranked model does not contain year as random effect and it matches the results in the previous run. We evaluate the goodness of fit of this top ranked model (j63)

# MacKenzie-Bailey Chi-square goodness of fit for j63 model 
gof_63 <- gof(j63)
# The model has a good fit: Chi-square = 15891.919, posterior predictive p = 0.925

# Leave one out cross validation
loo_j63 <- loo(j63)
# All Pareto k estimates are ok (k < 0.7).

########################################################################################
#                                         END OF SCRIPT                                #
########################################################################################
