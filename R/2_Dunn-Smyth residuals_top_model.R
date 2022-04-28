###################################################################
##  title:Dunn-Smyth residuals for occupancy-detection models    ##
##  modified from Warton et al. 2017                             ##
##  authors: "Ana Patricia Calderon & Aimara Planillo"           ##
##  date: "23/11/2021"                                           ##
###################################################################

# Load libraries
library(camtrapR)
library(unmarked)
library(tibble)
library(RPresence)
library(mgcv)

###############################################
#                   Prepare data              #
###############################################

# Load data
data <- read.csv("./data/jaguar_data.csv", header = TRUE, row.names = 1)

# Format everything for RPresence
s <- nrow(data[,1:14]) # No. of sites (s).
K <- ncol(data[,1:14]) # No. of surveys/sampling occasions (K)

# Create df with site covariates
unitcov <- data.frame(int=rep(1,s),
                      elev= data$elev,
                      dist_sett= data$dist_sett,
                      dist_riv= data$dist_riv,
                      dist_road= data$dist_road,
                      dist_pa= data$dist_pa,
                      npp= data$npp,
                      treecov= data$treecov,
                      agriculture= data$agriculture)

# Create dataframe with detection covariates 
eff_p <- data[, 25:38]
eff_mat <- as.matrix(eff_p)
team <- data[, 39:52]
team_mat <- as.matrix(team)
jaguar <- data[,1:14]
survcov <- data.frame(int=rep(1,s*K), eff= c(eff_mat), team= c(team_mat))

# Create RPresence object
mydata <- createPao(data = jaguar, unitcov= unitcov, survcov=survcov)

#############################################################################
#         Run top model j63 in Rpresence and store model residuals          #
#############################################################################

topmod <- occMod(model=list(psi~ npp + dist_sett + dist_riv, p~ eff + team), data=mydata, type="so")
summary(topmod)

# Store model residuals
topmod_residuals <- residuals(topmod, type = "DS")

# Store residuals for occupancy
topmod_res_occ<- topmod_residuals$occ
length(topmod_res_occ)

# Store residuals for detection and remove NA values
topmod_res_det<- topmod_residuals$det[which(!is.na(topmod_residuals$det))] 
length(topmod_res_det)

######################################################
#         Plot model's Dunn-Smyth residuals          #
######################################################

# Check the normality of the Occupancy model residuals
qqnorm(topmod_res_occ,main="",col="blue",pch=16,cex.axis=1.4,ylab="",xlab="")
abline(0,1,col="black")
mtext(text="sample quantiles",side=2,line=2.2,cex=1.2)
title(main="Occupancy normality",cex.main=1)

# Check the normality of the detection model residuals
qqnorm(topmod_res_det,main="",col="blue",pch=16,cex.axis=1.4,ylab="",xlab="")
abline(0,1,col="black")
mtext(text="sample quantiles",side=2,line=2.2,cex=1.2)
title(main="Detection normality",cex.main=1)

# Store the fitted values for occupancy probabilities 
topmod_x_occ <- topmod$real$psi[,1] 

# Store fitted values for detection probabilities
topmod_x_det<- apply(matrix(topmod$real$p[,1],ncol=K),1,sum,na.rm=T)[is.na(topmod_residuals$det)==FALSE] 

# Function to create Dunn-Smyth residual plots
requireNamespace("mgcv")
dsResidPlot <- function (x, y, ylim = c(-1, 1) * max(abs(y)), alpha = 0.05,
                         k = 5, ylab = "", xlab = "", ...)
{
  requireNamespace("mgcv")
  plot(x, y, pch = 16, cex = 1.2, col = "blue", ylim = ylim,
       ylab = ylab, xlab = xlab, ...)
  lsmod <- gam(y ~ s(x, k = k))
  lsmod.p <- predict(lsmod, se.fit = TRUE)
  z.crit <- qnorm(1 - alpha/2)
  upr <- lsmod.p$fit + (z.crit * lsmod.p$se.fit)
  lwr <- lsmod.p$fit - (z.crit * lsmod.p$se.fit)
  polygon(c(rev(sort(x)), sort(x)), c(rev(upr[order(x)]), (lwr[order(x)])),
          col = "grey80", border = NA)
  points(x, y, pch = 16, cex = 1.2, col = "blue")
  lines(sort(x), fitted(lsmod)[order(x)], lwd = 2)
  abline(h = 0, col = "black", lwd = 1)
}

# Get residual plot for occupancy
dsResidPlot(topmod_x_occ,topmod_res_occ)
title(main="Occupancy Dunn-Smyth residuals",cex.main=1)
mtext(text="fitted occupancy values",side=1,las=1,line=2.6,cex=1)

# Get residual plot for detection
dsResidPlot(topmod_x_det,topmod_res_det)
title(main="Detection Dunn-Smyth residuals",cex.main=1)
mtext(text=Sigma[j]~"fitted detection values",side=1,las=1,line=2.7,cex=1.2)

# Plot detection model residuals against eff  

# To work with Dunn-Smyth residuals we need one mean value per site. We use the effort values for sites without NAs 
mean_eff_st <- rowMeans(eff_p, na.rm = TRUE)
length(mean_eff_st) #388 values, one per site
topmod_eff_toplot <- mean_eff_st[is.na(topmod_residuals$det)==FALSE]

# Plot residuals against the detection covariate 'effort'
dsResidPlot(topmod_eff_toplot,topmod_res_det)
title(main="Detection Dunn-Smyth residuals - Effort",cex.main=1)
mtext(text=Sigma[j]~"fitted detection values",side=1,las=1,line=2.7,cex=1)

########################################################################################
#                                         END OF SCRIPT                                #
########################################################################################

