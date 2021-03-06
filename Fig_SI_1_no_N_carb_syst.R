#----------------------------------------------------------------------------#
# Fig. SI.1 Main program                                                     #
# Modeling the carbonate system in the algal boundary layer; no N source     #
#                                                                            #
# M. Lavoie 2018                                                             #
#----------------------------------------------------------------------------#

rm(list=ls())

#------------------------------#
# Load packages                #
#------------------------------#

library(ReacTran)

#-----------------------------------------------------------------------------#
# Function calculating the chemical species of the carbonate system           #
# in the boundary layer.            N uptake rate = 0                         #

# This function computes the concentration (mol/m^3) of CO2, HCO3-, CO32-, H+ and OH- as a function of the distance 
# from the cell surface. This also computes relative changes in the boundary layer.
# I : ionic strengh in mol/L
# R : radius (m)
# L : length of the layer surrounding the cell (m)
# pH : negative logarithm of H+ activity far from the cell

#--------------------------------------------------------------------------#

boundary_carb_syst_No_N<- function(I, R, L, pH) {
  
  
  #------------------------------#
  #            MODEL             #
  #------------------------------#
  
  #--- Partial derivative equations
  
  boundary_layer <- function(t, state, parms) {
    
    { with( as.list( c(t, state, parms) ), {
      
      # Reshape state variable as a 2D matrix
      S    <- matrix(nrow = X.grid$N, ncol = 5, data = state)
      
      # Initialize dC/dt matrix
      dCdt <- 0*S
      
      # Rate of change for CO2
      tran_1 <- tran.1D(C = S[,1], D = D, flux.up = F1, C.down = C[1], A = A.grid, dx = X.grid, full.output = T)
      prod1 <- (kb[1] * S[,4] + kb[2]) * S[,2]
      loss1 <- (kf[1] + kf[2] * S[,5]) * S[,1]
      dCdt[,1] <- tran_1$dC + prod1 - loss1
      
      # Rate of change for HCO3 
      tran_2 = tran.1D(C = S[,2], D = D, C.down = C[2], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod2 = (kf[1] * S[,1]) + (kf[2] * S[,1] * S[,5]) + (kf[4] * S[,3] * S[,4]) + (kb[5] * S[,3]) 
      loss2 = (kb[1] * S[,4] * S[,2]) + (kb[2] * S[,2]) + (kb[4] * S[,2]) + (kf[5] * S[,2] * S[,5])
      dCdt[,2] = tran_2$dC + prod2 - loss2
      
      # Rate of change for CO3
      tran_3 = tran.1D(C = S[,3], D = D, C.down = C[3], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod3 = (kb[4] * S[,2]) + (kf[5] * S[,2] * S[,5]) 
      loss3 = (kf[4] * S[,4] * S[,3]) + (kb[5] * S[,3])
      dCdt[,3] = tran_3$dC + prod3 - loss3
      # flux.up = 0 since N uptake rate is set to 0.
      # Rate of change for H
      tran_4 = tran.1D(C = S[,4], D = D, C.down = C[4], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod4 = (kb[4] - kb[1] * S[,4]) * S[,2] + (kf[1] * S[,1]) + kb[6]  
      loss4 = (kf[4] * S[,4] * S[,3]) + (kf[6] * S[,4] * S[,5])
      dCdt[,4] = tran_4$dC + prod4 - loss4
      
      # Rate of change for OH
      tran_5 = tran.1D(C = S[,5], D = D, C.down = C[5], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod5 = (kb[2] * S[,2]) + (kb[5] * S[,3]) + kb[6]
      loss5 = (kf[2] * S[,1] * S[,5]) + (kf[5] * S[,2] * S[,5]) + (kf[6] * S[,4] * S[,5])
      dCdt[,5] = tran_5$dC + prod5 - loss5
      
      return(list(dCdt = dCdt))
    } )
    }
  }
  
  #------------------------------#
  # Model grid definition        #
  #------------------------------#
  
  # Number of grid layer
  N    <- 10000
  
  # Radius of the cell (m)
  R    <- R
  
  # Model grid setup 
  X.grid <- setup.grid.1D(x.up = R, L = L, N = N) # x.up = radius of the cell
  
  # Interface area
  A.grid <- setup.prop.1D(grid = X.grid, func = function(r) 4*pi*r^2)
  
  
  #------------------------------#
  # Model parameters             #
  #------------------------------#
  options(digits=15)
  
  # Calculation of CO2 concentration at I = 0
  pCO2 <- 10^-3.51    # CO2 partial presure (atm)
  Kh <- 0.034    # Henry's constant (mol atm-1 L-1) at I =0 ; Correction of Kh at higher I can be neglected at I typical of freshwaters (See Weiss, 1974)
  CO2 <- pCO2 * Kh * 1000     # CO2 concentration (mol m-3)
  
  # Fixed H+ and OH- concentrations at a given pH (mol m-3)
  source("func_activity_I.R") 
  coeff1 <- coeff_act(charge = 1, I = I)$coeff_act
  
  H <- (1/coeff1) * 1000 * (10^-pH)
  OH <- (1/coeff1) * 1000 * 1E-14 / (10^-pH)
  
  # Calculations of HCO3- concentration
  K1 <- 10^6.35              # Thermodynamic equilibrium constant of H+ + HCO3- = CO2
  source("func_K_corr2_I.R")  # Loading the function for ionic strength correction
  K1_corr <- K_corr2(K0=K1, I=I, chargeA=1, chargeB=1, chargeC=0)
  K1_corr
  HCO3 <- (CO2 * 1000 * (1/K1_corr$K1))/H
  
  # Calculations of CO32- concentration
  K2 <- 10^10.33             # Thermodynamic equilibrium constant of H+ + CO32- = HCO3-
  K2_corr <- K_corr2(K0=K2, I=I, chargeA=1, chargeB=2, chargeC=1)
  K2_corr
  CO3 <- (HCO3 * 1000 * (1/K2_corr$K1)) / H
  
  # Concentration of CO2,HCO3,CO3,H,OH in bulk solution (mol m^-3) at I=0 and T=25 C and a given pH
  C <- c(CO2, HCO3, CO3, H, OH)
  
  # Growth rate (u in d-1)
  u <- 1.12 * ((R*1E+06)^-0.75)
  
  # Uptake flux of CO2 by the cell (mol m-^2 s^-1)
  F1 <- -(23230/3) * R * u / 86400
  
  # Diffusion coefficient of CO2,... (m^2 s^-1)
  D  <- 1.18E-09
  
  # Rate constants at I=0
  kb_I0 <- c(7.88029840776055E+01, 1.8E-04, NA, 2.33867570643599, 6.54216399387684E+05, 1.4)  # m^3 mol^-1 s^-1; s^-1; ; s^-1; s^-1; mol m^-3 s^-1
  kf_I0 <- c(3.52E-02, 8.04030465871734, NA, 5E+07, 3.06E+06, 1.4E+08)   # s-^1; m^3 mol^-1 s^-1; ; m^3 mol^-1 s^-1; m^3 mol^-1 s^-1; m^3 mol^-1 s^-1
  
  source("func_k_rate_cst_corr2_I.R")  # Load the function converting rate constant at a given I
  kwat_I0 <- 1E-14       # Ion product of water at 25 ?C
  kwat <- kwat_I0 / (coeff1 * coeff1)   # Conditional ion product of water
  kb1 <- k_corr2(k0 = kb_I0[1], I = I, chargeA = 1, chargeB = 1, chargeC = 0)$k1 # HCO3- + H+ -> CO2
  kb1
  kf1 <- 1000* kb1/K1_corr$K1   # CO2 + H2O -> HCO3- + H+
  kf1
  kf2 <- k_corr2(k0 = kf_I0[2], I = I, chargeA = 0, chargeB = 1, chargeC = 1)$k1 # CO2 + OH- -> HCO3-
  kf2
  kb2 <- kf2 * K1_corr$K1 * kwat * 1000 # HCO3- -> CO2 + OH-
  kb2
  kf4 <- k_corr2(k0 = kf_I0[4], I = I, chargeA = 2, chargeB = 1, chargeC = 1)$k1 # CO32- + H+ -> HCO3-
  kf4
  kb4 <- kf4 * (1/(K2_corr$K1 / 1000))    # HCO3- -> CO32- + H+
  kb4
  kf5 <- k_corr2(k0 = kf_I0[5], I = I, chargeA = 1, chargeB = 1, chargeC = 2)$k1  # HCO3- + OH- <-> CO32- + H2O   
  kf5
  kb5 <- kf5 * kwat * 1000 * 1000 * K2_corr$K1 / 1000      # CO32- + H2O -> HCO3- + OH-
  kb5
  kf6 <- k_corr2(k0 = kf_I0[6], I = I, chargeA = 1, chargeB = 1, chargeC = 0)$k1  # H+ + OH- -> H2O    
  kf6
  kb6 <- kf6 * kwat * 1000 * 1000
  kb6
  
  # Rate constants at I = I
  kb <- c(kb1, kb2, NA, kb4, kb5, kb6)  
  kf <- c(kf1, kf2, NA, kf4, kf5, kf6)
  
  # Define parameters + grid definition vector
  parms <- list(C=C, F1=F1, D=D, kb=kb, kf=kf, X.grid=X.grid, A.grid=A.grid)
  
  
  #------------------------------#
  # Model solution               #
  #------------------------------#
  
  # Numerical solution at steady state 
  
  Cini <- matrix(C, nrow=N, ncol=5, byrow=T)
  # I added pos = TRUE
  boundary <- steady.1D(y = Cini, func = boundary_layer, pos = TRUE, atol=1E-8, parms = parms, nspec = 5, names = c('CO2','HCO3','CO3','H','OH'))
  
  
  #------------------------------#
  # Plotting output              #
  #------------------------------#
  
  # Using S3 plot method of package rootSolve'
  
  plotmult <- plot(boundary, grid = X.grid$x.mid, xlab = 'distance from centre, m', ylab = 'mol/m3', main = c('CO2', 'HCO3', 'CO3', 'H', 'OH'), ask = F, mfrow = c(1,1))
  
  # Storing X-axis label
  xlabel <- X.grid$x.mid
  
  # Relative enrichment calculations
  bmat <- unlist(boundary$y)
  bmat1 <- as.matrix(bmat)
  CO2_res <- bmat1[,1]
  CO2enrich <- CO2_res / CO2
  
  HCO3_res <- bmat1[,2]
  HCO3enrich <- HCO3_res / HCO3
  
  CO3_res <- bmat1[,3]
  CO3enrich <- CO3_res / CO3
  
  H_res <- bmat1[,4]
  Henrich <- H_res / H
  
  OH_res <- bmat1[,5]
  OHenrich <- OH_res / OH
  
  return(list(plotmult = plotmult, CO2enrich = CO2enrich, HCO3enrich = HCO3enrich, CO3enrich = CO3enrich, Henrich = Henrich, OHenrich = OHenrich, xlabel = xlabel))
  
}

#--------------------#
# Plots of FigSI.1   #
#--------------------#

#------------------------------------------------#
# Calculations at different pH at R = 5 um      #
#------------------------------------------------#
carb_pH7 <- boundary_carb_syst_No_N(I = 0.001, R = 5E-06 , L = 150E-06, pH = 7)
carb_pH5 <- boundary_carb_syst_No_N(I = 0.001, R = 5E-06 , L = 150E-06, pH = 5)
carb_pH8 <- boundary_carb_syst_No_N(I = 0.001, R = 5E-06 , L = 150E-06, pH = 8)

# Plot of relative enrichment (C/Co) of each chemical species at different pHs.
tiff( "FigS1_AE_no_N.tiff", res = 100)
oldpar <- par(mfrow=c(3,2), mar=c(0,0,0,3), oma = c(4,4.1,4,0.4), las=1)  

# Panel A : CO2 concentrations (r = 5 um)
plot(carb_pH7$xlabel*1E+06-5, carb_pH7$CO2enrich, type = "l", lty = 1, xaxt = "n", xlim = c(0, 60), ylim = c(0.9, 1))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$CO2enrich, type = "l", lty = 2)
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$CO2enrich, type = "l", lty = 3)
legend("bottomright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3))
mtext(text = "A", side = 3, adj = 0.05, line = -1.4, font = 2)
mtext(text = "relative change", side=2, line = 2.7, outer=TRUE, las=0)
mtext(text = expression(paste(" distance from cell surface (", mu, "m)")), side = 1, line = 3, font = 2, outer=TRUE, las=1)

# Panel B : HCO3- concentrations (r = 5 um)
plot(carb_pH7$xlabel*1E+06-5, carb_pH7$HCO3enrich, type = "l", lty = 1, xaxt = "n", xlim = c(0, 60), ylim = c(0.95, 1.01))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$HCO3enrich, type = "l", lty = 2) 
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$HCO3enrich, type = "l", lty = 3) 
mtext(text = "B", side = 3, adj = 0.05, line = -1.4, font = 2)

# Panel C : CO32- concentrations (r = 5 um)
plot(carb_pH7$xlabel*1E+06-5, carb_pH7$CO3enrich, type = "l", lty = 1, xaxt = "n", xlim = c(0, 60), ylim = c(0.95, 1.03))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$CO3enrich, type = "l", lty = 2) 
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$CO3enrich, type = "l", lty = 3) 
mtext(text = "C", side = 3, adj = 0.05, line = -2.1, font = 2)

# Panel D : H+ concentrations (r = 5 um)
plot(carb_pH7$xlabel*1E+06-5, carb_pH7$Henrich, type = "l", lty = 1, xlim = c(0, 60), ylim = c(0.95, 1.01))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$Henrich, type = "l", lty = 2) 
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$Henrich, type = "l", lty = 3) 
mtext(text = "D", side = 3, adj = 0.05, line = -1.4, font = 2)

# Panel E : OH- concentrations (r = 5 um)
plot(carb_pH7$xlabel*1E+06-5, carb_pH7$OHenrich, type = "l", lty = 1, xlim = c(0, 60), ylim = c(0.99, 1.03))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$OHenrich, type = "l", lty = 2)
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$OHenrich, type = "l", lty = 3)
mtext(text = "E", side = 3, adj = 0.05, line = -2.1, font = 2)

par(oldpar)
dev.off()

#------------------------------------------------#
# Calculations at different pH at R = 30 um      #
#------------------------------------------------#
carb_pH7 <- boundary_carb_syst_No_N(I = 0.001, R = 30E-06 , L = 900E-06, pH = 7)
carb_pH5 <- boundary_carb_syst_No_N(I = 0.001, R = 30E-06 , L = 900E-06, pH = 5)
carb_pH8 <- boundary_carb_syst_No_N(I = 0.001, R = 30E-06 , L = 900E-06, pH = 8)

# Plot of relative enrichment (C/Co) of each chemical species at different pHs.
tiff( "FigS1_F_no_N.tiff", res = 100)
oldpar <- par(mfrow=c(3,2), mar=c(0,0,0,3), oma = c(4.1,4,4,0.4), las=1) 

# Panel F : CO2 concentrations (r = 30 um)
plot(carb_pH7$xlabel*1E+06-30, carb_pH7$CO2enrich, type = "l", lty = 1, xaxt = "n", xlim = c(0, 100), ylim = c(0, 1))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$CO2enrich, type = "l", lty = 2) 
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$CO2enrich, type = "l", lty = 3)
legend("bottomright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3))
mtext(text = "F", side = 3, adj = 0.05, line = -1.4, font = 2)
mtext(text = "relative change", side=2, line = 2.6, outer=TRUE, las=0)
mtext(text = expression(paste(" distance from cell surface (", mu, "m)")), side = 1, line = 3, font = 2, outer=TRUE, las=1)

# Panel G : HCO3- concentrations (r = 30 um)
plot(carb_pH7$xlabel*1E+06-30, carb_pH7$HCO3enrich, type = "l", lty = 1, xaxt = "n", xlim = c(0, 100), ylim = c(0.5, 1.01))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$HCO3enrich, type = "l", lty = 2)
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$HCO3enrich, type = "l", lty = 3)
mtext(text = "G", side = 3, adj = 0.03, line = -2.1, font = 2)

# Panel H : CO32- concentrations (r = 30 um)
plot(carb_pH7$xlabel*1E+06-30, carb_pH7$CO3enrich, type = "l", lty = 1, xaxt = "n", xlim = c(0, 100), ylim = c(0.85, 5.1))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$CO3enrich, type = "l", lty = 2)
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$CO3enrich, type = "l", lty = 3)
mtext(text = "H", side = 3, adj = 0.03, line = -1.4, font = 2)

# Panel I : H+ concentrations (r = 30 um)
plot(carb_pH7$xlabel*1E+06-30, carb_pH7$Henrich, type = "l", lty = 1, xlim = c(0, 100), ylim = c(0, 1.1))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$Henrich, type = "l", lty = 2)
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$Henrich, type = "l", lty = 3)
mtext(text = "I", side = 3, adj = 0.03, line = -1.4, font = 2)

# Panel J : OH- concentrations (r = 30 um)
plot(carb_pH7$xlabel*1E+06-30, carb_pH7$OHenrich, type = "l", lty = 1, xlim = c(0, 100), ylim = c(1, 5))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$OHenrich, type = "l", lty = 2)
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$OHenrich, type = "l", lty = 3)
mtext(text = "J", side = 3, adj = 0.03, line = -1.4, font = 2)

par(oldpar)
dev.off()

#-----------------------------#
# References                  #
#-----------------------------#
# Martell, A. E., Smith, R. M., Motekaitis, R. J. 2004. NIST critical stability constants of metal complexes, version 8. National Institute of Standards and Technology. Gaithersburg, MD. In Gaithersburg, MD, 2004.
# Weiss, R.F. 1974. Carbon dioxide in water and seawater : The solubility of a non-ideal gas. Marine Chemistry 2 : 203-215