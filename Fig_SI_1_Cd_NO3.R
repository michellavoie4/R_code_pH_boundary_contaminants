#-----------------------------------------------------------------------------#
# Fig. SI.1 Main program                                                      #
#                                                                             #
#                                                                             #
# Lavoie M. 2018                                                              #
#-----------------------------------------------------------------------------#

#------------------------------#
# Load packages                #
#------------------------------#

library(ReacTran)

#-----------------------------------------------------------------------------#
# Function calculating Cd speciation in the boundary layer 
# of chrysophyte cells.     
# Inclusion of the Cd 2+, Cd(OH)+, CdCO3 
# NO3- is the N source and R = 5 and 30 um ; pH = 7, 5 and 8

# This function computes the concentration (mol/m^3) of CO2, HCO3-, CO32-, H+, OH-, Cd2+, CdOH+ and CdCO3 as a function of the distance 
# from the cell surface assuming NO3- is the N source. This computes also relative changes in the boundary layer
# I : ionic strengh in mol/L
# R : radius (m)
# L : length of the layer surrounding the cell (m)
# pH : negative logarithm of H+ activity far from the cell

#--------------------------------------------------------------------------#

boundary_Cd_NO3<- function(I, R, L, pH) {
  
  
  #-----------------------------------------------------------------------------#
  # Part 1 : Bjerrum's ion-association model and metal speciation modeling      #
  #                                                                             #
  #-----------------------------------------------------------------------------#
  
  #------------#
  # Parameters #
  #------------#
  
  T <- 298    # Temperature in Kelvins (K)
  I <- I              # Ionic strength (mol/Kg)
  
  
  # Calculating Kos for Cadmium and OH-  (Cd2+ + OH- = CdOH-)
  Z1 <- 2                 # Absolute electric charge of the metal
  Z2 <- 1                 # Absolute electric charge of the ligand
  ion_rad <- 95           # Effective ionic radius of the metal (pm) (Marcus, 1988)
  water_diam <- 250       # Diameter of a water molecule (pm) (Shatzberg, 1967)
  ligand_rad <- 110       # Effective ionic radius of the ligand (pm)  (Sethi and Raghavan, 1988)
  a <- 1E-10* (ion_rad + water_diam + ligand_rad)  # minimum approach distance between the charge ions (cm)
  source('~/R software/R scripts/Functions/func_Bjerrum_ion_assoc.R')   # load the function calculating Kos using the Bjerrum's ion association model
  Kos_CdOH <- fun_Kos(Z1=Z1, Z2=Z2, a=a, T=T)$Kos
  
  # Calculating Kos for Cd and CO32-  (Cd2+ + CO32- = CdCO3)
  Z1 <- 2                 # Absolute electric charge of the metal
  Z2 <- 2                 # Absolute electric charge of the ligand
  ion_rad <- 95           # Effective ionic radius of the metal (pm) (Marcus, 1988)
  water_diam <- 250       # Diameter of a water molecule (pm) (Shatzberg, 1967)
  ligand_rad <- 300       # Effective ionic radius of the ligand (pm) (Sethi and Raghavan, 1998)
  b <- 1E-10* (ion_rad + water_diam + ligand_rad)  # minimum approach distance between the charge ions (cm)
  source('~/R software/R scripts/Functions/func_Bjerrum_ion_assoc.R')   # load the function calculating Kos using the Bjerrum's ion association model
  Kos_CdCO3 <- fun_Kos(Z1=Z1, Z2=Z2, a=b, T=T)$Kos
  
  # Calculating forward and backward rate constants for the equation : Cd2+ + OH- = CdOH+
  k_w <- 3E+08        # Water loss rate constants (in s-1) for Cd2+ (for the formation of unidentate complexes) (Stumm and Morgan, 1996; Wilkinson et al 2004)
  kf_CdOH <- Kos_CdOH * k_w
  kf_CdOH
  source("~/R software/R scripts/Functions/func_K_corr2_I.R")  # Loading the function for ionic strength correction of equilibrium constant (K0 to K1)
  K0_CdOH <- 10^3.91   # Thermodynamic equlibrium constant for CdOH+ formation
  K1_CdOH <- K_corr2(K0=K0_CdOH, I=I, chargeA=2, chargeB=1, chargeC=1)$K1
  K1_CdOH             # Thermodynamic equilibrium constant corrected at I = I
  log10(K1_CdOH)
  kb_CdOH <- kf_CdOH / K1_CdOH
  kb_CdOH
  
  # Calculating kf and kb for the equation : Cd2+ + CO32- = CdCO3
  k_w <- k_w    # Water loss rate constants (in s-1) for Cd2+ (for the aquo-ions)
  kf_CdCO3 <- Kos_CdCO3 * k_w     # in L mol-1 s-1
  kf_CdCO3
  source("~/R software/R scripts/Functions/func_K_corr2_I.R")  # Loading the function for ionic strength correction of equilibrium constant (K0 to K1)
  K0_CdCO3 <- 10^4.36   # Thermodynamic equlibrium constant for CdCO3 formation
  K1_CdCO3 <- K_corr2(K0=K0_CdCO3, I=I, chargeA=2, chargeB=2, chargeC=0)$K1
  K1_CdCO3
  log10(K1_CdCO3)
  kb_CdCO3 <- kf_CdCO3 / K1_CdCO3   # in s-1
  kb_CdCO3
  
  # Converting units in m^3
  kf_CdOH <- kf_CdOH / 1000       # Convert in m^3 mol-1 s-1
  kb_CdOH <- kb_CdOH              # no conversion (s-1)
  kf_CdCO3 <- kf_CdCO3 / 1000     # Convert in m^3 mol-1 s-1
  kb_CdCO3                        # no conversion (s-1)
  K1_CdOH <- K1_CdOH / 1000       # Convert in m^3 mol-1
  K1_CdCO3 <- K1_CdCO3 / 1000     # Convert in m^3 mol-1
  
  #-----------------------------------------------------------------------------#
  # Part 2 : Reaction diffusion modeling                                                                             #
  #-----------------------------------------------------------------------------#
  
  #------------------------------#
  #            MODEL             #
  #------------------------------#
  
  #--- Partial derivative equations
  
  boundary_layer <- function(t, state, parms) {
    
    { with( as.list( c(t, state, parms) ), {
      
      # Reshape state variable as a 2D matrix
      S    <- matrix(nrow = X.grid$N, ncol = 8, data = state)
      
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
      prod3 = (kb[4] * S[,2]) + (kf[5] * S[,2] * S[,5]) + (kb_CdCO3 * S[,8])
      loss3 = (kf[4] * S[,4] * S[,3]) + (kb[5] * S[,3]) + (kf_CdCO3 * S[,6] * S[,3])
      dCdt[,3] = tran_3$dC + prod3 - loss3
      # flux.up = NO3_H ou NH4_H
      # Rate of change for H
      tran_4 = tran.1D(C = S[,4], D = D, C.down = C[4], flux.up = NO3_H, A = A.grid, dx = X.grid, full.output = T)
      prod4 = (kb[4] - kb[1] * S[,4]) * S[,2] + (kf[1] * S[,1]) + kb[6]  
      loss4 = (kf[4] * S[,4] * S[,3]) + (kf[6] * S[,4] * S[,5])
      dCdt[,4] = tran_4$dC + prod4 - loss4
      
      # Rate of change for OH
      tran_5 = tran.1D(C = S[,5], D = D, C.down = C[5], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod5 = (kb[2] * S[,2]) + (kb[5] * S[,3]) + kb[6] + (kb_CdOH * S[,7])
      loss5 = (kf[2] * S[,1] * S[,5]) + (kf[5] * S[,2] * S[,5]) + (kf[6] * S[,4] * S[,5]) + (kf_CdOH * S[,6] * S[,5])
      dCdt[,5] = tran_5$dC + prod5 - loss5
      
      # Rate of change for Cd2+
      tran_6 = tran.1D(C = S[,6], D = D, C.down = C[6], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod6 = (kb_CdOH * S[,7]) + (kb_CdCO3 * S[,8])
      loss6 = (kf_CdOH * S[,6] * S[,5]) + (kf_CdCO3 * S[,6] * S[,3])
      dCdt[,6] = tran_6$dC + prod6 - loss6
      
      # Rate of change for CdOH
      tran_7 = tran.1D(C = S[,7], D = D, C.down = C[7], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod7 = (kf_CdOH * S[,6] * S[,5])
      loss7 = (kb_CdOH * S[,7])
      dCdt[,7] = tran_7$dC + prod7 - loss7
      
      # Rate of change for CdCO3
      tran_8 = tran.1D(C = S[,8], D = D, C.down = C[8], flux.up = 0, A = A.grid, dx = X.grid, full.output = T)
      prod8 = (kf_CdCO3 * S[,6] * S[,3])
      loss8 = (kb_CdCO3 * S[,8])
      dCdt[,8] = tran_8$dC + prod8 - loss8
      
      return(list(dCdt = dCdt))
    } )
    }
  }
  
  #------------------------------#
  # Model grid definition        #
  #------------------------------#
  
  # Number of grid layer
  N    <- 10000
  
  # Radius of the cell or boundary layer thickness (m)
  R    <- R 
  
  # Model grid setup 
  X.grid <- setup.grid.1D(x.up = R, L = L, N = N) # x.up = radius of the cell; radius = boundary layer thickness
  
  # Interface area
  A.grid <- setup.prop.1D(grid = X.grid, func = function(r) 4*pi*r^2)
  
  
  #------------------------------#
  # Model parameters             #
  #------------------------------#
  options(digits=15)
  I <- I        # Ionic strength (M)
  
  # Concentration of CO2,HCO3,CO3,H,OH in bulk solution (mol m^-3) at I=0 25 °C and pH = 7
  #C_I0  <- c(1.05070044705462E-02, 4.69330649964981E-02, 2.19522437871783E-05, 1E-04, 1E-04)
  
  # Calculation of CO2 concentration at I = 0 or I (independant of I)
  pCO2 <- 10^-3.51    # CO2 partial presure (atm)
  Kh <- 0.034    # Henry's constant (mol atm-1 L-1) at I =0 ; Correction of Kh at higher I can be neglected at low I (See Weiss, 1974)
  CO2 <- pCO2 * Kh * 1000     # CO2 concentration (mol m-3)
  
  # Fixed H+ and OH- concentrations at a given pH (mol m-3) and I = I
  source("~/R software/R scripts/Functions/func_activity_I.R") 
  coeff1 <- coeff_act(charge = 1, I = I)$coeff_act
  
  H <- (1/coeff1) * 1000 * (10^-pH)
  OH <- (1/coeff1) * 1000 * 1E-14 / (10^-pH)
  
  # Calculations of HCO3- concentration at I = I (mol m-3)
  K1 <- 10^6.35              # Thermodynamic euilibrium constant of H+ + HCO3- = CO2
  source("~/R software/R scripts/Functions/func_K_corr2_I.R")  # Loading the function for ionic strength correction
  K1_corr <- K_corr2(K0=K1, I=I, chargeA=1, chargeB=1, chargeC=0)
  K1_corr
  HCO3 <- (CO2 * 1000 * (1/K1_corr$K1))/H
  
  # Calculations of CO32- concentration at I = I (mol m-3)
  K2 <- 10^10.33             # Thermodynamic equilibrium constant of H+ + CO32- = HCO3-
  K2_corr <- K_corr2(K0=K2, I=I, chargeA=1, chargeB=2, chargeC=1)
  K2_corr
  CO3 <- (HCO3 * 1000 * (1/K2_corr$K1)) / H
  
  # Calculations of Cd, CdOH- , CdCO3 concentration at I = I (mol m-3)
  # Cdtot <- Cd + CdOH + CdCO3
  # K1_CdOH <- CdOH / (Cd * OH)    ------>  Cd * K1_CdOH - (1/OH ) * CdOH = 0
  # K1_CdCO3 <- CdCO3 / (Cd * CO3) --------> Cd * K1_CdCO3 - (1/CO3) * CdCO3 = 0
  #
  # Thus, the system matrix :    Cd +    CdOH         + CdCO3           = Cdtot
  #                      K1_CdOH Cd - (1/OH) CdOH     + 0               = 0
  #                    K1_CdCO3  Cd +   0           - (1/CO3) * CdCO3   = 0 
  Cdtot <- 1E-09 * 1000         # in mol m-3
  mat1 <- matrix(c(1, 1, 1, K1_CdOH, -1/OH, 0, K1_CdCO3, 0, -1/CO3), nrow=3, ncol=3, byrow = T)
  mat2 <- matrix(c(Cdtot, 0, 0), nrow=3, ncol=1)
  # Assuming that CdCO3 and CdOH << free CO32- and free OH, respectively
  # i.e., Assuming no depletion in free OH and free CO32-
  out <- solve(mat1,mat2)
  Cd <- out[1]
  CdOH <- out[2]
  CdCO3 <- out[3]
  
  # Concentration of CO2, HCO3-,CO3,H+,OH-, CdOH+, CdCO3 in bulk solution (mol m^-3) at I = I, T = 25 °C and a given pH
  C <- c(CO2, HCO3, CO3, H, OH, Cd, CdOH, CdCO3)
  
  # Growth rate (u in d-1)
  u <- 1.12 * ((R*1E+06)^-0.75)
  
  # Uptake flux of CO2 by the cell (mol m-^2 s^-1)
  F1 <- -(23230/3) * R * u / 86400
  
  # Acido-basic reactions and NO3- assimilation
  # With a Redfield ratio C:N of 106: 16, for each mol C assimilated, 0.15 mol H+ is removed from the medium if NO3- is the N source 
  NO3_H <- 0.15 * F1
  
  # Acido-basic reactions and NH4+ assimilation
  # With a Redfield ratio C:N of 106: 16, for each mol C assimilated, 0.15 mol H+ is produced in the medium if NH4+ is the N source 
  NH4_H <- -0.15 * F1
  
  # Diffusion coefficient of CO2,... (m^2 s^-1)
  D  <- 1.18E-09
  
  # Rate constants at I=0
  kb_I0 <- c(7.88029840776055E+01, 1.8E-04, NA, 6.13540983937499, 6.54216399387684E+05, 1.4)   # m^3 mol^-1 s^-1; s^-1; ; s^-1; s^-1; mol m^-3 s^-1
  kf_I0 <- c(3.52E-02, 8.04030465871734, NA, 1.31172736401427E+08, 3.06E+06, 1.4E+08)          # s-^1; m^3 mol^-1 s^-1; ; m^3 mol^-1 s^-1; m^3 mol^-1 s^-1; m^3 mol^-1 s^-1
  
  source('~/R software/R scripts/Functions/func_k_rate_cst_corr2_I.R')  # Load the function converting rate constant at a given I
  kwat_I0 <- 1E-14       # Ion product of water at 25 °C
  kwat <- kwat_I0 / (coeff1 * coeff1)    # Conditional ion product of water
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
  parms <- list(NO3_H=NO3_H, C=C, F1=F1, D=D, kb=kb, kf=kf, X.grid=X.grid, A.grid=A.grid)
  
  
  #------------------------------#
  # Model solution               #
  #------------------------------#
  
  # Numerical solution at steady state 
  
  Cini <- matrix(C, nrow=N, ncol=8, byrow=T)
  
  boundary <- steady.1D(y = Cini, func = boundary_layer, pos = TRUE, atol=1E-8, parms = parms, nspec = 8, names = c('CO2','HCO3','CO3','H','OH','Cd','CdOH','CdCO3'))
  
  
  #------------------------------#
  # Plotting output              #
  #------------------------------#
  
  # Using S3 plot method of package rootSolve'
  
  plotmult <- plot(boundary, grid = X.grid$x.mid, xlab = 'distance from centre, m', ylab = 'mol/m3', main = c('CO2', 'HCO3', 'CO3', 'H', 'OH', 'Cd', 'CdOH', 'CdCO3' ), ask = F, mfrow = c(1,1))
  
  # Storing X-axis label
  xlabel <- X.grid$x.mid
  
  # Relative enrichment calculations
  bmat <- unlist(boundary$y)
  bmat1 <- as.matrix(bmat)
  Cd_res <- bmat1[,6]
  Cdenrich <- Cd_res / Cd
  
  CdOH_res <- bmat1[,7]
  CdOHenrich <- CdOH_res / CdOH
  
  CdCO3_res <- bmat1[,8]
  CdCO3enrich <- CdCO3_res / CdCO3
  
  return(list(plotmult = plotmult, Cdenrich = Cdenrich, CdOHenrich = CdOHenrich, CdCO3enrich = CdCO3enrich, xlabel = xlabel))
  
}

#--------------------#
# Plots of Fig3      #
#--------------------#

#------------------------------------------------#
# Calculations at different pH at R = 5 um      #
#------------------------------------------------#
carb_pH7 <- boundary_Cd_NO3(I = 0.001, R = 5E-06 , L = 150E-06, pH = 7)
carb_pH5 <- boundary_Cd_NO3(I = 0.001, R = 5E-06 , L = 150E-06, pH = 5)
carb_pH8 <- boundary_Cd_NO3(I = 0.001, R = 5E-06 , L = 150E-06, pH = 8)

# Plot of relative change (C/Co) of each Cu species at different pHs.
tiff( "Fig3_A.tiff", res = 72)
oldpar <- par(mfrow=c(2,2)) #, mar=c(4,4,2,1), oma = c(0,0,0,0)) #, mgp=c(1,0,0))
xlab <- expression(paste(" distance from cell surface (", mu, "m)"))
plot(carb_pH7$xlabel*1E+06-5, carb_pH7$Cdenrich, type = "l", lty = 1, xlim = c(0, 60), ylim = c(0.95, 1.05), xlab = xlab, ylab = "relative change", main = expression(paste("A : Cd"^"2+", "  R = 5", mu, "m")))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$Cdenrich, type = "l", lty = 2) # pointillé
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$Cdenrich, type = "l", lty = 3) # points
legend("topright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3)) #, y.intersp = 0.2, bty = "n")

plot(carb_pH7$xlabel*1E+06-5, carb_pH7$CdOHenrich, type = "l", lty = 1, xlim = c(0, 60), ylim = c(0.7, 1.2), xlab = xlab, ylab = "relative change", main = expression(paste("B : CdOH"^"+", "  R = 5", mu, "m")))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$CdOHenrich, type = "l", lty = 2) # pointillé
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$CdOHenrich, type = "l", lty = 3) # points
legend("bottomright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3)) #, y.intersp = 0.2, bty = "n")

plot(carb_pH7$xlabel*1E+06-5, carb_pH7$CdCO3enrich, type = "l", lty = 1, xlim = c(0, 60), ylim = c(0.7, 1.2), xlab = xlab, ylab = "relative change", main = expression(paste("C : CdCO"[3], "  R = 5", mu, "m")))
lines(carb_pH5$xlabel*1E+06-5, carb_pH5$CdCO3enrich, type = "l", lty = 2) # pointillé
lines(carb_pH8$xlabel*1E+06-5, carb_pH8$CdCO3enrich, type = "l", lty = 3) # points
legend("bottomright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3)) #, y.intersp = 0.2, bty = "n")

par(oldpar)
dev.off()

#------------------------------------------------#
# Calculations at different pH at R = 30 um      #
#------------------------------------------------#
carb_pH7 <- boundary_Cd_NO3(I = 0.001, R = 30E-06 , L = 900E-06, pH = 7)
carb_pH5 <- boundary_Cd_NO3(I = 0.001, R = 30E-06 , L = 900E-06, pH = 5)
carb_pH8 <- boundary_Cd_NO3(I = 0.001, R = 30E-06 , L = 900E-06, pH = 8)

# Plot of relative change (C/Co) of each chemical species at different pHs.
tiff( "Fig3_D.tiff", res = 72)
oldpar <- par(mfrow=c(2,2)) #, mar=c(4,4,2,1), oma = c(0,0,0,0)) #, mgp=c(1,0,0))
plot(carb_pH7$xlabel*1E+06-30, carb_pH7$Cdenrich, type = "l", lty = 1, xlim = c(0, 100), ylim = c(0.95, 1.05), xlab = xlab, ylab = "relative change", main = expression(paste("D : Cd"^"2+", "  R = 30", mu, "m")))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$Cdenrich, type = "l", lty = 2) # pointillé
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$Cdenrich, type = "l", lty = 3) # points
legend("topright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3)) #, y.intersp = 0.2, bty = "n")

plot(carb_pH7$xlabel*1E+06-30, carb_pH7$CdOHenrich, type = "l", lty = 1, xlim = c(0, 100), ylim = c(0.8, 5), xlab = xlab, ylab = "relative change", main = expression(paste("E : CdOH"^"+", "  R = 30", mu, "m")))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$CdOHenrich, type = "l", lty = 2) # pointillé
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$CdOHenrich, type = "l", lty = 3) # points
legend("topright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3))#, y.intersp = 0.2, bty = "n")

plot(carb_pH7$xlabel*1E+06-30, carb_pH7$CdCO3enrich, type = "l", lty = 1, xlim = c(0, 100), ylim = c(0.8, 5), xlab = xlab, ylab = "relative change", main = expression(paste("F : CdCO"[3], "  R = 30", mu, "m")))
lines(carb_pH5$xlabel*1E+06-30, carb_pH5$CdCO3enrich, type = "l", lty = 2) # pointillé
lines(carb_pH8$xlabel*1E+06-30, carb_pH8$CdCO3enrich, type = "l", lty = 3) # points
legend("topright", legend= c("pH7", "pH5", "pH8"), lty = c(1,2,3))#, y.intersp = 0.2, bty = "n")

par(oldpar)
dev.off()

#-------------------------------------#
#   References                        #
#-------------------------------------#
# Marcus, Y., Volumes of aqueous hydrogen and hydroxide ions at 0 to 200 °C. J. Chem. Phys. 2012, 137, 154501.
# Martell, A. E., Smith, R. M., Motekaitis, R. J. 2004. NIST critical stability constants of metal complexes, version 8. National Institute of Standards and Technology. Gaithersburg, MD. In Gaithersburg, MD, 2004.
#	Schatzberg, P., On the Molecular Diameter of Water from Solubility and Diffusion Measurements. The Journal of Physical Chemistry 1967, 71, 4569-4570.
#	Sethi, M. S.; Raghavan, P. S., Concepts and problems in inorganic chemistry. Discovery Publishing House, 1998; p 425 p.
# Stumm, W., Morgan, J. J. 2006. Aquatic Chemistry Chemical Equilibria and Rates in Natural Waters. Third ed.; United States of America, 1996; p 1022.
# Wilkinson, K. J., Buffle, J. 2004.  Evaluation of Physicochemical Parameters and Processes for Modelling the Biological Uptake of Trace Metals in Environmental (Aquatic) Systems. In Physicochemical Kinetics and Transport at Biointerfaces, Leeuwen, H. P. v.; Köster, W., Eds. John Wiley and Sons Ltd: England, 2004; Vol. IUPAC Series on Analytical and Physical Chemistry of Environmental Systems. Volume 9.
