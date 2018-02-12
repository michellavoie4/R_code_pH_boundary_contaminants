#---------------------------------------------------------------------------------------------------------------------#
# Calculations of second order rate constant of a diffusion-limited reaction using the Smoluchowski-Debye equation
#---------------------------------------------------------------------------------------------------------------------#

# This function calculates the second order rate constant of a diffusion-limited reaction as a function of T as well as the ionic charge, the diffusion coefficients and the radii of the reactants
# When both reactants are charged, the long-range Coulombic interactions that affect diffusion is calculated and the Smoluchowski equation (valid for two neutral ions) is corrected for the presence of two ions. 
# Da and Db = Diffusion coefficient of reactant a and reactant b (cm^2/s)
# rab = distance of the closest approach at which the proton transfer occurs (7.5 A, Eigen and Maeyer, 1958; Laidler, 1964)
# epsi = relative permitivitty of the solvent or dielectric constant (for water at 20 C, epsi = 80.1) (for water at 25 C = 78.54)
# e = 3.7673E-10   Elementary charge in the electrostatic system of units within the CGS system (in StatCoulomb)
# k = 1.38064852E-16    # Boltzman's constant (cm^2 g s-2 K-1)
# kT = k * T      # Work
# Za = Electric charge of A (with the positive or the negative sign)
# Zb = Electric charge of B (with the positive or the negative sign)
# The function compute the maximum diffusion rate (k_diff) with the Debye correction (k_diff) and without the electrostatic correction (kd, for neutral reactants as well as the factor (proposed by Debye) (ratio) enhancing or decreasing Kd)

kdiff <- function(T=NULL, Da=NULL, Db=NULL, Za=NULL, Zb=NULL, epsi=NULL) {
  if(is.null(T))
    stop("T must be inputted")
  if(is.null(Da))               
    stop("Da must be inputted")
  if(is.null(Db))             
    stop("Db must be inputted")
  if(is.null(Za))                
    stop("Za must be inputted")
  if(is.null(Zb))               
    stop("Zb must be inputted")
  if(is.null(epsi))           
    stop("epsi must be inputted")
  if ( !is.vector(T) || !is.vector(Da) || !is.vector(Db) || !is.vector(Za) || !is.vector(Zb) || !is.vector(epsi) )
    stop("all arguments should be vectors of length = 1")  # this exclude matrix and arrays
  if ( (!length(T)==1) || (!length(Da)==1 ) || (!length(Db)==1) || (!length(Za)==1) || (!length(Zb)==1) || (!length(epsi)==1) )
    stop("all arguments should have length = 1")
  if ( !is.numeric(T) || (!is.numeric(Da)) || (!is.numeric(Db)) || (!is.numeric(Za)) || (!is.numeric(Zb)) || (!is.numeric(epsi)) )
    stop("all arguments should be numeric")   # This excludes characters
  if ( any(T <= 0) || any(Da <= 0) || any(Db <= 0) || any(epsi <= 0) ) 
    stop("T, Da, Db and epsi should be > 0")
  e = 4.803E-10 #4.7673E-10 # 4.803E-10 #3.7673E-10   # Elementary charge in the electrostatic system of units within the CGS system (in StatCoulomb)
  k = 1.38064852E-16    # Boltzman's constant (cm^2 g s-2 K-1)
  kT = k * T    # Work in g cm^2 s-1 
  rab = 750/(1E+10)   #750                     # For most ionic recombinations, the closest approach distance orr the average separation at which the proton transfer occur is often 7.5 A (Laidler, 1965)
  rab_A = (110+10)/(1E+10)              # in cm, since 1 A = 0.1 nm ou 100 pm ; 1 cm = 1E+10 pm, radius of OH-   + diameter of a water molecule  + radius H+
  N = 6.023E+23                           # Avogadro's number
  F = (Za * Zb * (e^2) ) / (epsi * kT * rab)
  kdiff = ( 4 * pi * N * (e^2) * Za * Zb * (Da + Db) ) / ( (1000 * epsi * kT) * (exp(F) - 1) ) # second order rate constant for a diffusion-controlled reaction involving two neutral reactants
  ratio = ( F / (exp(F) - 1))                               
  Rc = e^2 / (4 * pi * epsi * kT)     # Osanger distance
  kd = (Da + Db) * rab * ((4 * pi * N )/1000)
  return(list(kdiff = kdiff, ratio = ratio, kd = kd))
}

# Example for the couple H3O+ + OH- -> 2 H2O (measured rate constant = 1.4 x 10^11 L mol-1 s-1 at 25 °C, Laidler, 1965)
#kdiff_H2O <- kdiff(T = 298, Da = 9.31E-05, Db = 5.27E-05, Za = 1, Zb = -1, epsi = 78.54)  # 80.1
#kdiff_H2O
#log10(kdiff_H2O$kdiff)
#log10(kdiff_H2O$kd)
#log10(1.4E+11)
#erreur <- 100 * ((abs(kdiff_H2O$kdiff - 1.4E+11))/1.4E+11 )    # Error in percent
#erreur

# Example for the couple H3O+ + CH3COO- -> H2O + CH3COOH (measured rate constant = 4.5 x 10^10 L mol-1 s-1 or , Laidler, 1965)
#source("~/R software/R scripts/Functions/func_D_Stokes_Einstein.R")
#D_acet <- D_coeff(T = 298, n = 0.89004, V = 56)$D1
#D_acet
#kdiff_acet <- kdiff(T = 298, Da = 9.31E-05, Db = D_acet, Za = 1, Zb = -1, epsi = 78.54)
#kdiff_acet
#log10(kdiff_acet$kdiff)
#log10(kdiff_acet$kd)
#log10(4.5E+10)
#erreur <- 100*((abs(kdiff_acet$kdiff - 4.5E+10))/4.5E+10 )    # Error in percent
#erreur

# Example for the couple H3O+ + imidazole0 -> H2O + imidazolium+ ion (measured rate constant = 1.5 x 10^10 L mol-1 s-1, Laidler, 1965)
#D_imi <- D_coeff(T = 298, n = 0.89004, V = 70)$D1
#D_imi
#kdiff_imi <- kdiff(T = 298, Da = 9.31E-05, Db = D_acet, Za = 1, Zb = 0, epsi = 78.54)
#kdiff_imi
#log10(kdiff_imi$kdiff)
#log10(kdiff_imi$kd)
#log10(1.5E+10)
#erreur <- 100*((abs(kdiff_acet$kdiff - 1.5E+10))/1.5E+10 )    # Error in percent
#erreur

# Example for Clindamycin_neutral + H+ <-> Clindamycin_H_Plus  (kf = 10^10.75)
#D_Clind <- D_coeff(T = 298, n = 0.89004, V = 327.2)$D1
#D_Clind
#kdiff_Clind <- kdiff(T = 298, Da = 9.31E-05, Db = D_Clind, Za = 1, Zb = 0, epsi = 78.54)
#kdiff_Clind
#log10(kdiff_Clind$kdiff)
#log10(kdiff_Clind$kd)

# Example for Maleic hydrazide : Mal- + H+ <-> MalH (kf = 10^10.53 M-1 s-1)
#D_Mal <- D_coeff(T = 298, n = 0.89004, V = 98)$D1     
#D_Mal <- D_Mal
#kdiff_Mal <- kdiff(T = 298, Da = 9.31E-05, Db = D_Mal, Za = 1, Zb = 1, epsi = 78.54)
#kdiff_Mal
#log10(kdiff_Mal$kdiff)
#log10(kdiff_Mal$kd)

# Example for Clarithromycin : Clarith_neutral + H+ -> Clarith_H+  (kf = 10^10.74 M-1 s-1)
#D_Cla <- D_coeff(T = 298, n = 0.89004, V = 840)$D1    
#D_Cla <- D_Cla
#kdiff_Cla <- kdiff(T = 298, Da = 9.31E-05, Db = D_Cla, Za = 1, Zb = 0, epsi = 78.54)
#kdiff_Cla
#log10(kdiff_Cla$kdiff)
#log10(kdiff_Cla$kd)

# Example of H+ + CO32- -> HCO3-
# D of CO32- was calculated by Li et al 1974 and was 0.955 x 10-5 cm^2/s
#D_CO3 <- D_coeff(T = 298, n = 0.89004, V = 35)$D1
#D_CO3
#kdiff_CO3 <- kdiff(T = 298, Da = 9.31E-05, Db = D_CO3, Za = 1, Zb = -2, epsi = 78.54)
#kdiff_CO3
#log10(kdiff_CO3$kdiff)
#log10(kdiff_CO3$kd)

# Example of H3O+ + F- -> H2O + HF (measured rate cst = 1 x 10^11 L mol-2 s-1) (Laidler, 1965)
#D_F <- D_coeff(T = 298, n = 0.89004, V = 10.5)$D1
#D_F
#kdiff_D_F <- kdiff(T = 298, Da = 9.31E-05, Db = D_F, Za = 1, Zb = -1, epsi = 78.54)
#kdiff_D_F
#log10(kdiff_D_F$kdiff)
#log10(kdiff_D_F$kd)
#erreur <- 100*((abs(kdiff_D_F$kdiff - 1E+11))/1E+11 )    # Error in percent
#erreur

# Example for the reaction H3O+ + SO42- -> H2O + HSO4-  (Measured rate cst = 1 x 10^11 L mol-1 s-1) (Laidler, 1965)
#D_S <- D_coeff(T = 298, n = 0.89004, V = 63)$D1
#D_S
#kdiff_S <- kdiff(T = 298, Da = 9.31E-05, Db = D_S, Za = 1, Zb = -1, epsi = 78.54)
#kdiff_S
#log10(kdiff_S$kdiff)
#log10(kdiff_S$kd)
#erreur <- 100*((abs(kdiff_S$kdiff - 1E+11))/1E+11 )    # Error in percent
#erreur


#---------------#
# References    #
#---------------#
# Eigen M. and L. De Maeyer. 1958. Proc. Roy Soc. (London), A247, 505
# Eigen, M., and Eyring, E. M. 1962. J. Am. Chem. Soc., 84, 3254
# Eigen, M. 1960. Z. Elektrochem., 64. 1 15
# Laidler, K.J. 1964. Chemical kinetics. New York: McGraw-Hill Book Company. 566 p.
