#--------------------------------------------------------------------------------------------#
# Two functions calculating membrane permeability (Pm)
# 1) with the Gultknecht and Walter (1986) empirical equation                                #
# 2) with the Xiang and Andersen (1984) equation
#--------------------------------------------------------------------------------------------#

# 1) Gultknecht and Walter (1986) empirical equation  
# The empirical equation of Walter and Gutknecht (1986) relates Pm to Kow (log Pm = s log Kow + b )
# (s = 1.14 ; b = -0.32) for non-electrolytes compounds varying from water (18 g/mol) to codeine (300 g/mol)
# Kow from 4.0 x 10-2 (water) to 1.7 x 10^2 (salicylic acid)

# This function computes Pm and log(Pm) in cm^2/s
# WHERE :
# Kow : The octanol-water partition coefficient (dimensionless) (NOT in log units!)

Pm_Gultkecht <- function(Kow=NULL) {
  if(is.null(Kow))
    stop("Kow must be inputted")
  if ( !is.vector(Kow) )
    stop("Kow should be a vector of length = 1")  # this exclude matrix and arrays
  if ( !length(Kow)==1 )
    stop("Kow should have length = 1")
  if ( !is.numeric(Kow) )
    stop("Kow should be numeric")   # This excludes characters
  if ( any(Kow < 0)  ) 
    stop("Kow should not be < 0")
  s = 1.09
  b = -1.78
  LOGPm <- s * log10(Kow) + b        # LOG10(Pm) where Pm is in cm^2/s
  Pm = 10^(LOGPm)                    # Pm in cm^2/s
  return(list(LOGPm = LOGPm, Pm = Pm))
}

# Example for Maleic hydrazide (logKow = -0.14, Kow = 0.72, Chemspider) (molecular weight = 112.09 g/mol)
#Pm_Gultkecht(Kow = 0.72)
# Example for water (Kow = 4E-02, Walter and Gultnechkt, 1984)
#Pm_Gultkecht(Kow = 4E-02)
# Example for codeine (Kow = 1.6E+01)
#Pm_Gultkecht(Kow = 1.6E+01)
# Example of haxanoic acid (Kow = 7.6E+01)
#Pm_Gultkecht(Kow = 7.6E+01)
# Example for clarithromycin (Kow = 0.72) (but out of the range in MW and Kow consider by Gulterchnekt and Walter)
#Pm_Gultkecht(Kow = 10^3.16)
# Example for 4-OH-CB-107 (Kow = 10^6.91, V = 237)
#Pm_Gultkecht(Kow = 10^6.91)
# Example Fluoxetine
#Pm_Gultkecht(Kow = 10^4.09)
# Example Paroxetine
#Pm_Gultkecht(Kow = 10^3.89)

#----------------------------------------------------------------------------------------------#
# Pm calculations using the Xiang and Anderson (1984) empirical model                          #
#                                                                                              #
#----------------------------------------------------------------------------------------------#
# Xiang and Anderson (1984) tested two models (limiting cases) which ascribe the entire molecular size
# dependence in permeability coefficients to either sizedependent partitioning (n = 0) or size-dependent diffusion
# coefficients (a = 0).
# They used measured permeability coefficients across planar egg lecithin/decane bilayers and bulk hydrocarbon/water
# partition coefficients have been measured for 24 solutes with molecular volumes, V, 
# varying by a factor of 22 and Pm values varying by a factor of 107 to
# explore the chemical nature of the bilayer barrier and the effects of permeant size on permeability

# Although the coefficient of determination was highest when both size-dependent partitioning
# and diffusion were allowed, application of the F-test for equality of variances indicated no significant
# differences in the fits of the data regardless of the form of the size dependence relationship
# used.
# Their models were tested with 22 neutral solutes with a greater Kow and molecular weight range than those used in Gultknecht and Walter
# log Pm delta / Kow = log D0 - n log V -a V / 2.303
# the term " n log V" is the correction for the decrease in Pm due to an increase in diffusion coefficient
# the term " a V / 2.303 " is a correction term for the partitioning of molecules in the bilayer
# Across lecithine bilayer, assuming a = 0 :  log Pm delta / Kow = (- 1.4 log V ) - 3.6
# Across lecithine bilayer, assuming n = 0 :  log Pm delta / Kow = - (0.01 V / 2.303) - 5.9
# delta is assumed to be 18 A or 18 * 0.1 nm = 1.8 nm

# This function computes Pm and log(Pm) in cm/s
# WHERE :
# delta : membrane thickness in A (18 A)
# V : molar volume in cm^3/mol
# Kow : dimentionless (NOT IN LOG UNITS !!!)

Pm_Xiang <- function(Kow=NULL, V=NULL, D=NULL) {
  if(is.null(Kow))
    stop("Kow must be inputted")
  if(is.null(V))
    stop("V must be inputted")
  if(is.null(D))
    stop("D must be inputted")
  if ( !is.vector(Kow) || !is.vector(V) || !is.vector(D) )
    stop("All arguments should be a vector of length = 1")  # this exclude matrix and arrays
  if ( !length(Kow)==1 || !length(V)==1 || !length(D)==1 )
    stop("All arguments should have length = 1")
  if ( !is.numeric(Kow) || !is.numeric(V) || !is.numeric(D) )
    stop("All arguments should be numeric")   # This excludes characters
  if ( any(Kow <= 0) ||  any(V <= 0) || any(D <= 0) ) 
    stop("All arguments should be > 0")
  # Assuming a = 0 (Size-dependant diffusion coefficient is used to corrrect Pm)
  delta = 18 *0.1*1E-07    # 18 A is the assumed membrane thickness
  #LOGPm_a0 = log10(Kow) * ((-1.4 * log10(V)) - 3.6) / delta
  #Pm_a0 = 10^(LOGPm_a0)
  #LOGPm_a0 = (log10(D)  - (1.4 * log10(V)) - 3.6)  * (Kow/delta)
  LOG_x = (log10(D)  - (1.4 * log10(V)) - 3.6)
  Pm_a0 = 10^(LOG_x) * (Kow/delta)
  LOGPm_a0 = log10(Pm_a0)
  # Assuming n = 0 (size-dependant membrane partitioning is used to correct Pm)
  #LOGPm_n0 <- log10(Kow) * (((- 0.01 * V) / 2.303) - 5.9) / delta
  #Pm_n0 <- 10^(LOGPm_n0)
  LOG_y = log10(D) - ((0.01 * V) / 2.303) - 5.9
  Pm_n0 = 10^(LOG_y) * (Kow/delta)
  LOGPm_n0 = log10(Pm_n0)
  return(list(LOGPm_a0 = LOGPm_a0, Pm_a0 = Pm_a0, LOGPm_n0 = LOGPm_n0, Pm_n0 = Pm_n0))
}

# Example for water, but Pm is extremely sensitive to changes in V if V is small (Kow = 4E-02, Chemspider and Gultkecht adn Walter)
#Pm_Xiang(Kow = 4E-02, V = 18 , D = 1.27E-5)

# Example for maleic hydrazide (logKow = -0.14, kow = 0.72 and V = 72.7 for Maleic hydrazide) (Chemspider)
#Pm_Xiang(Kow = 0.72, V = 72.7, D = 7.2E-06)

# Example for clarithromycin (log Kow = 10^3.16) (V = 631.9 cm^3/mol)
#Pm_Xiang(Kow = 10^3.16, V = 631.9, D = 3.5E-06)

# Example for 4-OH-CB-107 (Kow = 10^6.91, V = 237)
#Pm_Xiang(Kow = 10^6.91, V = 237, D = 5.4E-06)

# Example Fluoxetine
#Pm_Xiang(Kow = 10^4.09, V = 266.7, D = 5.18E-06)

# Example Paroxetine
#Pm_Xiang(Kow = 10^3.89, V = 271.5, D = 5.15E-06)

# Example Duloxetine
#Pm_Xiang(Kow = 10^3.73, V = 256.8, D = 5.25E-06)

# Example Clindamycin
#Pm_Xiang(Kow = 10^1.83, V = 327.2, D = 4.84E-06)

# Example of alpha-naphtoic acid (LogKow = 3.13, Kow = 1348.963) (V = 136.1 cm^3/mol) (Chemspider)
#Pm_Xiang(Kow = 1348.963, V = 136.1)


#--------------#
# References   #
#--------------#
# Walter, A., and Gutknecht, J. (1986) Permeability of small nonelectrolytes through lipid bilayer membranes. J Membrane Biology 90: 207-217.
# Xiang, T.-X., and Anderson, B.D. (1994) The relationship between permeant size and permeability in lipid bilayer membranes. J Membrane Biology 140: 111-122.
