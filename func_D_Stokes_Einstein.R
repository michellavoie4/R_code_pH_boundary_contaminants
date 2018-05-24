#---------------------------------------------------------------------------------------------------------------------#
# Calculations of the diffusion coefficient in water using the Stokes-Einstein relation (See Wilke and Chang, 1955)
#---------------------------------------------------------------------------------------------------------------------#

# D coefficient at "infinite dilution" in water
# Assume rigid spherical molecules
# D1 and D2 = Diffusion coefficient (D1 in cm^2/s and D2 in m^2/s)
# T = the absolute temperature (kelvin) (273 + °C)
# n = the solution dynamic viscosity (in centipoise or cp) (cP = 0.89004 at 25 °C) (1 cP = 1 mPa s)
# V = molal volume of solutes at normal boiling point (at 1 atm) (c.c. / mole = cm^3/mol)
# V can be computed with the equation of Schotte (1992), the tables of Wilke and Chang (1955) or the Schroeder's method

D_coeff <- function(T=NULL, n=NULL, V=NULL) {
  if( is.null(T) || is.null(n) || is.null(V) )
    stop("T, n and V must be inputted")
  if ( !is.vector(T) || !is.vector(n) || !is.vector(V) )
    stop("all arguments should be vectors of length = 1")  # this exclude matrix and arrays
  if (( !length(T)==1) || (!length(n)==1 ) || (!length(V)==1) )
    stop("all arguments should have length = 1")
  if ( !is.numeric(T) || (!is.numeric(n)) || (!is.numeric(V)) )
    stop("all arguments should be numeric")   # This excludes characters
  if ( any(T < 0) || any(n < 0) || any(V < 0) ) 
    stop("arguments should not be < 0")
  D1 <- T / (n * 1.004E+07 * (V^(1/3)))  # D1 in cm^2/s
  D2 <- D1 / 10000                       # D2 in m^2/s
  return(list(D1 = D1, D2 = D2))
}

# Example for glucose (D measured = 6.7 x 10-6 cm2/s)
#D_Glu <- D_coeff(T = 298, n = 0.89004, V = 175)
#D_Glu
# Example for water (18 cm^3/mol , Chemspider)
#D_water <- D_coeff(T = 298, n = 0.89004, V = 18)
#D_water

# Fluoxetine (V = 266.7 cm^3/mol, Chemspider)
#D_fluo <- D_coeff(T = 298, n = 0.89004, V = 266.7)
#D_fluo

# Example for paroxetine (V = 271.5 cm^3/mol)
#D_paro <- D_coeff(T = 298, n = 0.89004, V = 271.5)
#D_paro

# Example for duloxetine ( V = 256.8 cm^3/mol)
#D_dulo <- D_coeff(T = 298, n = 0.89004, V = 256.8)
#D_dulo

# Clindamycin ( V = 327.2 cm^3/mol)
#D_clind <- D_coeff(T = 298, n = 0.89004, V = 327.2)
#D_clind

# Example for Maleic hydrazide
#D_Maleic <- D_coeff(T = 298, n = 0.89004, V = 98)
#D_Maleic

# Example for Clarithromycin
#D_Clarith <- D_coeff(T = 298, n = 0.89004, V = 840)
#D_Clarith

# Example for 4-OH-CB-107 (V = ????)
#D_PCB <- D_coeff(T = 298, n = 0.89004, V = 236.8)
#D_PCB

#------------#
# References #
#------------#
# Wilke, C.R., Chang, P. 1955. Correlation of diffusion coefficient in dilute solutions. A.I.Ch.E. Journal. 264-270
# Schotte W. 1992. Prediction of the molar volume at the normal boiling point. The Chemical Engineering Journal, 48 (1992) 167-172
# Baum, E.J. 1997. Chemical Property Estimation: Theory and Application. Gopgle book.
