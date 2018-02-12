#-----------------------------------------------------------------------------------------------------#
# Correction of equilibrium constant from I > 0 to I=0 (A + B -> C) 
#-----------------------------------------------------------------------------------------------------#
# This uses the Guntelberg equation at I < 0.1 M
# From the concentration-based equilibrium constant, this computes the thermodynamic equilibrium constant (K0)

# Input parameters :
# I is the ionic strength (in mol/L)
# K is the concentration-based equilibrium constant at I > 0 , K is not in log
# chargeA, chargeB or chargeC = 0, 1, 2, 3, or 4
K_corr <- function(K=NULL, I=NULL, chargeA=NULL, chargeB=NULL, chargeC=NULL) {
  if(is.null(K)) 
    stop("K must be inputted")
  if(is.null(I))              
    stop("I must be inputted")
  if(is.null(chargeA))                
    stop("chargeA must be inputted")
  if(is.null(chargeB))                
    stop("chargeB must be inputted")
  if(is.null(chargeC))                
    stop("chargeC must be inputted")
  if( I > 0.1 )
    stop("The Guntelberg equation is valid at I < 0.1 M")
  if ( (!chargeA==0) && (!chargeA==1) && (!chargeA==2) && (!chargeA==3) && (!chargeA==4) ) 
    stop("chargeA must be equal to 0, 1, 2, 3 or 4") 
  if ( (!chargeB==0) && (!chargeB==1) && (!chargeB==2) && (!chargeB==3) && (!chargeB==4) ) 
    stop("chargeB must be equal to 0, 1, 2, 3 or 4") 
  if ( (!chargeC==0) && (!chargeC==1) && (!chargeC==2) && (!chargeC==3) && (!chargeC==4) ) 
    stop("chargeC must be equal to 0, 1, 2, 3 or 4")
  if ( !is.vector(K) || !is.vector(I) || !is.vector(chargeA) || !is.vector(chargeB) || !is.vector(chargeC) )
    stop("all arguments should be vectors of length = 1")  # this exclude matrix and arrays
  if (( !length(K)==1) || (!length(I)==1 ) || (!length(chargeA)==1) || (!length(chargeB)==1) || (!length(chargeC)==1))
    stop("all arguments should have length = 1")
  if ( !is.numeric(K) || (!is.numeric(I)) || (!is.numeric(chargeA)) || (!is.numeric(chargeB)) || (!is.numeric(chargeC)))
    stop("all arguments should be numeric")   # This excludes characters
  if ( any(K < 0) || any(I < 0) || any(chargeA < 0) || any(chargeB < 0) || any(chargeC < 0) ) 
    stop("arguments should not be < 0")
  
  # Temperature = 25 °C and dielectric constant of water
  T <- 25 + 273
  epsilon <- 78.54   # At 25 °C, this is the dielectric constant of water
  A <- 1.82E+06 * (epsilon * T)^(-3/2)
  
  # Guntelbeg equation :
  if(chargeA == 1 || chargeA == 2 || chargeA == 3 || chargeA == 4 )  {
    gammaA = 10^((-A * chargeA^2) * (sqrt(I) / (1 + sqrt(I))))
  }
  if(chargeB == 1 || chargeB == 2 || chargeB == 3 || chargeB == 4 ) {
    gammaB = 10^((-A * chargeB^2) * (sqrt(I) / (1 + sqrt(I))))
  }
  if(chargeC == 1 || chargeC == 2 || chargeC == 3 || chargeC == 4 ) {
    gammaC = 10^((-A * chargeC^2) * (sqrt(I) / (1 + sqrt(I))))
  }
  if(chargeA == 0) {
    gammaA = 1
  }
  if(chargeB == 0) {
    gammaB = 1
  } 
  if(chargeC == 0) {
    gammaC = 1
  }
  
  K0 <- K * ( gammaC / (gammaA * gammaB))
  logK0 <- log10(K0)
  return(list(K0 = K0, logK0 = logK0, gammaA = gammaA, gammaB = gammaB, gammaC = gammaC))
}

# Example of results
#CO2 -> HCO3- + H+ 
#At I = 0 and T = 25 °C
#K1 <- 10^6.35
#K1 <- 10^16.45
#K_corr(K = K1, I = 0.1, chargeA = 2, chargeB = 1, chargeC = 0)

#----------------#
# References     #
#----------------#
# Stumm, W., and Morgan, J.J. (1996) Aquatic Chemistry Chemical Equilibria and Rates in Natural Waters. United States of America.


