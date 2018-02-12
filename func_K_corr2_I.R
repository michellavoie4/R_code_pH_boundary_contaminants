#---------------------------------------------------------------------------------------------#
# Correction of equilibrium constant from infinite dilution (I=0) to I > 0
#--------------------------------------------------------------------------------------------#
# This function uses the Guntelbrg equation
# This computes the value of a concentration-based equilibrium constant using a thermodynamic equilibrium constant
# A + B -> C
# K0 in 10^x (not in logarithm) : thermodynamic equilibrium constant
# Input parameters :
# I in mol/L (>0)
# charge = 1 to 4 (0 if no charge)
K_corr2 <- function(K0=NULL, I=NULL, chargeA=NULL, chargeB=NULL, chargeC=NULL) {
  if(is.null(K0)) 
    stop("K0 must be inputted")
  if(is.null(I))                
    stop("ion_str must be inputted")
  if( I > 0.1)
    stop("The Guntelberg equation is applicable at I < 0.1 M")
  if(is.null(chargeA))                
    stop("chargeA must be inputted")
  if(is.null(chargeB))                
    stop("chargeB must be inputted")
  if(is.null(chargeC))              
    stop("chargeC must be inputted")
  if ( (!chargeA==0) && (!chargeA==1) && (!chargeA==2) && (!chargeA==3) && (!chargeA==4) ) 
    stop("chargeA must be equal to 0, 1, 2, 3 or 4") 
  if ( (!chargeB==0) && (!chargeB==1) && (!chargeB==2) && (!chargeB==3) && (!chargeB==4) ) 
    stop("chargeB must be equal to 0, 1, 2, 3 or 4") 
  if ( (!chargeC==0) && (!chargeC==1) && (!chargeC==2) && (!chargeC==3) && (!chargeC==4) ) 
    stop("chargeC must be equal to 0, 1, 2, 3 or 4")
  if ( !is.vector(K0) || !is.vector(I) || !is.vector(chargeA) || !is.vector(chargeB) || !is.vector(chargeC) )
    stop("all arguments should be vectors of length = 1")  # this exclude matrix and arrays
  if (( !length(K0)==1) || (!length(I)==1 ) || (!length(chargeA)==1) || (!length(chargeB)==1) || (!length(chargeC)==1))
    stop("all arguments should have length = 1")
  if ( !is.numeric(K0) || (!is.numeric(I)) || (!is.numeric(chargeA)) || (!is.numeric(chargeB)) || (!is.numeric(chargeC)))
    stop("all arguments should be numeric")   # This excludes characters
  if ( any(K0 < 0) || any(I < 0) || any(chargeA < 0) || any(chargeB < 0) || any(chargeC < 0) ) 
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
  
  K1 <- K0 / ( gammaC / (gammaA * gammaB))
  K1_log <- log10(K1)
  return(list(K1=K1, K1_log=K1_log, gammaA = gammaA, gammaB = gammaB, gammaC = gammaC))
}

# example
#K0 <- 4.47 * 10^15
#K_corr2(K0, I=0.1, chargeA=2, chargeB=4, chargeC=0)  
  
#----------------#
# References     #
#----------------#
# Stumm, W., and Morgan, J.J. (1996) Aquatic Chemistry Chemical Equilibria and Rates in Natural Waters. United States of America.

