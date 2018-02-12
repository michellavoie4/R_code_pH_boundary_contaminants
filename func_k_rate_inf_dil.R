#-----------------------------------------------------------------------------#
# Correction of rate constants from I > 0 TO infinite dilution (I=0)          #
#-----------------------------------------------------------------------------#
# This function uses the Guntelberg equation to determine activity coefficients

# A + B -> C
# k in 10^x (not in logarithm) : rate constants at I > 0
# I in mol/L (>0)
# chargeA, chargeB or chargeC = 1 to 4 (0 if no charges) (absolute values, no + or -)
# charge C is the charge (absolute value) of the activated complex (charge of A + charge of B = charge C). If A2+ + B2- -> C, chargeC = 0, not 4
# This computes k0 (rate constants at I = 0) and three activity coefficients (gammaA, gammaB and gammaC) from k at I >0
# Increasing the ionic strength lowers the reaction rate between a cation and an anion; increases the reaction rate between like-charged species; has little or no effect on reaction rate when one or both of the reactants is uncharged
k_corr <- function(k=NULL, I=NULL, chargeA=NULL, chargeB=NULL, chargeC=NULL) {
  if(is.null(k))
    stop("k must be inputted")
  if(is.null(I))                
    stop("ion_str must be inputted")
  if(I > 0.1)
    stope("The Guntelberg equation is applicable at I < 0.1 M")
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
  if ( !is.vector(k) || !is.vector(I) || !is.vector(chargeA) || !is.vector(chargeB) || !is.vector(chargeC) )
    stop("all arguments should be vectors of length = 1")  # this exclude matrix and arrays
  if (( !length(k)==1) || (!length(I)==1 ) || (!length(chargeA)==1) || (!length(chargeB)==1) || (!length(chargeC)==1))
    stop("all arguments should have length = 1")
  if ( !is.numeric(k) || (!is.numeric(I)) || (!is.numeric(chargeA)) || (!is.numeric(chargeB)) || (!is.numeric(chargeC)))
    stop("all arguments should be numeric")   # This excludes characters
  if ( any(k < 0) || any(I < 0) || any(chargeA < 0) || any(chargeB < 0) || any(chargeC < 0) ) 
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
  
  k0 <- k * ( gammaC / (gammaA * gammaB))
  k0_log <- log10(k0)
  return(list(k0=k0, k0_log=k0_log, gammaA = gammaA, gammaB = gammaB, gammaC = gammaC))
}

# Example
#ex <- k_corr(k=6E+09, I=0.1, chargeA=1, chargeB=1, chargeC=2)
#ex
#----------------#
# References     #
#----------------#
# Stumm, W., and Morgan, J.J. (1996) Aquatic Chemistry Chemical Equilibria and Rates in Natural Waters. United States of America.
