#---------------------------------------------------------------------------------#
# Calculations of the activity coefficient as a function of the ionic strength    #
#---------------------------------------------------------------------------------#

# This function calculates the activity coefficient as a function of I
# Where charge is the absolute charge of the ion (always positive)
# Where I is the ionic strength (I in mol/L)

coeff_act <- function(charge = NULL, I = NULL) {
  if ( is.null(charge) ) 
    stop("charge must be inputted")
  if ( !is.null(charge) ){
  if ( (!charge==0) && (!charge==1) && (!charge==2) && (!charge==3) && (!charge==4) ) 
    stop("charge must be equal to 0, 1, 2, 3 or 4") 
    }
  if ( is.null(I) )              
    stop("ionic strength must be inputted")
  if ( I > 0.1 )
    stop("The Guntelberg equation is applicable at I < 0.1 M")
  if ( !is.vector(charge) || (!is.vector(I)))
    stop("charge and ionic strength should be vectors of length = 1")  # this exclude matrix and arrays
  if (( !length(charge)==1) || (!length(I)==1 ))
    stop("charge and ionic strength should have length = 1")
  if (!is.numeric(charge) || (!is.numeric(I)))
    stop("charge and ionic strength should be numeric")   # This excludes characters
  if (any(charge < 0) || any(I < 0)) 
    stop("charge and ion_str should not be < 0")
 
  # Temperature = 25 °C and dielectric constant of water
  T <- 25 + 273
  epsilon <- 78.54   # At 25 °C, this is the dielectric constant of water
  A <- 1.82E+06 * (epsilon * T)^(-3/2)
  
  # Guntelbeg equation :
  if(charge == 1 || charge == 2 || charge == 3 || charge == 4 )  {
    gamma = 10^((-A * charge^2) * (sqrt(I) / (1 + sqrt(I))))
  }
  if(charge == 0) {
    gamma = 1
  }
  coeff_act <- gamma
  return(list(coeff_act = coeff_act))
  
}

# Example
#coeff_act(charge = 4, I = 0.001)

#----------------#
# References     #
#----------------#
# Stumm, W., and Morgan, J.J. (1996) Aquatic Chemistry Chemical Equilibria and Rates in Natural Waters. United States of America.

