#-----------------------------------------------------------------#
# Function yielding Kos following Bjerrum's ion association model #
#-----------------------------------------------------------------#

# This equation is described in Kowalak et al (1967) and Wilkins (1970)
# Valid at I < 5 mM since it assumes acitivity coefficients = 1
# Input : Z1: Absolute electric charge of the metal
# Z2 : Absolute electric charge of the ligand
# a : minimum approach distance between the charged ions (in cm )  (= 1E-10* (ion_rad + water_diam + ligand_rad)), where on_rad and water_rad and ligand_rad are in picometer
# T : temperature in kelvins (Temperature in Kelvins = 273 + temeprature in celsius)
fun_Kos <- function(Z1=NULL, Z2=NULL, a=NULL, T=NULL) {
  if ( is.null(Z1) || is.null(Z2) || is.null(a) || is.null(T) )
    stop("Z1, Z2, a and T must be inputted")
  if ( !is.vector(Z1) || !is.vector(Z2) || !is.vector(a) || !is.vector(T) )
    stop("all arguments should be vectors of length = 1")  # this exclude matrix and arrays
  if (( !length(Z1)==1) || (!length(Z2)==1 ) || (!length(a)==1) || (!length(T)==1) )
    stop("all arguments should have length = 1")
  if ( !is.numeric(Z1) || (!is.numeric(Z2)) || (!is.numeric(a)) || (!is.numeric(T)) )
    stop("all arguments should be numeric")   # This excludes characters
  if ( any(Z1 < 0) || any(Z2 < 0) || any(a < 0) || any(T < 0) ) 
    stop("arguments should not be < 0")
  
  k <- 1.38064852E-23    # Boltzman's constant (m^2 Kg s-2 K-1)
  D <- 78                # Dielectric constant
  kT <- k * T  * 10^7      # Work in Erg or g cm^2 s-1 (1 Erg = 10-7)
  gamma_plus <- 1          # activity coefficient (1 is valid at I < 5 mM)
  e <- 3.7673E-10         # Elementary charge in the electrostatic system of units within the CGS system (in StatCoulomb)
  N <- 6.02E+23           # Avogadro's number (mol-1)
  Kos <- ((4 * pi * N * (a^3)) / 3000 ) * exp(((Z1 * Z2) * (e^2))/ (a * D * kT)) * gamma_plus^2
  log_Kos <- log10(Kos)
  return(list(Kos=Kos, log_Kos=log_Kos))
}

# Example
#fun_Kos(Z1=0, Z2=2, a=300E-10, T=298)
#----------------#
# References     #
#----------------#
# Kowalak, A., Kustin, K., Pasterna, R.F., and Petrucci, S. (1967) Steric effects in fast metal complex substitution reactions. 2. Journal of the American Chemical Society 89: 3126-&.
# Wilkins, R.G. (1970) Mechanisms of ligand replacement in octahedral nickel(II) complexes. Accounts of Chemical Research 3: 408-&.
