# Polacheck Hilborn and Punt 1993
# Canadian Journal Fisheries and Aquatic Science
# SP Models
# Effort averaging - do not use
# Process error - problematic
# Obs error - use this one
# Here we implement their model

# By+1 = By + g(By) - Cy    (1)
# Iy = q By                 (2)

# Pella-Tomlinson
# g(B) = (r / p) * B * (1 - (B/K)^p)
# p = 1, Schaeffer
# p -> 0, Fox

# Obs error
# Assume pop dynamics (1) is deterministic and all error occurs in relationship
# between stock biomass and index of abundance
# Project forward from Binitial under historic catches.
# Assume error in (2) is multiplicative and lognormal with a constant coefficient
# of variation
# Iy = q By e^E
# E ~ N(0; sigma2)
# Estimates of Binitial, r, q, and K are obtained by maximizing the likelihood
# L = product [exp( - vhaty^2 / (2 sigmahat2) ) / (sqrt(2 pi) sigmahat ]
# product over all years we have data for
# vhat = log(C / E)y - log(C / E)haty
# sigma2hat = sum (vhat 2 / n)
# Ey is fishing effort: q By = Iy = Cy / Ey
# C / E)haty is the model predicted catch rate for year y
# Value of q which minimises L is:
# qhat = exp(1/n sum(log (Iy / By)))

# Also, better to assume that Binitial = K
# set p = 1





