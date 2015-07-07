import math
import coolmath
import scipy.special

# ==============================================================================
def uniform_sample(xmin,xmax,x):
  y = xmin + (xmax-xmin)*x
  return y
# ==============================================================================

# ==============================================================================
def gaussian_sample(mu,sigma,x):
  roottwo = 1.414213562373095
  y = mu + roottwo*sigma*coolmath.inverf(-1.0+2.0*x)
  return y
# ==============================================================================

# ==============================================================================
def cappedgaussian_sample(xmin,mu,sigma,x):
  roottwo = 1.414213562373095
  temp = math.erfc( (-xmin+mu)/(roottwo*sigma) ) + \
         x*( 2.0 - math.erfc( (-xmin+mu)/(roottwo*sigma) ) )
  y = mu + roottwo*sigma*coolmath.inverf(-1.0+temp)
  return y
# ==============================================================================

# ==============================================================================
def truncatedgaussian_sample(xmin,xmax,mu,sigma,x):
  roottwo = 1.414213562373095
  temp = math.erfc( (-xmin+mu)/(roottwo*sigma) ) + \
         x*( math.erfc( (-xmax+mu)/(roottwo*sigma) ) - \
         math.erfc( (-xmin+mu)/(roottwo*sigma) ) )
  y = mu + roottwo*sigma*coolmath.inverf(-1.0+temp)
  return y
# ==============================================================================

# ==============================================================================
def modjeffreys_sample(xinflection,xmax,x):
  y = -( xinflection**(1.0-x) )*( xinflection**x - ( xinflection+xmax )**x )
  return y
# ==============================================================================

# ==============================================================================
def beta_sample(alpha,beta,x):
  y=scipy.special.betaincinv(alpha, beta, x)
  return y
# ==============================================================================
