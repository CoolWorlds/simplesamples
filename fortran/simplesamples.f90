MODULE simplesamples

use coolmath		! $compiler -c coolmath.f90
use inversebeta		! $compiler -c inversebeta.f90

 implicit none

 CONTAINS

 ! =============================================================================
 REAL(8) FUNCTION uniform_sample(xmin,xmax,x)

  implicit none

  REAL(8), INTENT(IN) :: xmin, xmax
  REAL(8), INTENT(IN) :: x

  uniform_sample = xmin + (xmax-xmin)*x

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION gaussian_sample(mu,sigma,x)

  implicit none

  REAL(8), INTENT(IN) :: mu, sigma
  REAL(8), INTENT(IN) :: x
  REAL(8), PARAMETER :: roottwo = 1.414213562373095D0

  gaussian_sample = mu + roottwo*sigma*inverf(-1.0D0+2.0D0*x)

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION cappedgaussian_sample(xmin,mu,sigma,x)

  implicit none

  REAL(8), INTENT(IN) :: xmin, mu, sigma
  REAL(8), INTENT(IN) :: x
  REAL(8) :: temp
  REAL(8), PARAMETER :: roottwo = 1.414213562373095D0

  temp = DERFC( (-xmin+mu)/(roottwo*sigma) ) + &
         x*( 2.0D0 - DERFC( (-xmin+mu)/(roottwo*sigma) ) )
  cappedgaussian_sample = mu + roottwo*sigma*inverf(-1.0D0+temp)

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION truncatedgaussian_sample(xmin,xmax,mu,sigma,x)

  implicit none

  REAL(8), INTENT(IN) :: xmin, xmax, mu, sigma
  REAL(8), INTENT(IN) :: x
  REAL(8) :: temp
  REAL(8), PARAMETER :: roottwo = 1.414213562373095D0

  ! Truncated Gaussian: Rsol=mean, Rdel=stdev, Rmin=min, Rmax=max
  temp = DERFC( (-xmin+mu)/(roottwo*sigma) ) + &
         x*( DERFC( (-xmax+mu)/(roottwo*sigma) ) - &
         DERFC( (-xmin+mu)/(roottwo*sigma) ) )
  truncatedgaussian_sample = mu + roottwo*sigma*inverf(-1.0D0+temp)

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION jeffreys_sample(xmin,xmax,x)

  implicit none

  REAL(8), INTENT(IN) :: xmin, xmax
  REAL(8), INTENT(IN) :: x

  jeffreys_sample = ( xmax**x )*( xmin**(1.0D0-x) )

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION modjeffreys_sample(xinflection,xmax,x)

  implicit none

  REAL(8), INTENT(IN) :: xinflection, xmax
  REAL(8), INTENT(IN) :: x

  modjeffreys_sample = -( xinflection**(1.0D0-x) )*( xinflection**x &
                       - ( xinflection+xmax )**x )

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION beta_sample(alpha,beta,x)

  implicit none

  REAL(8), INTENT(IN) :: alpha, beta
  REAL(8), INTENT(IN) :: x
  REAL(8) :: beta_log
  INTEGER :: ifault

  beta_log = alngam(alpha,ifault) + alngam(beta,ifault) &
             - alngam(alpha+beta,ifault)

  beta_sample = xinbta(alpha,beta,beta_log,x,ifault )

 END FUNCTION
 ! =============================================================================

END MODULE simplesamples
