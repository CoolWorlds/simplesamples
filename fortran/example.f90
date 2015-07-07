! g95 -O3 -c inversebeta.f && g95 -O3 -c coolmath.f90 && g95 -O3 -c simplesamples.f90 && g95 -O3 -o example example.f90 inversebeta.o coolmath.o simplesamples.o

PROGRAM example

use simplesamples ! $compiler -c simplesamples.f90

 REAL(8) :: x

 write(*,*) 'x = ?'
 read(*,*) x
 !write(*,*) 'inverf(x) = ',inverf(x)
 !write(*,*) 'erf[inverf(x)] = ',erf(inverf(x))
 write(*,*) beta_sample(2.0D0,3.0D0,x)

END PROGRAM example
