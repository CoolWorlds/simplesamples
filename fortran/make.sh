#!/bin/bash
#

compiler='gfortran'
rm -f *.o *.mod example *~
$compiler -O3 -c inversebeta.f
$compiler -O3 -c coolmath.f90
$compiler -O3 -c simplesamples.f90
$compiler -O3 -o example example.f90 inversebeta.o coolmath.o simplesamples.o
