#!/bin/bash

#cpp vlasov_poisson.F90 > vlasov_poisson_prep.F90
#f2py -c vlasov_poisson_prep.F90 -m vp

#cpp vlasov_advection.F90 > vlasov_advection_prep.F90
#f2py -c vlasov_advection_prep.F90 -m af

#cpp vlasov_poisson.F90 > vlasov_poisson_prep.F90
#gfortran -std=f2008 vlasov_poisson_prep.F90 -o vlasov-poisson

#cpp vlasov_advection.F90 > vlasov_advection_prep.F90
#gfortran -std=f2008 vlasov_advection_prep.F90 -o vlasov-advection


#clearing
if [ -e a.out ]
then 
  rm a.out
fi

if [ -e Emax.dat ]
then
  rm Emax.dat
fi

if [ -e E.dat ]
then
  rm E.dat
fi

if [ -e f.dat ]
then
  rm f.dat
fi

if [ -e f0.dat ]
then
  rm f0.dat
fi

if [ -e Esq.dat ]
then
  rm Esq.dat
fi  

if [ -e Phi.dat ]
then
  rm Phi.dat
fi 

if [ -e Rho.dat ]
then 
  rm Rho.dat
fi

#compiling
gfortran -ffree-line-length-512 AFmain.F90

#executing
./a.out

#plotting
python3 plot.py

