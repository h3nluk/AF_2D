#include "constants.h"


subroutine geterror2D(f, f0, sizex, sizev, err2D)

  implicit none
  
  integer :: sizex, sizev
  real(kind=DTYPE), dimension(-2*B:sizex+2*B, -2*B:sizev+2*B) :: f, f0
  real(kind=DTYPE) :: err2D
  
  integer :: i, j
  real(kind=DTYPE) :: denom, numer
  
  numer = 0.
  denom = 0.
  
  do i=0, sizex
  do j=0, sizev
    
    numer = numer + abs(f0(i,j) - f(i,j))
    denom = denom + abs(f0(i,j))
    
  enddo
  enddo
  
  err2D = numer/denom
  
end subroutine geterror2D
