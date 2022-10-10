#include "constants.h"

subroutine boundary(u,sizex,sizev)
  !periodic boundary condition
  implicit none

  integer :: sizex, sizev
  real(kind=DTYPE) :: u(-2*B:sizex+2*B,-2*B:sizev+2*B)

  integer :: ix,iy

  !x-direction
  do ix=1,2*B
    u(-ix,:) = u(sizex-ix,:)
    u(sizex+ix,:) = u(ix,:)
  enddo
  
  !y-direction
  do iy=1,2*B
    u(:,-iy) = u(:,sizev-iy)
    u(:,sizev+iy) = u(:,iy)
  enddo
  
end subroutine boundary

subroutine advection_super_simple(u,sizex,speed,dt,dx)
		  implicit none
		  
		  integer :: sizex
		  real(kind=DTYPE) :: speed,dt,dx
		  real(kind=DTYPE) :: u(-2*B:sizex+2*B)
		  
		  integer :: i
		  real(kind=DTYPE) :: nu
		  real(kind=DTYPE), dimension(-2*B:sizex+2*B) :: uold
		  
		  uold(:) = u(:)
		  !uhalf(:) = u(:)
		  
		  !evolution
		  nu = dt*abs(speed)/dx !CFL
		  
		  !positive advection velocity
		  if (speed >= 0.) then 
		  	
		  	!interface update formula
		  	do i=0, sizex, 2
		  	
		  	  u(i) = nu * (3.*nu - 2.) * uold(i-2) &
		  	       + 6.*nu*(1. - nu) * uold(i-1) &
		  	       + (1.-nu)*(1. - 3.*nu) * uold(i)
		  	
		  	enddo
		  	
		  	!average update formula
		  	do i=1, sizex-1, 2
		  	
		  	  u(i) = nu**2 * (nu-1.) * uold(i-3) &
		  	       + nu**2 * (3.-2.*nu) * uold(i-2) & 
		  	       + nu * (1.-nu) * uold(i-1) &
		  	       + (1.-nu)**2 * (1.+2.*nu) * uold(i) &
		  	       - nu*(1.-nu)**2 * uold(i+1) 
		  	    
		  	enddo
		  	
		  !negative advection velocity
		  else if (speed < 0.) then 
		  
		  do i=0, sizex, 2
		  
		    u(i) = (1.-nu) * (1. - 3.*nu) * uold(i) &
		         + 6.*nu*(1. - nu) * uold(i+1) &
		         + nu * (3.*nu - 2.) * uold(i+2)
		  
		  enddo
		  
		  do i=1, sizex-1, 2
		  
		    u(i) = nu**2 * (nu-1.) * uold(i+3) &
		  	     + nu**2 * (3.-2.*nu) * uold(i+2) & 
		  	     + nu * (1.-nu) * uold(i+1) &
		  	     + (1.-nu)**2 * (1.+2.*nu) * uold(i) &
		  	     - nu*(1.-nu)**2 * uold(i-1) 
		  
		  enddo
		  
		  endif
		  	
end subroutine advection_super_simple
