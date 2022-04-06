#include "constants.h"

subroutine init(u,sizex,sizev,xb,xe,vb,ve,dx,dv)

  implicit none 
  
  integer :: sizex, sizev
  real(kind=DTYPE) :: xb, xe, vb, ve, dx, dv
  real(kind=DTYPE) :: u(-2*B:sizex+2*B, 0:sizev)
  
  integer :: i, j
  real(kind=DTYPE) :: x, y, z
  
  do i=0, sizex
  do j=0, sizev
    
    x = xb + i*(0.5*dx)
    y = vb + j*(0.5*dv)
    
    z = exp(-0.5*(x**2 + y**2))
    !z = sin(x)*cos(y)
    u(i,j) = z
    
  enddo
  enddo

end subroutine init


subroutine init_averaging(u,sizex,sizev,dx,dv)
  
  implicit none 
  
  integer :: sizex,sizev
  real(kind=DTYPE) :: dx,dv
  real(kind=DTYPE) :: u(-2*B:sizex+2*B, 0:sizev)
  
  integer :: i, j
  real(kind=DTYPE) :: avg
  
  do i=1, sizex-1, 2
  do j=1, sizev-1, 2
    
    !naive averaging
    !avg = (1./9.)*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+u(i+1,j+1)+u(i+1,j-1)+u(i-1,j+1)+u(i-1,j-1)+u(i,j))
    
    !Simpson averaging
    avg = (1./36.)*(u(i-1,j-1)+u(i-1,j+1)+u(i+1,j+1)+u(i+1,j-1)) + 4.*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)) + 16.*u(i,j)
    
    u(i,j) = avg
    
  enddo
  enddo
  
end subroutine init_averaging
