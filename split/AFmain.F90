#include "constants.h"

#include "init.F90"
#include "parameters.F90"
#include "advection2D.F90"
#include "error2D.F90"


program main
  
  use parameters
  
  implicit none 
  
  real(kind=DTYPE), dimension(:,:), allocatable :: fe, fehalf, feold, f0
  real(kind=DTYPE), dimension(:), allocatable :: vx, vy
  
  real(kind=DTYPE) :: t = 0
  real(kind=DTYPE) :: adv_err2D
  integer :: i, ix, iy
  
  allocate(fe    (-2*B:sizex+2*B,-2*B:sizev+2*B))
  allocate(f0    (-2*B:sizex+2*B,-2*B:sizev+2*B))
  allocate(fehalf(-2*B:sizex+2*B,-2*B:sizev+2*B))
  allocate(feold (-2*B:sizex+2*B,-2*B:sizev+2*B))
  
  allocate(vx(-2*B:sizev+2*B))
  allocate(vy(-2*B:sizex+2*B))
  
  call init(fe,sizex,sizev,xb,xe,vb,ve,dx,dv)
  call boundary(fe,sizex,sizev)
  
  call init_averaging2(fe,sizex,sizev,dx,dv)
  call boundary(fe,sizex,sizev)
  
  f0(:,:) = fe(:,:)
  
  !vx
!~   do i=0, sizev
!~   	vx(i) = vb + i * 0.5 * dv
!~   enddo
  
!~   !vy
!~   do i=0, sizex
!~   	vy(i) = 1.*sin(0.2*PI*(xb+i*0.5*dx))
!~   	!vy(i) = 0.
!~   	!vy(i) = sin(PI*(xb+i*0.5*dx))
!~   enddo
  
  vx(:) = ax
  vy(:) = ay
  
  do while(t < tmax)
    
    !stepY
    do ix=0, sizex
      
      call advection_super_simple(fe(ix,:), sizex, vy(ix), dt, dx)
      
    enddo
    call boundary(fe, sizex, sizev)
    
    
    !stepX
    do iy=0, sizev
      
      call advection_super_simple(fe(:,iy), sizev, vx(iy), dt, dv)
      
    enddo
    call boundary(fe, sizex, sizev)
    
    !vy(:) = sin(2.*PI*t) * vy(:)
    
    t = t + dt
    
  enddo
  
  !Output
  do i=-2*B, sizex+2*B 
    open(10, file="f.dat", status ="unknown", position="append")
    write(10,*) fe(i,:)
    close(10)
  end do
  
  do i=-2*B, sizex+2*B 
    open(11, file="f0.dat", status ="unknown", position="append")
    write(11,*) f0(i,:)
    close(11)
  end do
  
  write(*,*) "time : ", t
  write(*,*) "max : ", maxval(fe(:,:))
  write(*,*) "min : ", minval(fe(:,:))
  
  call geterror2D(fe,f0,sizex,sizev,adv_err2D)
  write(*,*) "error2D : "
  write(*,*) dimX, dimV,  adv_err2D
  
  write(*,*) "max : "
  write(*,*) dimX, dimV, maxval(fe(:,:))
  
  write(*,*) "CFL : ", CFL
  write(*,*) "dt : ", dt
  write(*,*) "ax : ", ax
  write(*,*) "ay : ", ay
  
  
  deallocate(fe)
  deallocate(f0)
  deallocate(fehalf)
  deallocate(feold)
  deallocate(vx)
  deallocate(vy)

end program main
