#include "constants.h"

#include "init.F90"
#include "parameters.F90"
#include "advection2D.F90"

program main
  
  use parameters
  
  implicit none 
  
  real(kind=DTYPE), dimension(:,:), allocatable :: fe, fehalf, feold
  real(kind=DTYPE), dimension(:), allocatable :: v
  real(kind=DTYPE), dimension(:), allocatable :: E
  
  real(kind=DTYPE) :: t = 0
  integer :: i
  
  allocate(fe(-2*B:sizex+2*B,0:sizev))
  allocate(fehalf(-2*B:sizex+2*B,0:sizev))
  allocate(feold(-2*B:sizex+2*B,0:sizev))
  allocate(v(0:sizev))
  allocate(E(-2*B:sizex+2*B))
  
  call init(fe,sizex,sizev,xb,xe,vb,ve,dx,dv)
  call boundary(fe,sizex,sizev)
  
  !advection velocities
  do i=0, sizev
  	v(i) = vb + i * 0.5 * dv
  enddo
  
  do i=-2*B, sizex+2*B
  	E(i) = 0.35*sin((xb+i*0.5*dx)*0.25)
  	!E(i) = 0.
  enddo
  
  !E(:) = 0.
  
  do while(t < tmax)
    
    feold(:,:) = fe(:,:)
    fehalf(:,:) = fe(:,:)
    
    !Interface-Update
    call Evolution(fe,sizex,sizev,v,E,dx,dv,dt)
    call boundary(fe,sizex,sizev)
    call Evolution(fehalf,sizex,sizev,v,E,dx,dv,0.5*dt)
    call boundary(fehalf,sizex,sizev)
    
    !Average-Update
    call Conservation(fe,fehalf,feold,v,E,E,E,dt,dx,dv,sizex,sizev)
    
    t = t + dt 
  enddo
  
  !Output
  do i=-2*B, sizex+2*B 
    open(10, file="f.dat", status ="unknown", position="append")
    write(10,*) fe(i,:)
    close(10)
  end do
  
  write(*,*) "time: ", t
  write(*,*) "max: ", maxval(fe(:,:))
  write(*,*) "min: ", minval(fe(:,:))
  
  
  deallocate(fe)
  deallocate(fehalf)
  deallocate(feold)
  deallocate(v)
  deallocate(E)

end program main
