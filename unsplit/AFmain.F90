#include "constants.h"

#include "init.F90"
#include "parameters.F90"
#include "advection2D.F90"
#include "error2D.F90"


program main
  
  use parameters
  
  implicit none 
  
  real(kind=DTYPE), dimension(:,:), allocatable :: fe, fehalf, feold, f0
  real(kind=DTYPE), dimension(:), allocatable :: v
  real(kind=DTYPE), dimension(:), allocatable :: E
  
  real(kind=DTYPE) :: t = 0
  real(kind=DTYPE) :: adv_err2D
  integer :: i
  
  allocate(fe    (-2*B:sizex+2*B,-2*B:sizev+2*B))
  allocate(f0    (-2*B:sizex+2*B,-2*B:sizev+2*B))
  allocate(fehalf(-2*B:sizex+2*B,-2*B:sizev+2*B))
  allocate(feold (-2*B:sizex+2*B,-2*B:sizev+2*B))
  
  allocate(v(-2*B:sizev+2*B))
  allocate(E(-2*B:sizex+2*B))
  
  call init(fe,sizex,sizev,xb,xe,vb,ve,dx,dv)
  call boundary(fe,sizex,sizev)
  
  call init_averaging(fe,sizex,sizev,dx,dv)
  
  f0(:,:) = fe(:,:)
  
  !advection velocities
!~   do i=0, sizev
!~   	v(i) = vb + i * 0.5 * dv
!~   enddo
  
!~   do i=-2*B, sizex+2*B
!~   	E(i) = 0.35*sin((xb+i*0.5*dx)*0.25)
!~   	!E(i) = 0.
!~   enddo
  
  v(:) = ax
  E(:) = ay
  
  do while(t < tmax)
    
    feold(:,:)  = fe(:,:)
    fehalf(:,:) = fe(:,:)
    
    !Interface-Update
    call Evolution(fe,sizex,sizev,v,E,dx,dv,dt)
    call boundary(fe,sizex,sizev)
    
    call Evolution(fehalf,sizex,sizev,v,E,dx,dv,0.5*dt)
    call boundary(fehalf,sizex,sizev)
    
    !Average-Update
    call Conservation(fe,fehalf,feold,v,E,E,E,dt,dx,dv,sizex,sizev)
    call boundary(fe,sizex,sizev)
    
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
  deallocate(v)
  deallocate(E)

end program main
