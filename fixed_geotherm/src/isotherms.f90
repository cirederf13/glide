subroutine isotherms(matrices,params,nx,ny,elev,lon,lat)

!  This routine calculates the perturbation of isotherms due to topography.

!  This is solved in frequency space so the code proceeds as follows...
!  ... inbed topography in a power of 2 grid
!  ... apply taper to reduce edge effects
!  ... take the discreet fourier transform in the first dimension and then the 2nd dimension

!  ... loop through the number of systems used
!  ... apply continuation function 
!  ... inverse fourier transform (in both dimensions) 

use definitions
use omp_lib

implicit none

type (parm) params
type (matr) matrices

integer i,j,k,dir,size,cnt,sys,nsys
integer nx,ny,nx_2,ny_2,nx1,ny1,ij,ii,jj,i1,j1
integer,dimension(:),allocatable:: systems
double precision u,t,xmin,xmax,ymin,ymax,xlat,xlon
double complex,dimension(:,:), allocatable:: a,a_ffted
double precision lon(nx,ny),lat(nx,ny),elev(nx,ny)
double precision,dimension(:,:,:),allocatable:: s
double precision dx,dy,fi,fj,dist,A_o,xl,yl,xstep,ystep,kdx,kdy,pi
double precision decay,depth,edot,kappa,A_k,manck,pec,factor
        
xmin=minval(lon)
xmax=maxval(lon)
ymin=minval(lat)
ymax=maxval(lat)

xlat=ymin
xlon=xmin

ny_2=2048

print*,'in isotherms'
!minimum power of 2 that's suitable
nx_2 = 2**float(ceiling(log(float(nx*2))/log(2.)))
ny_2 = 2**float(ceiling(log(float(ny*2))/log(2.)))
ny_2=max(nx_2,ny_2)
nx_2=ny_2

!needed to imbed topo
nx1 = (nx_2-nx)/2
ny1 = (ny_2-ny)/2

size=ny_2

allocate(a(size,size),a_ffted(size,size))
 
!imbed topography in a 2**nx_2 x 2**nx_y grid                        
a=0.
do j=1,ny
  do i=1,nx
    ij=(j-1)*nx+i
    ii=nx1+i-1
    jj=ny1+j-1 
    if ((ii-1)*(ii-nx_2).gt.0 .or. (jj-1)*(jj-ny_2).gt.0) then
      print*,ii,jj,i,j,nx1,ny1
      stop 'error in zc calculation, make m in 2**m bigger '
    endif
    a(ii,jj)=cmplx(elev(i,j),0.0)
  enddo 
enddo

!taper topography to reduce edge effects
do i=1,nx1
  a(i,:)=a(nx1,:)*float(i)/float(nx1)
enddo
do i=nx_2-nx1,nx_2
  a(i,:)=a(nx_2-nx1,:)*float(nx_2-i)/float(nx1)
enddo
do j=1,ny1
  a(:,j)=a(:,ny1)*float(j)/float(ny1)
enddo
do j=ny_2-ny1,ny_2
  a(:,j)=a(:,ny_2-ny1)*float(ny1-j)/float(ny1)
enddo

dir=1

!grid spacing in degrees
dx=lon(nx/2,ny/2)-lon(nx/2-1,ny/2)
dy=lat(nx/2,ny/2)-lat(nx/2,ny/2+1)

!size of region
yl=dy*(ny-1)*111.11*1.d3
xl=dx*(nx-1)*111.11*1.d3*cos(((minval(lat)+maxval(lat)/2.)*3.141592654/180.))

!grid spacing meters
xstep=xl/(nx-1)
ystep=yl/(ny-1)

!size of FFT grid
xl=(nx_2-1)*xstep
yl=(ny_2-1)*ystep

!Do first dimension 

! The following line is the OMP directive and is not commented (this is true for all OMP commands)
!$OMP PARALLEL DO SCHEDULE (RUNTIME) 
    do i=1,size
         call four1(a(:,i),size,dir)
        !call four1(a(i,:),size,isign)
    enddo
!$OMP END PARALLEL DO 

!$OMP PARALLEL DO SCHEDULE (RUNTIME) 
    do i=1,size
!                 call four1(a(:,i),size,dir)
        call four1(a(i,:),size,dir)
    enddo
!$OMP END PARALLEL DO

!save the fourier transform of the topography
a_ffted=a

dir=-1

factor=1/(size**2.)
        
!grid spacing in meters
dx=xstep
dy=ystep

!grid spacing in cycles/m
kdx=1./(nx_2*dx)
kdy=1./(ny_2*dy)
pi=atan(1.)*4.

nsys=7
allocate(systems(nsys))

do i=1,params%n
if ((matrices%isys(i)).eq.1) systems(1)=1
if ((matrices%isys(i)).eq.2) systems(2)=1
if ((matrices%isys(i)).eq.3) systems(3)=1
if ((matrices%isys(i)).eq.4) systems(4)=1
if ((matrices%isys(i)).eq.5) systems(5)=1
if ((matrices%isys(i)).eq.6) systems(6)=1
if ((matrices%isys(i)).eq.7) systems(7)=1
enddo

allocate(s(nx,ny,nsys))

edot=params%edot_mean
kappa=35.
!convert to m
edot=edot*1000.
kappa=kappa*(1000.**2.)
pec=edot/(2.*kappa)

! Loops through systems
nsys=7
do sys=1,nsys
!print*,'looping through the systems in isotherms',sys

!use saved transform
a=a_ffted

if (systems(sys).eq.0) then 
   s(:,:,sys)=0. 
   continue
else

depth=matrices%ages(sys,2)*1000.
A_o=(matrices%ages(sys,4)-6)/matrices%ages(sys,5)

!$OMP PARALLEL DO SCHEDULE (RUNTIME) PRIVATE(i,k,fi,fj,decay,A_k,dist,manck)
    do j=1,size
      do i=1,size
        a(i,j)=a(i,j)*factor
        call kvalue(i,j,ny_2,nx_2,kdy,kdx,fj,fi)
        fj=(fj)**2.
        fi=(fi)**2.
        decay=sqrt(fj+fi)      
        !apply filter
        manck=pec+sqrt(((pec)**2.)+ ((2.*pi*decay)**2.))
        if (decay.eq.0.) manck=2.*pec
        A_k=A_o*exp(-1.*(manck*depth))
        a(i,j)=a(i,j)*A_k
      enddo
    enddo
!$OMP END PARALLEL DO 

!Inverse fourier transform

!$OMP PARALLEL DO SCHEDULE (RUNTIME) 
    do i=1,size
         call four1(a(:,i),size,dir)
        !call four1(a(i,:),size,isign)
    enddo
!$OMP END PARALLEL DO 

!$OMP PARALLEL DO SCHEDULE (RUNTIME) 
    do i=1,size
       !  call four1(a(:,i),size,dir)
        call four1(a(i,:),size,dir)
    enddo
!$OMP END PARALLEL DO

!extract filtered topography.            
cnt=1       
do j=1,ny
  do i=1,nx
    ij=(j-1)*nx+i
    ii=nx1+i-1
    jj=ny1+j-1
    if ((ii-1)*(ii-nx_2).gt.0 .or. (jj-1)*(jj-ny_2).gt.0) then
      print*,ii,jj,i,j,nx1,ny1
      stop 'error in zc calculation'
    endif
    s(i,j,sys)=real(a(ii,jj))
  enddo 
enddo

endif

enddo

!Output closure isotherms
!open(25,file=params%run//"/stuff/closure.xyz",status="unknown")
!do j=1,ny
!  do i=1,nx
!    write(25,'(10f12.4)') lon(i,j),lat(i,j),elev(i,j),s(i,j,1),s(i,j,2),s(i,j,3),s(i,j,4),s(i,j,5),s(i,j,6),s(i,j,7)
!  enddo
!enddo
!close(25)

!grid spacing in degrees
dx=lon(nx/2,ny/2)-lon(nx/2-1,ny/2)
dy=lat(nx/2,ny/2)-lat(nx/2,ny/2+1)

!calculate the effective perturbation at the points where there are ages
print*,'calculating perturbation at data points'
do i=1,params%n
 i1=floor((matrices%x(i)-xlon)/dx)
 if (i1.eq.nx) i1=nx-1
 j1=floor((matrices%y(i)-xlat)/dy)
 if (j1.eq.ny) j1=ny-1
 j1=ny-1-j1
! it works on a regular grid. if the grid has been modified, then dx, dy in
! lat/lon varies spacially
 do j=1,nx
   if (lon(j,1).gt.matrices%x(i)) then
     i1=j-1
     if (i1.eq.nx) i1=nx-1
     exit
   endif
  enddo
 do j=1,ny
   if (lat(1,j).gt.matrices%y(i)) then
     j1=j-1
     if (j1.eq.ny) j1=ny-1
     exit
   endif
 enddo
 if (j1.gt.ny) then
 print*,'too big',j1,matrices%y(i),minval(lat),maxval(lat)
 stop
 elseif (j1.le.0) then
 print*,'too small',j1,matrices%y(i),minval(lat),maxval(lat)
 endif
 u = (matrices%x(i)-lon(i1,1))/(lon(i1+1,1)-lon(i1,1))
 t = (matrices%y(i)-lat(1,j1+1))/(lat(1,j1)-lat(1,j1+1))
 if (matrices%isys(i).gt.0) then
 matrices%zc(i)=((1.-t)*(1.-u)*s(i1+1,j1+1,(matrices%isys(i)))) + & 
        (t*(1.-u)*s(i1+1,j1,(matrices%isys(i))))+(t*u*s(i1,j1,(matrices%isys(i))))+&
        ((1.-t)*u*s(i1,j1+1,(matrices%isys(i))))
 matrices%depth1(i)=(1.-t)*(1.-u)*elev(i1+1,j1+1)+ t*(1.-u)*elev(i1+1,j1)+t*u*elev(i1,j1)+&
        (1.-t)*u*elev(i1,j1+1) - matrices%elev(i)
 endif
enddo
print*,'apatite fission track',minval(s(:,:,1)),maxval(s(:,:,3))
print*,'apatite helium',minval(s(:,:,3)),maxval(s(:,:,3))
print*,'data zc',minval(matrices%zc),maxval(matrices%zc)

deallocate(a,a_ffted)
        
return
end subroutine isotherms

