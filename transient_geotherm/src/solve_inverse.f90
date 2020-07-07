subroutine solve_inverse(matrices,params)

! In this subroutine the maximum likelihood estimates of the exhumation rates are obtained.

! edot = edot_pr + CG'(GCG'+Cee) (Zc - (Aedot_pr))

! First ACA'+Cee is computed adn stored in Y1
! Y1 is inverted using LUD and backsubstitution

! G'(GCG'+Cee) is stored in Y2
! zz stores (Zc - (Aedot_pr))
 
! The total CG'(GCG'+Cee) (Zc- (Aedot_pr)) is stored in BB

! Next the dummy points are computed, for the dummy points .....

! edot = edot_pr + C~^G'(ACG'+Cee) (Zc - (Aedot_pr))
!~ represents a parameter at a control point
!^ represents a parameter at a data point
!Cmm is a row of C^~ dscribing covariance of control point with data points

! H is the inverse operator used to calculate the exhumation rates at the data points
! this is used in other subroutines later

use definitions
use omp_lib

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
integer i,j,k,INFO,ik,jk
integer id,nproc,chunk,first,last
double precision dist,resid,resim
double precision,dimension(:),allocatable::Cmm,zp

allocate(matrices%Y1(params%n,params%n))
allocate(matrices%B(params%n),matrices%II(params%n,params%n),matrices%work(100*params%n)) 

! computes Y1=GCG'+Cee
! firstly it computes the exhumation at the dummy, and then at data points

! note Y2=CA'
!print*,"calculating Y2..."

matrices%Y2=0.d0
call dgemm('N','T',params%n*params%m_max,params%n,params%n*params%m_max,1.d0,matrices%cov, & 
           params%n*params%m_max,matrices%G,params%n,0.d0,matrices%Y2,params%n*params%m_max)

matrices%Y1=0.d0
call dgemm('N','N',params%n,params%n,params%n*params%m_max,1.d0,matrices%G,params%n,matrices%Y2,& 
           params%n*params%m_max,0.d0,matrices%Y1,params%n)

do i=1,params%n
if (params%iflag_error.eq.1) matrices%Y1(i,i)=matrices%Y1(i,i)+(1.d0*0.35d0)**2.
if (params%iflag_error.eq.2) matrices%Y1(i,i)=matrices%Y1(i,i)+(1.d0*0.7d0)**2.
if (params%iflag_error.eq.3) matrices%Y1(i,i)=matrices%Y1(i,i)+(matrices%a_error(i)*0.35d0)**2.
if (params%iflag_error.eq.4) matrices%Y1(i,i)=matrices%Y1(i,i)+(matrices%a_error(i)*params%edot_mean)**2.
enddo

matrices%B=0

call dgetrf(params%n,params%n,matrices%Y1,params%n,matrices%B,INFO)
call dgetri(params%n,matrices%Y1,params%n,matrices%B,matrices%work,100*params%n,INFO)

matrices%Y2=0.d0
call dgemm('T','N',params%n*params%m_max,params%n,params%n,1.d0,matrices%G,params%n, &
           matrices%Y1,params%n,0.d0,matrices%Y2,params%n*params%m_max)

matrices%zz=0.d0
call dgemv('N',params%n,params%n*params%m_max,1.d0,matrices%A,params%n,matrices%edot_pr,1,0.d0,matrices%zz,1)

matrices%zz=(matrices%zc+matrices%elev)-(matrices%zz)

allocate(matrices%BB(params%n*params%m_max))

matrices%BB=0.d0
call dgemv('N',params%n*params%m_max,params%n,1.d0,matrices%Y2,params%n*params%m_max,matrices%zz,1,1.d0,matrices%BB,1)

matrices%edot=0.d0

!$omp parallel private(i,j,k,ik,jk,id,first,last,chunk,dist,Cmm)

  id=omp_get_thread_num()
  nproc=omp_get_num_threads()
  chunk=floor((params%m_max*params%dummy)/float(nproc))

  first=(id*chunk)+1
  last=(id+1)*chunk
  if (id.eq.nproc-1) last=last+mod((params%m_max*params%dummy),nproc)

  allocate(Cmm(params%n*params%m_max))
 
  do i=first,last

  Cmm=0.
  ik = floor(float(i-1)/float(params%m_max)) +1
  do j=1,params%n*params%m_max
    jk = floor(float(j-1)/float(params%m_max))+1
    if (mod(i,params%m_max).eq.mod(j,params%m_max)) then
      dist=sqrt((matrices%x_dum(ik)-matrices%x(jk))**2+(matrices%y_dum(ik)-matrices%y(jk))**2)
      Cmm(j)=params%sigma2*exp(-(dist/params%xL))    
    endif
  enddo

  do k=1,params%n*params%m_max
    matrices%edot(i) = matrices%edot(i)+Cmm(k)*matrices%BB(k)
  enddo
 
 matrices%edot(i)=matrices%edot(i)+matrices%edot_pr_dum(i)
enddo

  deallocate(Cmm)

!$omp end parallel

where(matrices%edot.lt.0.) matrices%edot=0.

matrices%H=0.d0
call dgemm('N','N',params%n*params%m_max,params%n,params%n*params%m_max,1.d0,matrices%cov,&
           params%n*params%m_max,matrices%Y2,params%n*params%m_max,0.,matrices%H,params%n*params%m_max)

! save output for weight Matric H
open(890,file=params%run//"/res/inverse_matrix.txt",status="unknown")
do i=1,params%m_max*params%n
     write(890,*) matrices%H(i,:)
enddo
close(890)
open(891,file=params%run//"/res/data_loc.txt",status="unknown")
do i=1,params%n
    write(891,*) matrices%x_true(i),matrices%y_true(i)
enddo
close(891)

matrices%zz=0.d0
call dgemv('N',params%n,params%n*params%m_max,1.d0,matrices%A,params%n,matrices%edot_pr,1,0.d0,matrices%zz,1)

allocate(zp(params%n))
zp=matrices%zz
matrices%zz = (matrices%zc + matrices%elev) - (matrices%zz)
matrices%edot_dat=0.d0
call dgemv('N',params%n*params%m_max,params%n,1.d0,matrices%H,params%n*params%m_max,matrices%zz,1,1.d0,matrices%edot_dat,1)

matrices%edot_dat=matrices%edot_dat+matrices%edot_pr
!matrices%edot_dat=matrices%edot_pr

where(matrices%edot_dat.lt.0.) matrices%edot_dat=0.

deallocate(matrices%Y1,matrices%B,matrices%II,matrices%BB,matrices%work)

return
end subroutine solve_inverse

