subroutine posterior_full(matrices,params)
!call posterior(matrices,params)

!This routine calculates the posterior variance at the control points
!and temporal resolution

!Cee is the posterior variance at the control points
!~ represents a parameter at a control point
!^ represents a parameter at a data point
!C~~ is the covariance at control points

!G series of steps is used to calculate Cee outlined below

!Cee = C~~ - C~^G'(GCG'+Ce)**-1 GC^~
!                 [    =Y1    ]
!Cee = C~~ - C~^ G'Y1 G C^~
!               [ =Y2]
!Cee = C~~ - C~^ Y2 G C^~
!                [Y3]

!In order not to build a full Cee matrix, which could be enormous ...
!  we only calculate the diagonal terms of Cee

!Cm is a row of C^~ dscribing covariance of control point with data
!points

!Cee(i,i) = C(i,i) - Ci^  Y3 C^i
!where Ci^ is the covariance of cont point i with data points

!Cee(i,i) = C(i,i) - Ci^  Y4
!Cee(i,i) = C(i,i) - Y5

!Different processors handle different control points

!For the resolution matrix

!R = C~^G'(GCG'+Ce)**-1 G 
!similarly a row of the resolution matrix

!R = Ci^ Y3

!Temporal resolution is calculated by integrating in space across time

use definitions
use omp_lib

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
double precision dist,tmp1,cdum,jtime,rest,rest2
integer i,j,ik,jk,k,itime,kk
integer ind,id,nproc,chunk,first,last
double precision,dimension(:,:),allocatable::Y4,Cm,Cee_pr
!real*4,dimension(:,:),allocatable::Y4,Cm,Cee_pr

!Y2 has been calculated in solve inverse
!Y2 = G' Y1
!Y3 = Y2 G
allocate(matrices%Y3(params%n*params%m_max,params%n*params%m_max))

allocate(matrices%Cee_full(params%dummy*params%m_max,params%dummy*params%m_max))
allocate(Cee_pr(params%dummy*params%m_max,params%dummy*params%m_max))

matrices%Y3=0.d0
call dgemm('N','N',params%n*params%m_max,params%n*params%m_max,params%n,1.d0,matrices%Y2,& 
           params%n*params%m_max,matrices%G,params%n,0.d0,matrices%Y3,params%n*params%m_max)

allocate(Cm(params%dummy*params%m_max,params%n*params%m_max))

  do i=1,params%dummy*params%m_max
  ik = floor(float(i-1)/float(params%m_max)) +1
  ind=0
  do j=1,params%n*params%m_max
    jk = floor(float(j-1)/float(params%m_max)) +1
    if (mod(i,params%m_max).eq.mod(j,params%m_max)) then
     if (ind.eq.0) ind = j   
       dist=sqrt((matrices%x_dum(ik)-matrices%x(jk))**2+(matrices%y_dum(ik)-matrices%y(jk))**2)
       Cm(i,j)=params%sigma2*exp(-(dist/params%xL))    
     endif
  enddo
  enddo

!print*,"done Cm"
  
allocate(Y4(params%dummy*params%m_max,params%n*params%m_max))
Y4=0.d0
call dgemm('N','N',params%dummy*params%m_max,params%n*params%m_max,params%n*params%m_max,1.d0,Cm,& 
           params%dummy*params%m_max,matrices%Y3,params%n*params%m_max,0.d0,Y4,params%dummy*params%m_max)

!print*,"Y4 done"

matrices%Cee_full=0.d0
call dgemm('N','T',params%dummy*params%m_max,params%dummy*params%m_max,params%n*params%m_max,1.d0,Y4,& 
           params%dummy*params%m_max,Cm,params%dummy*params%m_max,0.d0,matrices%Cee_full,params%dummy*params%m_max)
   
!print*,"Cee done"

do i=1,params%dummy*params%m_max
ik = floor(float(i-1)/float(params%m_max)) +1
do j=i,params%dummy*params%m_max 
   cdum=0.d0 
   jk = floor(float(j-1)/float(params%m_max)) +1
   if (mod(i,params%m_max).eq.mod(j,params%m_max)) then
      dist=sqrt((matrices%x_dum(ik)-matrices%x_dum(jk))**2+(matrices%y_dum(ik)-matrices%y_dum(jk))**2)
      cdum=params%sigma2*exp(-(dist/params%xL))    
   endif 
   matrices%Cee_full(i,j) = (cdum - matrices%Cee_full(i,j))
   matrices%Cee_full(j,i) = matrices%Cee_full(i,j) 
   Cee_pr(i,j) = cdum
   Cee_pr(j,i) = Cee_pr(i,j)
   enddo
enddo

deallocate(Cm,Y4)

!print*,"posterior matrices computed!"

open(88,file=params%run//"/full_cpost.txt",status="unknown")
open(91,file=params%run//"/full_cpostOFFDIAG.txt",status="unknown")
open(89,file=params%run//"/full_locs.txt",status="unknown")
open(87,file=params%run//"/full_edot.txt",status="unknown")

do k=params%m_max-2,params%m_max,2
!do k=params%m_max,params%m_max

        rest=sum(matrices%tsteps(1:params%m_max-k+1))
        rest2=sum(matrices%tsteps(1:params%m_max-k)) 

!print*,rest,"-",rest2,"Ma"

do i=1,params%dummy
if (k.eq.params%m_max) write(89,*) matrices%x_dum_true(i),matrices%y_dum_true(i)

write(87,*) matrices%edot(k+(i-1)*params%m_max),params%edot_mean
  do j=1,params%dummy
   write(88,*) matrices%Cee_full(k+(i-1)*params%m_max,k+(j-1)*params%m_max),Cee_pr(k+(i-1)*params%m_max,k+(j-1)*params%m_max)
  enddo
  if (k.eq.params%m_max-2) j=params%m_max
  if (k.eq.params%m_max) j=params%m_max-2
  write(91,*) matrices%Cee_full(k+(i-1)*params%m_max,j+(i-1)*params%m_max),Cee_pr(k+(i-1)*params%m_max,j+(i-1)*params%m_max)
enddo
enddo

deallocate(matrices%Cee_full,Cee_pr)
deallocate(matrices%Y3)

return
end subroutine posterior_full
