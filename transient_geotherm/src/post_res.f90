subroutine posterior(matrices,params)

!This routine calculates the posterior variance at the control points and temporal resolution

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

!Cm is a row of C^~ dscribing covariance of control point with data points

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
double precision dist,tmp1,Y5,cdum,jtime
integer i,j,ik,jk,k,itime
integer ind,id,nproc,chunk,first,last
double precision,dimension(:),allocatable::Y4,Cee,Cm,t_res,R_row,time_row

!Y2 has been calculated in solve inverse
!Y2 = G' Y1
!Y3 = Y2 G
allocate(matrices%Y3(params%n*params%m_max,params%n*params%m_max))

matrices%Y3=0.d0
call dgemm('N','N',params%n*params%m_max,params%n*params%m_max,params%n,1.d0,matrices%Y2, & 
           params%n*params%m_max,matrices%G,params%n,0.d0,matrices%Y3,params%n*params%m_max)

cdum=params%sigma2*exp(-sqrt((0.d0/params%xL)**2.))

!Y4 = C\tilde{m}\hat{m} Y3
!Cmm = C\hat{m}\tilde{m}

allocate(Cee(params%dummy*params%m_max),t_res(params%dummy*params%m_max))

matrices%sf=0.d0
t_res=0.d0

! The next line is required for the parallel implementation
!$omp parallel private(i,j,k,ik,jk,ind,id,first,last,R_row,time_row,itime,jtime,tmp1,Y5,Cm,Y4) shared(Cee,cdum,t_res)

  allocate(Cm(params%n*params%m_max))
  allocate(Y4(params%n*params%m_max))
  allocate(R_row(params%n*params%m_max))
  allocate(time_row(params%m_max))
   
  id=omp_get_thread_num()
  nproc=omp_get_num_threads()
  chunk=floor((params%m_max*params%dummy)/float(nproc))

  first=(id*chunk)+1
  last=(id+1)*chunk
  if (id.eq.nproc-1) last=last+mod((params%m_max*params%dummy),nproc)

  ! first and last are the indices of the first and last control points ...
  ! ... that a processor has to `look after'
  do i=first,last
 
  ! Cm contains the covariance values for the parameter (i) with all the parameters at the data points 
  Cm=0.d0
  ik = floor(float(i-1)/float(params%m_max)) +1
  ind=0
  do j=1,params%n*params%m_max
    jk = floor(float(j-1)/float(params%m_max)) +1
    if (mod(i,params%m_max).eq.mod(j,params%m_max)) then
     if (ind.eq.0) ind = j   
       dist=sqrt((matrices%x_dum(ik)-matrices%x(jk))**2+(matrices%y_dum(ik)-matrices%y(jk))**2)
       Cm(j)=params%sigma2*exp(-(dist/params%xL))    
     endif
  enddo

  ! First multiplication to get Y4, see line 27 
  Y4=0.d0
  ! loop in steps of params%m_max, cross timestep terms are 0
  do j=ind,params%n*params%m_max,params%m_max
    do k=1,params%n*params%m_max
      Y4(j)= Y4(j) + matrices%Y3(k,j)*Cm(k)
    enddo
  enddo
   
  ! calculate row of resolution matrix
  R_row=0.
  do k=1,params%m_max*params%n
    do j=ind,params%m_max*params%n,params%m_max
      R_row(k) = R_row(k) + Cm(j)*matrices%Y3(k,j)
    enddo
  enddo
   
  ! and integrate in space for each time intervals
  time_row=0.
  itime = i-((ik-1)*params%m_max)
  do j=1,params%m_max!-1 
    jk = floor(float(j-1)/float(params%m_max))+1      
    jtime = j-((jk-1)*params%m_max)
      do k=1,params%n
        time_row(j) = time_row(j) + R_row(j+(k-1)*params%m_max)
      enddo
    t_res(i) = time_row(itime) 
  enddo

  !Calculate Y5 see, line 28 
  Y5=0.d0
  !Y5 = Y4 Cmm'
   do k=ind,params%n*params%m_max,params%m_max
     Y5 = Y5 + Cm(k)*Y4(k)
   enddo

   Cee(i) = sqrt(cdum - Y5)

  enddo

 deallocate(Cm,Y4)
 deallocate(R_row,time_row)

!$omp end parallel

matrices%eps_dum=Cee
matrices%sf=t_res

print*,"posterior matrices computed!"

deallocate(Cee,t_res,matrices%Y3)

return
end subroutine posterior

