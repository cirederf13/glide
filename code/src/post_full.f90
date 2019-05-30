subroutine posterior_full(matrices,params)
!call posterior(matrices,params)

! this routine calculates the full posterior covariance matrix

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

!Y2 has been calculated in solve inverse
!Y2 = G' Y1
!Y3 = Y2 G
allocate(matrices%Y3(params%n*params%m_max,params%n*params%m_max))

allocate(matrices%Cee_full(params%dummy*params%m_max,params%dummy*params%m_max))
allocate(Cee_pr(params%dummy*params%m_max,params%dummy*params%m_max))

matrices%Y3=0.
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
  
allocate(Y4(params%dummy*params%m_max,params%n*params%m_max))
Y4=0.
call dgemm('N','N',params%dummy*params%m_max,params%n*params%m_max,params%n*params%m_max,1.d0,Cm,& 
           params%dummy*params%m_max,matrices%Y3,params%n*params%m_max,0.d0,Y4,params%dummy*params%m_max)

matrices%Cee_full=0.
call dgemm('N','T',params%dummy*params%m_max,params%dummy*params%m_max,params%n*params%m_max,1.d0,Y4,& 
           params%dummy*params%m_max,Cm,params%dummy*params%m_max,0.d0,matrices%Cee_full,params%dummy*params%m_max)

do i=1,params%dummy*params%m_max
ik = floor(float(i-1)/float(params%m_max)) +1
do j=i,params%dummy*params%m_max 
   cdum=0. 
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

print*,"posterior matrices computed!"

deallocate(matrices%Cee_full,Cee_pr)
deallocate(matrices%Y3)

return
end subroutine posterior_full
