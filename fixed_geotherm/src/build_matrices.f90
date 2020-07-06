subroutine build_matrices(matrices,params)

! This subroutine builds the forward operator, A and 
! the covariance matrix, cov
! It also allocates matrices used in solve_inverse

use definitions

implicit none

type (matr) matrices
type (parm) params

integer i,j,k,m,ik,jk
double precision rest,dist,summ

allocate(matrices%A(params%n,params%n*params%m_max),&
         matrices%G(params%n,params%n*params%m_max),&
         matrices%edot(params%dummy*params%m_max),&
         matrices%edot_dat(params%n*params%m_max),&
         matrices%eps_dum(params%dummy*params%m_max),&
         matrices%eps(params%n*params%m_max),&
         matrices%sf(params%dummy*params%m_max),&
         matrices%cov(params%n*params%m_max,params%n*params%m_max),&
         matrices%cpost(params%n*params%m_max,params%n*params%m_max))
!         matrices%cov_eps(params%n), &

allocate(matrices%eps_dat(params%n*params%m_max),&
         matrices%eps_res(params%n*params%m_max),&
         matrices%H(params%n*params%m_max,params%n),&
         matrices%Y2(params%n*params%m_max,params%n))

! Build the matrix A (i.e. the discretized ages)

matrices%A=0.d0
matrices%G=0.d0

do i=1,params%n
  m=1
  k=1
  summ=0.d0
    !calculate the number of time steps required and the remainder term
    do while (summ.lt.matrices%ta(i))
      summ = summ + matrices%tsteps(m)
      m=m+1
    enddo
    m=m-1
    if (m.eq.0) then
       summ = 0.d0
       rest = 0.d0
    else
       summ = summ - matrices%tsteps(m)
       rest = matrices%ta(i) - summ 
    endif
    k=0
    do j=(params%m_max-m)+(i-1)*(params%m_max)+1,(i-1)*(params%m_max-1)+params%m_max+(i-1)
      k=k+1
      matrices%A(i,j)=matrices%tsteps(m-k+1)
      if (k.eq.1) matrices%A(i,j)=rest
    enddo
enddo
      matrices%G=matrices%A

! It builds the covariance matrices defining the spatial correlation
matrices%cov=0.d0
ik=0
!loop through lines of matrix
do i=1,params%n*params%m_max
  ik=(i+params%m_max-1)/params%m_max
  jk=ik-1
  do j=i,params%n*params%m_max,params%m_max
    jk=jk+1
    !calculate distances between points
    dist=sqrt((matrices%x(ik)-matrices%x(jk))**2+(matrices%y(ik)-matrices%y(jk))**2)
    matrices%cov(i,j)=params%sigma2*exp(-(dist/params%xL))
    matrices%cov(j,i)=matrices%cov(i,j)
  enddo
enddo

return
end subroutine build_matrices
