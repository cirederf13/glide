subroutine posterior_dat(matrices,params,post)

use definitions

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
integer i,j,k,post
double precision,dimension(:,:),allocatable::Y4,Y5,Cres,Cdata,Y6,Cee

allocate(matrices%Y3(params%n,params%n*params%m_max))

!Y3=GC

matrices%Y3=0.
call dgemm('N','N',params%n,params%n*params%m_max,params%n*params%m_max,1.d0,matrices%G,params%n, &
           matrices%cov,params%n*params%m_max,0.d0,matrices%Y3,params%n)

!Cpost=C-HGC
!     = C-HY3
matrices%cpost=0.
call dgemm('N','N',params%n*params%m_max,params%n*params%m_max,params%n,1.d0,matrices%H,params%n*params%m_max, &
           matrices%Y3,params%n,0.d0,matrices%cpost,params%n*params%m_max)

matrices%cpost=matrices%cov-matrices%cpost

matrices%eps=0.
do i=1,params%n*params%m_max
  matrices%eps(i)=sqrt(matrices%cpost(i,i))
enddo

! deallocates the working tables...
deallocate(matrices%Y3)

!post=1

if (post.eq.1) then !flag to calculate the matrices

!Next i have to compute Cres = (HG-I)C(HG-I)'
!Y4 = HG-I 
print*,"calculating Y4 ... "
allocate(Y4(params%n*params%m_max,params%n*params%m_max))

do j=1,params%n*params%m_max
 do k=1,params%n
  Y4(:,j)=Y4(:,j)+matrices%G(k,j)*matrices%H(:,k)
 enddo
enddo 

do i=1,params%n*params%m_max
 Y4(i,i)=Y4(i,i)-1.
enddo 

!where (Y4.lt.0) Y4=0.
print*,"calculating Y5 ... "
!!!! Y5 = (HG-I)C=Y4C
allocate(Y5(params%n*params%m_max,params%n*params%m_max))

Y5=0.
do j=1,params%n*params%m_max
 do k=1,params%n*params%m_max
  Y5(:,j)=Y5(:,j)+matrices%cov(k,j)*Y4(:,k)
 enddo
enddo

print*,"calculating Cres ... "
!!Cres=C(HG-I)'=CY4
allocate(Cres(params%n*params%m_max,params%n*params%m_max))

do j=1,params%n*params%m_max
 do k=1,params%n*params%m_max
  Cres(:,j)=Cres(:,j)+Y4(j,k)*Y5(:,k) !note I use (j,k) beause it is the transpose
 enddo
enddo

!Cdata=HCeH'
allocate(Y6(params%n*params%m_max,params%n),cdata(params%n*params%m_max,params%n*params%m_max),Cee(params%n,params%n))


print*,"calculating Y6... "
Cee=0.
do i=1,params%n
 Cee(i,i) = (matrices%a_error(i)*params%edot_mean)**2
enddo

!subroutine SGEMM(n,n,rowsG,columnB,common,1,G(LDG,:),LDG,B(LDG,:),LDB,)
Y6=0.
do j=1,params%n
 do k=1,params%n
  Y6(:,j)=Y6(:,j) + Cee(k,j)*matrices%H(:,k)
 enddo
enddo

print*,"calculating Cdata ... "
Cdata=0.
do j=1,params%n*params%m_max
 do k=1,params%n
  Cdata(:,j) = Cdata(:,j)+matrices%H(j,k)*Y6(:,k) !H is transpose!
 enddo
enddo

matrices%eps_res=0.
matrices%eps_dat=0.
do i=1,params%n*params%m_max
 matrices%eps_res(i)=sqrt(Cres(i,i)) !purple
 matrices%eps_dat(i)=sqrt(Cdata(i,i)) !blue
enddo

deallocate(Y6,cdata,Cee,Cres,Y4,Y5)

open(98,file=params%run//"/res/data_resolution.txt",status="unknown")
do i=1,params%n*params%m_max
 write(98,*) matrices%eps_res(i)
enddo
close(98)

endif

return
end subroutine posterior_dat
