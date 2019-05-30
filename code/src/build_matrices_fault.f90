subroutine build_matrices(matrices,params)
! This subroutine first builds matrices A, cov, cov_eps

use definitions

implicit none

type (matr) matrices
type (parm) params

integer i,j,k,m,ik,jk,nfault,kk
real rest,dist,fist,fdis,a1,a2,b1,b2,cross1,cross2,angle,cos1,sin1,summ
real,dimension(:),allocatable::flat,flon,xold,yold,xolder,yolder
real,dimension(:,:),allocatable::test

!fault file, of length nfault lon,lat
!open(73,file="data/corsica/fault.gmt")
!nfault = 111
nfault=0

allocate(matrices%A(params%n,params%n*params%m_max),&
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
allocate(xold(params%n),yold(params%n),xolder(params%n),yolder(params%n))

allocate(flat(nfault),flon(nfault),matrices%distm(params%n,params%n))

! It builds the matrix A (i.e. the discretized ages)

matrices%A=0.
print*,"building A"

     do i=1,params%n

     m=1
     k=1
     summ=0.
      do while (summ.lt.matrices%ta(i))
       summ = summ + matrices%tsteps(m)
       m=m+1
       enddo
       m=m-1
       if (m.eq.0) then
         summ = 0.
         rest = 0.
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

!Faults can be coded like this but better if they are straight
!do i=1,nfault
! read(73,*) flon(i),flat(i) 
! flon(i) = (flon(i)-params%lon1)*111.11*cos(((params%lat1+params%lat2)/2.)*3.141592654/180.)
! flat(i) = (flat(i)-params%lat1)*111.11
!enddo

!do i=1,params%n
!  xolder(i)=matrices%x(i)
!  yolder(i)=matrices%y(i)
!enddo

xolder=matrices%x
yolder=matrices%y

if(params%angle.ne.0.or.params%aspect.ne.0.) then

	matrices%x=matrices%x-(maxval(matrices%x)-minval(matrices%x))/2.
	matrices%y=matrices%y-(maxval(matrices%y)-minval(matrices%y))/2.
	
	!do i=1,params%n
	!  xold(i)=matrices%x(i)
	!  yold(i)=matrices%y(i)
	!enddo
	
        xold=matrices%x
        yold=matrices%y

	!angle used to test define anistropy
	!params%angle=17.
	print*,"anistropy angle=",params%angle,"can be changed in build_matrices.f90"
	params%angle=params%angle*(3.141592654/180.)
         
        cos1=cos(params%angle)
        sin1=sin(params%angle)
!	do i=1,params%n
!	!	rotating 
!		matrices%x(i) = (xold(i)*cos1) - (yold(i)*sin1)
!	    matrices%y(i) = (xold(i)*sin1) + (yold(i)*cos1)
!	!   stretching 
!	    matrices%y(i)=matrices%y(i)*params%aspect
!	enddo	
        !rotating
        matrices%x=(cos1*xold)-(sin1*yold)
        matrices%y=(sin1*xold)+(cos1*yold)
        !stretching
        matrices%y=matrices%y*params%aspect
else
	continue
endif

! It builds the covariance matrices defining the spatial correlation

!!$OMP PARALLEL PRIVATE(i,j,ik,jk)
!!$OMP DO SCHEDULE(STATIC,500)
matrices%cov=0.0
ik=0
do i=1,params%n*params%m_max
     ik=(i+params%m_max-1)/params%m_max
      jk=ik-1
      do j=i,params%n*params%m_max,params%m_max
       jk=jk+1
       dist=sqrt((matrices%x(ik)-matrices%x(jk))**2+(matrices%y(ik)-matrices%y(jk))**2)
       	
       	!go along fault, find which point on the fault is between     

!        fist=1000000.
!        fdis=0.
!       	kk=0
!        do k=1,nfault
!         fdis=(sqrt((matrices%x(ik)-flon(k))**2.+(matrices%y(ik)-flat(k))**2.)+sqrt((matrices%x(jk)-flon(k))**2.+(matrices%y(jk)-flat(k))**2.))
!         if (fdis.lt.fist) then
!          kk=k
!		  fist=fdis
!		 endif
!        enddo
!
!        !take the cross product 
!        if (kk.gt.1) then
!        a1=flon(kk)-flon(kk-1)
!        a2=flat(kk)-flat(kk-1)
!        b1=matrices%x(ik)-flon(kk)
!        b2=matrices%y(ik)-flat(kk)
!        cross1 = a1*b2-a2*b1
!        b1=matrices%x(jk)-flon(kk)
!        b2=matrices%y(jk)-flat(kk)
!        cross2 = a1*b2-a2*b1
!        else
!        a1=flon(kk)-flon(kk+1)
!        a2=flat(kk)-flat(kk+1)
!        b1=matrices%x(ik)-flon(kk)
!        b2=matrices%y(ik)-flat(kk)
!        cross1 = a1*b2-a2*b1
!        b1=matrices%x(jk)-flon(kk)
!        b2=matrices%y(jk)-flat(kk)
!        cross2 = a1*b2-a2*b1
!        endif
!        
!        !if the cross products have different signs, set dist far away
!        if (sign(cross1,cross2).ne.cross1) dist=max(dist,2.*params%xL)
!	
       !number over dist is the corelation distance
!       matrices%cov(i,j)=params%sigma2*exp(-sqrt((dist/params%xL)**2.))
!       matrices%cov(i,j)=params%sigma2*exp(-(dist/params%xL)**2.)
        matrices%cov(i,j)=params%sigma2*exp(-(dist/params%xL))
!       if(i.lt.params%n*params%m_max) matrices%cov(i+1,j)=0.5*matrices%cov(i,j)
!       if(i.gt.1) matrices%cov(i-1,j)=0.5*matrices%cov(i,j)
!       matrices%cov(i,j)=matrices%distm(ik,jk) 
!       if(j.lt.params%n*params%m_max) matrices%cov(j+1,i)=0.5*matrices%cov(i,j)
!       if(j.gt.1) matrices%cov(j-1,i)=0.5*matrices%cov(i,j)
       matrices%cov(j,i)=matrices%cov(i,j)
      enddo
     enddo
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

print*,"covariance matrix complete"!"mincov",minval(matrices%cov),"maxcov",maxval(matrices%cov)

!these lines are required for anisotropic correlation
matrices%x=xolder
matrices%y=yolder

deallocate(flat,flon,xold,yold,xolder,yolder,matrices%distm)

return
end subroutine build_matrices
