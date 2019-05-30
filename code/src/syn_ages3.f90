subroutine syn_ages (matrices,params,misfit)

use definitions
use omp_lib

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
integer i,j,k,p,mz,nt,iflag,first,last,chunk
integer id,nproc,pos,ntsteps

double precision Ts,Tb,kappa,hp,Tc
double precision dz,zl,xint
double precision tstart,tend,dt,xtime,misfit
double precision sample_age,cooling,frac,dist
double precision,dimension(:,:),allocatable::temp,temperature,age,diag,sup,inf,f,deep
double precision,dimension(:),allocatable::exhum,geotherm,ax,bx,cx,lambda,alpha,query,b_flux

allocate(geotherm(params%n))


!$omp parallel private(i,j,k,nt,id,first,last,chunk,p,nproc,temp,exhum,xtime,&
!$omp & iflag,Tc,sample_age,xint,temperature,age,diag,sup,inf,f,alpha,lambda,&
!$omp & ax,bx,cx,query,pos,ntsteps,b_flux,deep,frac,dist)
!shared(geotherm,tstart,tend,dt,mz

  ntsteps=100
   
  !get thread number, and the total number of threads
  id=omp_get_thread_num()
  nproc=omp_get_num_threads()
  !chunk is the number of samples each proc 
  chunk=floor((params%n)/float(nproc))

  first=(id*chunk)+1
  last=(id+1)*chunk
  if (id.eq.nproc-1) last=last+mod((params%n),nproc)
  
  chunk=last-first+1
  
  mz=211
  dz=params%zl/float(mz-1)

allocate(temp(mz,chunk),exhum(chunk),ax(chunk),bx(chunk),cx(chunk),query(ntsteps))
allocate(diag(mz,chunk),sup(mz,chunk),inf(mz,chunk),f(mz,chunk),alpha(chunk),lambda(chunk),b_flux(chunk))

! It solves the heat transfer equation in 1D to include a transient solution
! it does it at each data point using a crank nicholson finite differencing scheme

!dt to give accurate results
dt = (dz**2.)/params%kappa/4.1

tstart=0.
tend=params%t_total
nt=int(tend/dt)+1
temp=0.

do i=1,mz
temp(i,:)=params%Ts+float(i-1)/float(mz-1)*(params%Tb-params%Ts)
enddo

b_flux=(temp(mz,:)-temp(mz-1,:))

! first tracks back the inital depth of each sample

!matrices%depth(first:last)=0.
matrices%depth(first:last)=-1.*matrices%elev(first:last)

xtime=tend
j=0
 do while(xtime.gt.tstart)
 j=j+1
   p=1
   query=xtime
   query=matrices%tsteps_sum-query
   where(query.lt.0.) query=-9999.
   pos=minloc((query)**2.,1)
   do i=first,last
!     exhum(p)=matrices%edot_dat(floor(xtime/params%deltat)+1+(i-1)*params%m_max)
     exhum(p)=matrices%edot_dat(pos+(i-1)*params%m_max)    
     matrices%depth(i)=matrices%depth(i)+exhum(p)*dt
    p=p+1
   enddo
   xtime=xtime-dt
 enddo
nt=j

allocate(temperature(nt,chunk),age(nt,chunk),deep(nt,chunk))

! goes forward in time and solves the heat transfer equation
! and tracks the Tt path
xtime=tstart
j=1
!!print*,'computing thermal model'
 do while (xtime.lt.tend)
    
    query=xtime
    query=matrices%tsteps_sum-query
    where(query.lt.0.) query=-9999.
    pos=minloc((query)**2.,1)
  !  print*,pos,floor(xtime/params%deltat)+1 
   !p is used as a dummy counter
   p=1
   do i=first,last
   ! exhum(p)=matrices%edot_dat(floor(xtime/params%deltat)+1+(i-1)*params%m_max)
    exhum(p)=matrices%edot_dat(pos+(i-1)*params%m_max)
    p=p+1
   enddo

   lambda(:)=params%kappa*dt/(2.*(dz**2))
   alpha(:)=exhum(:)*dt/(4.*dz)
   ax=(alpha-lambda)
   bx=(1+(2.*lambda))
   cx=-1.*(lambda+alpha)

    !diag=diag
    !sup=above
    !inf=below
    do p=2,mz-1
      sup(p,:)=cx 
      inf(p,:)=ax
      diag(p,:)=bx
      f(p,:)=(lambda-alpha)*temp(p-1,:) + &
           (1-(2.*lambda))*temp(p,:) + &
           (lambda+alpha)*temp(p+1,:) + params%hp*dt
    enddo
    
    diag(1,:)=1.
    sup(1,:)=0.
    f(1,:)=params%Ts !0.
    diag(mz,:)=1.
    inf(mz,:)=0.
    f(mz,:)=temp(mz-1,:)+b_flux !params%Tb
    !f(mz,:)=params%Tb
     
    call tridag_par(inf,diag,sup,f,temp,mz,chunk)
    
    !tridag called, and finished
    temp(1,:)=params%Ts   
    temp(mz,:)=temp(mz-1,:)+b_flux !temp(mz-1,:)+((dz*20.)/3.2)
    !temp(mz,:)=params%Tb

  age(j,:)=tend-xtime

! extract temperature for each data at each time step
  p=1
  do i=first,last
    matrices%depth(i)=matrices%depth(i)-exhum(p)*dt
    deep(j,p)=matrices%depth(i)
    k=max(1,int(matrices%depth(i)/dz)+1)
   ! k=int(matrices%depth(i)/dz)+1
   ! if (k.lt.0) !print*,'negative depth...',i,k
    if (k.ge.mz) then
      temperature(j,p)=params%Tb
    elseif (k.le.0) then
       stop !!print*,i,k
    else
      xint=(matrices%depth(i)-dz*float(k-1))/dz
      temperature(j,p)=temp(k+1,p)*xint+temp(k,p)*(1.-xint)
    endif
    p=p+1
  enddo
   
  xtime=xtime+dt
  j=j+1
enddo


! It goes through the solution to find the closure depth
p=1
do j=first,last
 iflag=matrices%isys(j)
 call dodson (temperature(:,p),age(:,p),nt,sample_age,Tc,cooling,iflag)
  do i=1,nt
    if (temperature(i,p).lt.Tc) then
      if (i.eq.1) then 
        matrices%zcp(j)=deep(i,p)
      else
        frac = (Tc-temperature(i,p))/(temperature(i-1,p)-temperature(i,p))
        matrices%zcp(j) = deep(i,p) + frac*(deep(i-1,p)-deep(i,p))
      endif 
      goto 1234  
    endif
  enddo   
1234 continue
geotherm(j)=(temp(2,p)-temp(1,p))/dz
p=p+1
enddo

deallocate(temp,temperature,age,exhum,ax,bx,cx,diag,sup,inf,f,alpha,lambda,query,b_flux,deep)

!$omp end parallel

matrices%zz = max(0.01,matrices%elev + matrices%zcp)
!calculate synthetic ages based on the closure depths

matrices%syn_age=0.

do i=1,params%n
  dist=0.
  j=params%m_max
  do while (dist.lt.matrices%zz(i))
    matrices%syn_age(i) = matrices%syn_age(i) + matrices%tsteps(params%m_max+1-j)
    dist = dist + matrices%tsteps(params%m_max+1-j)*matrices%edot_dat(j+(i-1)*params%m_max)
    j=j-1
  enddo
  j=j+1
  dist = dist - matrices%tsteps(params%m_max+1-j)*matrices%edot_dat(j+(i-1)*params%m_max)
  matrices%syn_age(i) = matrices%syn_age(i) - matrices%tsteps(params%m_max+1-j)
  frac = (matrices%zz(i) - dist)/matrices%edot_dat(j+(i-1)*params%m_max)
  matrices%syn_age(i) = matrices%syn_age(i) + frac
enddo

print*,"min syn_age=",minval(matrices%syn_age),"max syn_age=",maxval(matrices%syn_age)

!deallocate(geotherm)

matrices%misfits=abs(matrices%syn_age-matrices%ta)

misfit=0.
misfit=sum((matrices%misfits**2.)/(matrices%a_error**2.))
misfit=sqrt(misfit/(float(params%n)))

return
end subroutine syn_ages

