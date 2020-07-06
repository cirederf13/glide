subroutine find_zc(matrices,params)

!This subroutine calculates closure depths 

! The routine proceeds as follows

! For each sample do  
! Set up initial conditions, a linear increase of temperature with depth ...
!  ... use crank-nicholson finite differencing to step through time  
!  ... at the time equivalent to the measured age record the materical derivatives
!  ... use these cooling rates to estimate closure temperatures
!  ... find location in depth where closure depth is equal to temperature

use definitions

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
integer i,j,mz,nt,iflag,kk,skip,exit_flag,pos

double precision,dimension(:),allocatable::temp,temp2
double precision Tc,b_flux,tdot_dist
double precision dz,xjunk,M,P,exhum
double precision tstart,tend,dt,xtime,ax,bx,cx,alpha,lambda
double precision,dimension(:),allocatable::diag,sup,inf,f,temp_pr,temp_age,tdot,closure,query

! parameters defining the geometry
mz=131
dz=params%zl/float(mz-1)

allocate(temp(mz),temp2(mz))

allocate(diag(mz),sup(mz),inf(mz),f(mz),temp_pr(mz),temp_age(mz),tdot(mz))

allocate(closure(mz),query(30))

skip=minloc(matrices%ta,1)

! It solves the heat transfer equation in 1D to include a transient solution
! it does it at each data point
do kk=1,params%n

call random_seed()

! goes forward in time and solves the heat transfer equation

dt = (dz**2.)/params%kappa/4.1

tstart=0.
tend=params%t_total 
nt=int(tend/dt)+1
temp=0.
temp2=0.

xtime=tstart
j=0

!set up initial conditions
do i=1,mz
  temp(i)=params%Ts+float(i-1)/float(mz-1)*(params%Tb-params%Ts)
enddo

!basal flux
b_flux=(temp(mz)-temp(mz-1))

!exit if exit_flag=1
exit_flag=0

!loop through time
do while (xtime.lt.tend)

   query=xtime
   query=matrices%tsteps_sum-query
   where(query.lt.0.) query=-9999.
   pos=minloc((query)**2.,1)
   exhum=matrices%edot_dat(pos+(kk-1)*params%m_max)

   !constants used by tridag
   lambda=params%kappa*dt/(2.*(dz**2))
   alpha=exhum*dt/(4.*dz)
   ax=(alpha-lambda)
   bx=(1+(2.*lambda))
   cx=-1.*(lambda+alpha)

   !diag=diag
   !sup=above
   !inf=below
    do i=2,mz-1
      diag(i)=bx
      sup(i)=cx
      inf(i)=ax
      f(i)=(lambda-alpha)*temp(i-1) + &
           (1-(2.*lambda))*temp(i) + &
           (lambda+alpha)*temp(i+1) + params%hp*dt
    enddo
     
    !apply boundary conditions
    diag(1)=1.
    sup(1)=0.
    f(1)=params%Ts
    diag(mz)=1.
    inf(mz)=0.
    f(mz)=temp(mz-1)+b_flux

    ! solve linear equations
!    call tridag(inf,diag,sup,f,temp,mz)

    temp(1)=params%Ts
    temp(mz)=temp(mz-1)+b_flux
         
    ! increase time
    xtime=xtime+dt
    j=j+1
    
    !when the sample age of the model equals the measured age, exit
    if (exit_flag.eq.1) then
      temp_age=temp
      exit
    endif  

    !time step before the age is reached, change dt
    if (tend-xtime-dt.lt.matrices%ta(kk)) then
      dt=tend-xtime-matrices%ta(kk)
      lambda=params%kappa*dt/(2.*(dz**2))
   
      !update exhum
      query=xtime
      query=matrices%tsteps_sum-query
      where(query.lt.0.) query=-9999.
      pos=minloc((query)**2.,1)
      exhum=matrices%edot_dat(pos+(kk-1)*params%m_max)

      !define constants for timstepping 
      alpha=exhum*dt/(4.*dz)
      ax=(alpha-lambda)
      bx=(1+(2.*lambda))
      cx=-1.*(lambda+alpha)
      temp_pr=temp   
      exit_flag=1
    endif

enddo 

tdot_dist=exhum*dt

!calculate material derivative,tdot
do i=1,mz-1
  tdot(i)=temp_age(i)-temp_pr(i)+(temp_pr(i+1)-temp_pr(i))*((tdot_dist)/dz)
  tdot(i)=tdot(i)/dt
enddo

!force tdot to be positive
!tdot=max(tdot,0.01)

!calculate tdot dependant closure depths
closure=0.
iflag=matrices%isys(kk)
call closure_temps (tdot,temp_age,closure,mz,iflag)

do i=1,mz-1
  !Solve the simultaneous equations to calculate, z_c,T_c
  if (temp_age(i).gt.closure(i)) then
     M=1.*(dz/(closure(i)-closure(i-1)))
     P=1*(dz/(temp_age(i)-temp_age(i-1)))
     Tc=((M*closure(i-1)) - (P*temp_age(i-1)))/(M-P)
     xjunk=M*(Tc-closure(i))
   
      matrices%zc(kk)=dz*float(i-1)+xjunk
    exit
   endif
enddo

!ending loop through samples
enddo

return
end subroutine
