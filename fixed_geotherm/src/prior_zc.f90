subroutine prior_zc(matrices,params)

!This subroutine calculates closure depths 

! The routine proceeds as follows

! For each sample do  
! Set up initial conditions, a linear increase of temperature with depth ...
!  ... use crank-nicholson finite differencing to step through time  
!  ... at the time equivalent to the measured age record the materical derivatives
!  ... use these cooling rates to estimate closure temperatures
!  ... find location in depth where closure depth is equal to temperature

! Determine average closure depth for each system, ie AFT
! also record other parameters for use in isotherms

use definitions

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices
integer i,j,mz,nt,iflag,kk,skip,exit_flag,pos

double precision,dimension(:),allocatable::temp,temp2
double precision Tc,tdot_dist,b_flux
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

open(123,file=params%run//"/stuff/Tz.txt",status="unknown")
open(124,file=params%run//"/stuff/heat_flux.txt",status="unknown")
open(125,file=params%run//"/stuff/cross_over.txt",status="unknown")
open(357,file=params%run//"/stuff/anal.txt",status="unknown")

!matrices%ages contains closure depths, geothermal gradients at depth and at the surface ....
! see lines 191:196 for info  
matrices%ages=0.

!loop through data points
do kk=1,params%n

call random_seed()

! goes forward in time and solves the heat transfer equation
! at the time equivalent to the cooling age, the material derivatives 
! for all depths are stored in temp_pr and temp 

!dt to give accurate results
dt = ((dz**2)/params%kappa)/4.

tstart=0.d0
tend=params%t_total! 12. 
nt=int(tend/dt)+1
temp=0.d0
temp2=0.d0

xtime=tstart
j=0

do i=1,mz
temp(i)=params%Ts+float(i-1)/float(mz-1)*(params%Tb-params%Ts)
enddo

if (kk.eq.skip) print*, "Geotherm at start of model",(temp(2)-temp(1))/dz

b_flux=(temp(mz)-temp(mz-1))

exit_flag=0

do while (xtime.lt.tend)

   query=xtime
   query=matrices%tsteps_sum-query
   where(query.lt.0.) query=-9999.
   pos=minloc((query)**2,1)
   exhum=matrices%edot_pr(pos+(kk-1)*params%m_max)

   !constants used by tridag
   lambda=params%kappa*dt/(2.d0*(dz**2))
   alpha=exhum*dt/(4.d0*dz)
   ax=(alpha-lambda)
   bx=(1.d0+(2.d0*lambda))
   cx=-1.d0*(lambda+alpha)

   !diag=diag
   !sup=above
   !inf=below
    do i=2,mz-1
      diag(i)=bx
      sup(i)=cx
      inf(i)=ax
      f(i)=(lambda-alpha)*temp(i-1) + &
           (1.d0-(2.d0*lambda))*temp(i) + &
           (lambda+alpha)*temp(i+1) + params%hp*dt
    enddo
    
    diag(1)=1.d0
    sup(1)=0.d0
    f(1)=params%Ts 
    diag(mz)=1.d0
    inf(mz)=0.d0
    f(mz)=temp(mz-1)+b_flux
    !f(mz)=params%Tb

!    call tridag(inf,diag,sup,f,temp,mz)

    temp(1)=params%Ts
    temp(mz)=temp(mz-1)+b_flux
    !temp(mz)=params%Tb
         
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
      where(query.lt.0.) query=-9999.d0
      pos=minloc((query)**2,1)
      exhum=matrices%edot_pr(pos+(kk-1)*params%m_max)

      alpha=exhum*dt/(4.d0*dz)
      ax=(alpha-lambda)
      bx=(1.d0+(2.d0*lambda))
      cx=-1.d0*(lambda+alpha)
      temp_pr=temp   
      exit_flag=1
    endif

    !write out tranisent solution for youngest age   
    if (kk.eq.skip) then
      if(mod(j,500).eq.0) then
        !writing geotherm for gmt
        write(123,'(A,4X)',advance='NO') ">"
        write(123,'(A,4X)',advance='NO') "-Z"
        write(123,'(F10.3)',advance='NO') xtime
        write(123,*)
        do i=1,mz
          write(123,*) temp(i),float(i-1)*dz,xtime
        enddo
        !2.4 = thermal conductivity of crust
        write(124,*) xtime,2.4*(temp(2)-temp(1))/dz,xtime
      endif
    endif
               
!ending time stepping
enddo 
    
if (kk.eq.skip) print*, "Geotherm at end of model",(temp(2)-temp(1))/dz

!tdot_dist is the total distance travelled between when temp_pr and temp 
! were calculated
tdot_dist=exhum*dt

!calculate material dervative,tdot
do i=1,mz-1
  tdot(i)=temp_age(i)-temp_pr(i)+(temp_pr(i+1)-temp_pr(i))*((tdot_dist)/dz)
  tdot(i)=tdot(i)/dt
enddo

!calculate tdot dependant closure depths
closure=0.d0
iflag=matrices%isys(kk)
call closure_temps (tdot,temp_age,closure,mz,iflag)
!call linear_Tc (tdot,temp_age,closure,mz,iflag)

do i=1,mz-1
  !Solve the simultaneous equations to calculate, z_c,T_c
  if (temp_age(i).gt.closure(i)) then
    
     ! M and P are coefficients for the simultaneous equations
     M=1.d0*(dz/(closure(i)-closure(i-1)))
     P=1.d0*(dz/(temp_age(i)-temp_age(i-1)))
     Tc=((M*closure(i-1)) - (P*temp_age(i-1)))/(M-P)
     xjunk=M*(Tc-closure(i))
   
     matrices%ages(iflag,1) = matrices%ages(iflag,1) + Tc                                !closure temp
     matrices%ages(iflag,2) = matrices%ages(iflag,2) + dz*float(i-1)+xjunk               !closure depth
     matrices%ages(iflag,3) = matrices%ages(iflag,3) + matrices%ta(kk)                 !age
     matrices%ages(iflag,4) = matrices%ages(iflag,4) + (temp_age(2)-temp_age(1))/dz      !dT/dz,z=0
     matrices%ages(iflag,5) = matrices%ages(iflag,5) + (temp_age(i)-temp_age(i-1))/dz    !dT/dz,z=z_c
     matrices%zc(kk)=dz*float(i-1)+xjunk
      write(125,*) Tc,dz*float(i-1)+xjunk,(temp_age(2)-temp_age(1))/dz
    exit
   endif
enddo

!ending loop through samples
enddo

!average ages matrix for each system, ie, obtain an average closure depth for AFT
do i=1,7
  matrices%ages(i,:)=matrices%ages(i,:)/float(matrices%nsystems(i))
enddo

open(355,file=params%run//"/stuff/closures.txt",status="unknown")
do j=1,7
  write(355,*) matrices%ages(j,1),matrices%ages(j,2),matrices%ages(j,3),matrices%ages(j,4),matrices%ages(j,5)
enddo

deallocate(temp,temp2)
deallocate(diag,sup,inf,f,temp_pr,temp_age,tdot)
deallocate(closure,query)

return
end subroutine
