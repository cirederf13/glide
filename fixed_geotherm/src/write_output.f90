subroutine write_output(matrices,params)

use definitions

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices

integer i,iter,j,itime,loc,k,run
double precision rest,tmp4,tmp1,tmp2,tmp3,rest2
character cs*4,ct*4
double precision,dimension(:),allocatable::rgb

! It creates outputs to vizualise the results

allocate(rgb(3))

iter=1

open(555,file=params%run//"/stuff/cont.xy",status="unknown")
open(55,file=params%run//"/stuff/edot2",status="unknown")
open(99,file=params%run//"/stuff/unc2",status="unknown")

k=130
109 format(7X,2F6.2,A,F8.2,A,F8.2,A)
do j=1,params%m_max
rest=0.
        rest=sum(matrices%tsteps(1:params%m_max-j+1))
        rest2=sum(matrices%tsteps(1:params%m_max-j)) 
    
        itime = nint(rest*10)
        write(cs,'(i4)') itime
        if (itime.lt.10) cs(1:3)='000'
        if (itime.lt.100) cs(1:2)='00'
        if (itime.lt.1000) cs(1:1)='0'

        if (iter.eq.1) then
        open(k,file=params%run//"/"//cs,status="unknown")
        open(k+1,file=params%run//"/unc/"//cs,status="unknown")
        open(15,file=params%run//"/"//cs//".txt",status="unknown")
           
        if (params%lon1.lt.0) then
          write(15,109) params%lon1+360,params%lat2,"0 20 0 1 BL ",rest," -",rest2," Ma"
        else
          write(15,109) params%lon1,params%lat2,"0 20 0 1 BL ",rest," -",rest2," Ma"
        endif
        endif

        !do i=1,params%n
        !   write(k,*) sngl(matrices%x(i)),sngl(matrices%y(i)),sngl(matrices%edot(j+(i-1)*params%m_max))
        !enddo
        
        do i=1,params%dummy
            write(k,*) sngl(matrices%x_dum_true(i)),sngl(matrices%y_dum_true(i)),sngl(matrices%edot(j+(i-1)*params%m_max))
            write(k+1,*) sngl(matrices%x_dum_true(i)),sngl(matrices%y_dum_true(i)),&
                sngl((matrices%eps_dum(j+(i-1)*params%m_max)**2)/params%sigma2),sngl(matrices%sf(j+(i-1)*params%m_max))
            !write(k+1,*) sngl(matrices%x_dum_true(i)),sngl(matrices%y_dum_true(i)),&
            !    ((matrices%eps_dum(j+(i-1)*params%m_max)**2)/params%sigma2),(matrices%sf(j+(i-1)*params%m_max))
        enddo
        !print*,"written a timestep"
k=k+2
enddo

!Output individual exhumation histories for control points

!Give each separate control point a different colour, read one in or output one
run=0
if (run.eq.1) then
open(433,file="colours.txt",status="old")
else
open(432,file="colours.txt",status="unknown")
endif

call random_seed()
do j=params%dummy-params%contr+1,params%dummy
  if (run.ne.1) then 
   call random_number(rgb)
   rgb=rgb*255
   write(432,*) rgb
  else
   read(433,*) rgb 
  endif
  write(55,'(A,4X)',advance='NO') ">" 
  write(55,'(A,I3.3,A,I3.3,A,I3.3)',advance='NO') "-Wthicker,",int(rgb(1)),"/",int(rgb(2)),"/",int(rgb(3))
  write(55,*)
  write(99,'(A,4X)',advance='NO') ">" 
  write(99,'(A,I3.3,A,I3.3,A,I3.3)',advance='NO') "-Wthicker,",int(rgb(1)),"/",int(rgb(2)),"/",int(rgb(3))
  write(99,*)
    do i=1,params%m_max!-1
    rest=sum(matrices%tsteps(1:params%m_max-i))       
       loc= i+(j-1)*params%m_max
       write(55,*) rest,matrices%edot(loc)
       write(99,*) rest,((matrices%eps_dum(loc))**2.)/params%sigma2
       if (i.lt.params%m_max) then 
           write(99,*) rest,((matrices%eps_dum(loc+1))**2.)/params%sigma2
           write(55,*) rest,matrices%edot(loc+1)
       endif
     enddo
write(555,'(A,4X)',advance='NO') ">"
write(555,'(A,I3.3,A,I3.3,A,I3.3)',advance='NO') "-G",int(rgb(1)),"/",int(rgb(2)),"/",int(rgb(3))
write(555,*)
write(555,*) matrices%x_dum_true(j),matrices%y_dum_true(j)
enddo

iter=1
if (iter.eq.1) then
open(k,file=params%run//"/stuff/params",status="unknown")
write(k,*) "edot_pr=",params%edot_mean,"?",sqrt(params%sigma2)
write(k,*) "thermal",params%kappa,params%hp,params%Tb
write(k,*) "time stuff",params%t_total,"dt",params%deltat
write(k,*) "correlation",params%xL,params%angle,params%aspect
endif

open(123,file=params%run//"/stuff/aposterior_zc.txt",status="unknown")
open(234,file=params%run//"/stuff/misfits.txt",status="unknown")
do i=1,params%n
write(123,*) matrices%x_true(i),matrices%y_true(i),matrices%zcp(i),&
             matrices%isys(i),matrices%elev_true(i)
write(234,*) sngl(matrices%ta(i)),sngl(matrices%a_error(i)),&
             sngl(matrices%syn_age(i)),&
             matrices%isys(i),sngl(matrices%elev_true(i))
enddo
close(123)

return 
end subroutine

