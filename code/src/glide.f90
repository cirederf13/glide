program glide

! Gaussian Linear Inversion of Data to Exhumation rate
! Matthew Fox and Frédéric Herman, 2009-2012
! no warranty on the code.
! the code cannot be used for publication without consent from the authors

use definitions

implicit none

! definition of each type are found in module_definitions.f90

type (matr) matrices
type (parm) params

integer iter,post,i
double precision,dimension(:),allocatable::mis

!a priori and post chi**2
allocate(mis(2))

mis=0.

print*,"initializing parameters" 
call initialize_parameters(matrices,params)

print*,"building matrices" 
call build_matrices(matrices,params)
   
!set exhumation rates to prior model
matrices%edot = params%edot_mean
matrices%edot_dat = matrices%edot_pr      

print*,"calculating misfit with prior model" 
call syn_ages (matrices,params,mis(1))
print*,"prior misfit=",mis(1)

!post defines which posterior matrices to calclate, see posterior_dat   
post=2

iter=0

do while (iter.lt.3)

   iter=iter+1
   print*,"solving inverse, iteration =",iter
   call solve_inverse(matrices,params)
   
   print*,"computing closure depths"
   call find_zc(matrices,params)

enddo
   
print*,"calculating misfit with posterior model" 
call syn_ages(matrices,params,mis(2))

print*,"calculating posteriori matrix for control points"
call posterior(matrices,params)
print*,"calculating FULL posteriori matrix for control points"
call posterior_full(matrices,params)

print*,"calculating posteriori matrix for data points"
call posterior_dat(matrices,params,post)

print*,"writing output"
call write_output(matrices,params) 

print*,"a priori chi=",mis(1),"post chi=",mis(2)

print*,"deallocating arrays"
call clean_up(matrices,params)

deallocate(mis)

end program GLIDE

