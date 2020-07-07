subroutine clean_up(matrices,params)

!this just deallocates everything

use definitions

implicit none

! definition of each type are found in module_definitions.f90

type (parm) params
type (matr) matrices

!the purpose of this routine is to deallocate everything

!allocated in initialise parameters
deallocate(matrices%ta,matrices%a_error,matrices%zc,matrices%zcp,matrices%x,&
           matrices%x_dum,matrices%y_dum,matrices%elev,matrices%misfits,matrices%syn_age,&
           matrices%isys,matrices%depth1,matrices%y,matrices%depth,matrices%zz,matrices%elev_true)
          
deallocate(matrices%x_true,matrices%y_true,matrices%x_dum_true,matrices%y_dum_true)

deallocate(matrices%tsteps,matrices%tsteps_sum,matrices%edot_pr,matrices%edot_pr_dum)

deallocate(matrices%nsystems,matrices%ages)

!allocated in build matrices
deallocate(matrices%A,matrices%edot,matrices%edot_dat,matrices%eps_dum,matrices%eps,&
           matrices%cov,matrices%cpost,matrices%eps_dat,matrices%eps_res,matrices%H,matrices%Y2)


!allocated in resolution
!deallocate(matrices%R,matrices%kernel)

return
end subroutine clean_up



