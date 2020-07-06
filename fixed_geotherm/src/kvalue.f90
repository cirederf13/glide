      subroutine kvalue(i,j,nx,ny,dkx,dky,kx,ky)
!  Subroutine KVALUE finds the wavenumber coordinates of one 
!  element of a rectangular grid from subroutine FOURN.  
!
!  Input parameters:
!    i  - index in the ky direction.
!    j  - index in the kx direction.
!    nx - dimension of grid in ky direction (a power of two).
!    ny - dimension of grid in kx direction (a power of two).
!    dkx - sample interval in the kx direction.
!    dky - sample interval in the ky direction.
!
!  Output parameters:
!    kx - the wavenumber coordinate in the kx direction.
!    ky - the wavenumber coordinate in the ky direction.
!
      double precision kx,ky,dkx,dky
      integer i,j,nx,ny,nxqx,nyqy
      nxqx=nx/2+1
      nyqy=ny/2+1
      !print*,"kvalue,value",i,j
      if(j.le.nxqx)then
           kx=(j-1)*dkx
         else
           kx=(j-nx-1)*dkx
      end if
      if(i.le.nyqy)then
           ky=(i-1)*dky
         else
            ky=(i-ny-1)*dky
          end if
      return
      end
