      subroutine closure_temps (cooling,temp,closure,mz,iflag)

      !calculates closure temperatures at for the material derviatives t_dot 

      implicit double precision (a-h,o-z)

      double precision cooling(mz),temp(mz),closure(mz)
      double precision,dimension(:),allocatable::tau
      integer iflag
     
      allocate(tau(mz))


!these values are taken from reiners' review paper and correspond to Ea (x1000.) and omega (x secinyr)      
      if (iflag.eq.1) then
!here are Ea and D0/a2 values for AFT from ketcham. 1999
!taken from reiners 2004
      energy=147.d3
      geom=1.d0
      diff=2.05e6*3600.*24.*365.25e6      
      elseif (iflag.eq.2) then
!here are Ea and D0/a2 values for ZFT from reiners2004
!taken from reiners 2004
      energy=208.d3!/(4.184*1.e3)
      geom=1.d0
      diff=4.0e8*3600.*24.*365.25e6
      !energy=224.d3
      !geom=1.d0
      !diff=1.24e8*3600.*24.*365.25e6
diff=1.24e8*3600.*24.*365.25e6
      elseif (iflag.eq.3) then
!here are Ea and D0/a2 values for AHe from Farley et al. 2000
!taken from reiners 2004
      energy=138.d3
      geom=1.d0
      diff=7.64e7*3600.*24.*365.25e6
      elseif (iflag.eq.4) then
!here are Ea and D0/a2 values for ZHe from reiners2004
!taken from reiners 2004
      energy=169.d3!/(4.184*1.e3)
      geom=1.d0
      diff=7.03e5*3600.*24.*365.25e6
      !energy=178.d3
      !geom=1.d0
      !diff=7.03d5*3600.d0*24.d0*365.25d6
      
      !the following are for argon argon, might be a bit much for most, 51,52,53 in glide
      elseif (iflag.eq.5) then
!here are Ea and D0/a2 values for hbl from harrison81
!taken from reiners 2004
      energy=268.d3
      geom=1.d0
      diff=1320*3600.*24.*365.25e6
      elseif (iflag.eq.6) then
!here are Ea and D0/a2 values for mus from hames&bowring1994,robbins72
!taken from reiners 2004
      energy=180.d3
      geom=1.d0
      diff=3.91*3600.*24.*365.25e6
      elseif (iflag.eq.7) then
!here are Ea and D0/a2 values for bio from grove&harrison1996
!taken from reiners 2004
      energy=197.d3
      geom=1.d0
      diff=733.*3600.*24.*365.25e6
      endif
   
      r=8.314d0
        ! note that cooling cannot be nil and is therefore forced to
! be at least 1deg/10My
        cooling=max(cooling,1.d0/10.d0)
        tau=r*(temp+273.d0)**2/energy/cooling
        closure=energy/r/log(geom*tau*diff)-273.d0

      return
      end

