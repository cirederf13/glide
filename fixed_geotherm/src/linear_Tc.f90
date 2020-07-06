      subroutine linear_tc (cooling,temp,closure,mz,iflag)

      implicit double precision (a-h,o-z)

      double precision cooling(mz),temp(mz),closure(mz)
      double precision,dimension(:),allocatable::log_cooling
      double precision secinyr,Tc_10,E,A,r
      
      integer iflag,mz
     
      allocate(log_cooling(mz))

      energy=energy_cal*4.184*1.e3

!these values are taken from reiners' review paper and correspond to Ea (x1000.) and omega (x secinyr)      
      if (iflag.eq.1) then
!here are Ea and D0/a2 values for AFT from ketcham. 1999
!taken from reiners 2004
      energy=147.d3
      geom=1.d0
      diff=2.05e6 
      Tc_10=116.   
      elseif (iflag.eq.2) then
!here are Ea and D0/a2 values for ZFT from brandon1999
!taken from reiners 2004
      energy=208.d3
      geom=1.d0
      diff=1.0e8
      Tc_10=232.
      elseif (iflag.eq.3) then
!here are Ea and D0/a2 values for AHe from Farley et al. 2000
!taken from reiners 2004
      energy=138.d3
      geom=1.d0
      diff=7.64e7
      Tc_10=67.
      elseif (iflag.eq.4) then
!here are Ea and D0/a2 values for ZHe from reiners2004
!taken from reiners 2004
      energy=169.d3
      geom=1.d0
      diff=7.03e5
      Tc_10=183.
      
      !the following are for argon argon, might be a bit much for most, 51,52,53 in glide
      elseif (iflag.eq.5) then
!here are Ea and D0/a2 values for hbl from harrison81
!taken from reiners 2004
      energy=268.d3
      geom=1.d0
      diff=1320.
      Tc_10=553.
      elseif (iflag.eq.6) then
!here are Ea and D0/a2 values for mus from hames&bowring1994,robbins72
!taken from reiners 2004
      energy=180.d3
      geom=1.d0
      diff=3.91
      Tc_10=380.
      elseif (iflag.eq.7) then
!here are Ea and D0/a2 values for bio from grove&harrison1996
!taken from reiners 2004
      energy=197.d3
      geom=1.d0
      diff=733.
      Tc_10=348.
      endif
      
      Tc_10=Tc_10+273.16
      
      r=8.314
      secinyr=3600.*24.*365.25e6  
      diff=diff*secinyr
      E=energy/r
      A=diff*r/energy   
      A=log(A)
      log_cooling=log(cooling)
      
      closure=(E+2.*Tc_10)/(A+(2.*log(Tc_10))+2.-log_cooling)
      
      closure=closure-273.16
      
      return
      end

