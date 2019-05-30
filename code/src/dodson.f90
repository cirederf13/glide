      subroutine dodson (temp,time,nstep,age,closure,cooling,iflag)

      implicit double precision (a-h,o-z)

      double precision temp(nstep),time(nstep)
      integer iflag

!these values are taken from reiners' review paper and correspond to his Ea (x1000.) and omega (x secinyr)      
      if (iflag.eq.1) then
      energy=147.d3
      geom=1.d0
      diff=2.05e6*3600.*24.*365.25e6      
      elseif (iflag.eq.2) then
      energy=208.d3
      geom=1.d0
      diff=1.24e8*3600.*24.*365.25e6
      elseif (iflag.eq.3) then
      energy=138.d3
      geom=1.d0
      diff=7.64e7*3600.*24.*365.25e6
      elseif (iflag.eq.4) then
      energy=169.d3
      geom=1.d0
      diff=7.03e5*3600.*24.*365.25e6
      elseif (iflag.eq.5) then
      energy=268.d3
      geom=1.d0
      diff=1320*3600.*24.*365.25e6
      elseif (iflag.eq.6) then
      energy=180.d3
      geom=1.d0
      diff=3.91*3600.*24.*365.25e6
      elseif (iflag.eq.7) then
      energy=197.d3
      geom=1.d0
      diff=733.*3600.*24.*365.25e6
      endif
   
      r=8.314
      age=time(1)
        do i=nstep,2,-1
          if (i.eq.1) then
          cooling=(temp(i+1)-temp(i))/(time(i+1)-time(i))
          elseif (i.eq.nstep) then
          cooling=(temp(i)-temp(i-1))/(time(i)-time(i-1))
          else
          cooling=(temp(i+1)-temp(i-1))/(time(i+1)-time(i-1))
          endif
! note that cooling cannot be nil and is therefore forced to
! be at least 1deg/10My
        cooling=max(cooling,1.d0/10.d0)
        tau=r*(temp(i)+273.d0)**2/energy/cooling
        closure=energy/r/log(geom*tau*diff)-273.d0
          if (temp(i).gt.closure) then
          ratio=(closurep-tempp)/(closurep-tempp+temp(i)-closure)
          age=time(i)+(time(i-1)-time(i))*ratio
          return
          endif
        closurep=closure
        tempp=temp(i)
        enddo

      return
      end

