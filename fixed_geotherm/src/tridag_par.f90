      SUBROUTINE tridag_par(a,b,c,r,u,n,chunk)
      
      integer n,chunk
      double precision a(n,chunk),b(n,chunk),c(n,chunk),r(n,chunk),u(n,chunk)
      integer j 
      double precision,dimension(:),allocatable::bet
      double precision,dimension(:,:),allocatable::gam
      !print*,"in tridag"
      !pause
      allocate (gam(n,chunk),bet(chunk))
      !if(b(1,:).eq.0.) stop 'in tridag'
      !print*,"in tridag"
      bet=b(1,:)
      u(1,:)=r(1,:)/bet(:)
      do 11 j=2,n
        gam(j,:)=c(j-1,:)/bet(:)
        bet(:)=b(j,:)-a(j,:)*gam(j,:)
       ! if(bet.eq.0.) then
       ! this bit should be reinserted
       ! print*,'tridag failed'
       ! stop
       ! endif
        u(j,:)=(r(j,:)-a(j,:)*u(j-1,:))/bet(:)
11    continue
      do 12 j=n-1,1,-1
        u(j,:)=u(j,:)-gam(j+1,:)*u(j+1,:)
12    continue
      deallocate (gam,bet)
      return
      END

