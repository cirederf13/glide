      SUBROUTINE FOUR1(DATA,NN,ISIGN)
!------------------------------------
!--- computes Discrete Complex Fast Fourier Transform
!--- taken from 'Numerical Recipes', Press et al., 1992
!    see also N/R p. 501.
!--- routine called: none
      double precision WR,WI,WPR,WPI,WTEMP,THETA
      double precision DATA(2*NN)
      double precision TEMPR,TEMPI
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
 1      IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
 11   CONTINUE
      MMAX=2
 2    IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959/(ISIGN*MMAX)
        WPR=-2.*SIN(0.5*THETA)**2
        WPI=SIN(THETA)
        WR=1.0
        WI=0.0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
 12       CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
 13     CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END SUBROUTINE FOUR1
      
