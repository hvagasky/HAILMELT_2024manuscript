  SUBROUTINE TERMINL_RH87X_Re(DENSA,DENSE,D,TC,VT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!
  !!!! INTERP: Calculate terminal velocity of the hailstone
  !!!!
  !!!! INPUT: DENSA  density of updraft air (kg/m3)
  !!!!        DENSE  density of hailstone
  !!!!        D      diameter of hailstone (m)
  !!!!        TC     updraft temperature (K)
  !!!! OUTPUT:VT     hailstone terminal velocity (m/s)
  !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      
      REAL*8 D
!f2py REAL*8, INTENT(IN) :: D
      REAL*8 DENSA, DENSE, TC, VT
!f2py REAL, INTENT(IN) :: DENSA, DENSE, TC
!f2py REAL, INTENT(OUT) :: VT
      REAL*8 GMASS, GX, RE, W, Y
      REAL, PARAMETER :: PI = 3.141592654, G = 9.78956
      REAL*8 ANU
      
      !Mass of stone in kg
      GMASS = (DENSE * PI * (D**3.)) / 6.
      
      !Dynamic viscosity
      ANU = (0.00001718)*(273.155+120.)/(TC+120.)*(TC/273.155)**(1.5)
      
      !CALC THE BEST NUMBER, X AND REYNOLDS NUMBER, RE 
      GX=(8.0*GMASS*G*DENSA)/(PI*(ANU*ANU))
      !RE=(GX/0.6)**0.5

      !SELECT APPROPRIATE EQUATIONS FOR TERMINAL VELOCITY DEPENDING ON 
      !THE BEST NUMBER (Rasmussen and heymsfield, 1989: equations B1-B4)
      IF (GX.LT.550) THEN
        W=LOG10(GX)
        Y= -1.7095 + 1.33438*W - 0.11591*(W**2.0)      
        RE=10**Y
        !VT=ANU*RE/(D*DENSA)
      ELSE IF (GX.GE.550.AND.GX.LT.1800) THEN
        W=LOG10(GX)
        Y= -1.81391 + 1.34671*W - 0.12427*(W**2.0) + 0.0063*(W**3.0)
        RE=10**Y
        !VT=ANU*RE/(D*DENSA)
      ELSE IF (GX.GE.1800.AND.GX.LT.3.45E08) THEN
        RE=0.4487*(GX**0.5536)
        !VT=ANU*RE/(D*DENSA)
      ELSE 
        RE=(GX/0.6)**0.5
        !VT=ANU*RE/(D*DENSA)
      ENDIF

      ! Calculate vertical velocity
      VT=ANU*RE/(D*DENSA)
      
  END SUBROUTINE TERMINL_RH87X_Re
