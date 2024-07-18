  SUBROUTINE TERMINL_HW2014X_Re(DENSA,DENSE,D,TC,VT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!
  !!!! INTERP: Calculate terminal velocity of the hailstone using the
  !!!!  (6a) and (6b) from Heymsfield and Wright (2014)
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
      REAL*8 GMASS, GX, RE
      REAL, PARAMETER :: PI = 3.141592654, G = 9.78956
      REAL*8 ANU
      
      !Mass of stone in kg
      GMASS = (DENSE * PI * (D**3.)) / 6.
      
      !Dynamic viscosity
      ANU = (0.00001718)*(273.155+120.)/(TC+120.)*(TC/273.155)**(1.5)
      
      !CALC THE BEST NUMBER (GX)
      GX=(8.0*GMASS*G*DENSA)/(PI*(ANU*ANU))

      !Calculate Re given the Best number GX using HW2014's eq. 6a and 6b:
      IF (GX.LT.6.77E4) THEN
        RE=0.106*GX**0.693 !just graupel
      ELSE
        RE=0.55*GX**0.545 !graupel and hail
      ENDIF

      !Calculate vertical velocity
      VT=ANU*RE/(D*DENSA)
      
  END SUBROUTINE TERMINL_HW2014X_Re
