  SUBROUTINE TERMINL_H2018X_Re(DENSA,DENSE,D,TC,VT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!
  !!!! INTERP: Calculate terminal velocity of the hailstone using the
  !!!!  (5a) and (5b) from Heymsfield et al. 2018)
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

      !Calculate Re given the Best number GX using H2018's eq. 5a and 5b:
      IF (GX.LT.1.E6) THEN
        RE=0.29*GX**0.59 !just graupel
      ELSE
        RE=1.16*GX**0.49 !graupel and hail
      ENDIF

      !Calculate vertical velocity
      VT=ANU*RE/(D*DENSA)
      
  END SUBROUTINE TERMINL_H2018X_Re
