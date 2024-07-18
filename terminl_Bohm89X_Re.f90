! File containing the melt subroutine from HAILCAST, modified
! to be variable by height

!f2py -c terminl_Bohm89X_Re.f90 -m terminl_Bohm89X_Re

SUBROUTINE TERMINL_Bohm89X_Re(DENSA,DENSE,D,TC,VT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!
  !!!! INTERP: Calculate terminal velocity of the hailstone using (11)
  !!!!  from Bohm 1989 JAS
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
      REAL*8 DENSA, DENSE, TC
!f2py REAL, INTENT(IN) :: DENSA, DENSE, TC
      REAL*8 VT
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

      !Calculate Re given the Best number GX using Bohm's eq. 11:
      RE = 8.5 * ( (1+0.1519*GX**0.5)**0.5 -1. )**2.

      !Calculate vertical velocity
      VT=ANU*RE/(D*DENSA)

  END SUBROUTINE TERMINL_Bohm89X_Re
