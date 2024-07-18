  !SUBROUTINE TERMINL_HEYMSFIELD2018(D,P,VT)
  SUBROUTINE TERMINL_H2018V_D(D,P,VT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!
  !!!! INTERP: Calculate terminal velocity of the hailstone
  !!!!
  !!!! INPUT: D      diameter of hailstone (m)
  !!!!        P      atmospheric pressure (Pa)
  !!!! OUTPUT:VT     hailstone terminal velocity (m/s)
  !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      
      REAL*8 D, P
!f2py REAL*8, INTENT(IN) :: D, P
      REAL*8 VT
!f2py REAL*8, INTENT(OUT) :: VT
      
      !Calculate terminal velocity at 1000 hPa (m s-1)
      D = D * 100 ! Requires cm
      IF (D .lt. 0.5) THEN ! Graupel
              VT = 6.35 * D**0.87
      ELSEIF ( D .ge. 0.5) THEN ! Hail
              VT = 6.1 * D**0.72
      ENDIF

      !Convert from assumed 1000 hPa pressure to actual pressure
      P = P / 100
      VT = VT * (1000/P)**0.55

  END SUBROUTINE TERMINL_H2018V_D  
