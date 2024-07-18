! File containing the melt subroutine from HAILCAST, modified
! to be variable by height

!f2py -c melt.f90 -m melt

SUBROUTINE MELT_G1969orig(DIN, TLAYER, PLAYER, RLAYER, LDEPTH, VT, DOUT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  This is a spherical hail melting estimate based on the Goyer 
  !!!  et al. (1969) eqn (3).  The depth of the warm layer, estimated 
  !!!  terminal velocity, and mean temperature of the warm layer are 
  !!!  used.  DRB.  11/17/2003.
  !!!
  !!!  Note from RAS: In my opinion, eq. (3) of Goyer et al. (1969) has
  !!!  error. They include vapor diffusion and evaporation impacts
  !!!  as separate processes, one an opposite sign than the other,
  !!!  which doesn't make sense to me. This code uses the eq. 
  !!!  as published, except as per the original coding by DRB it uses
  !!!  the wet bulb temperature when calculating the temperature 
  !!!  difference for the heat diffusion, instead of the air temp. 
  !!!  MELT_G1969fix uses the equation, with the air temp, but with the 
  !!!  erroneous extra term removed. RAS 20231218.
  !!!
  !!!  INPUT:  TLAYER   mean layer temperature (K)
  !!!          PLAYER   mean layer pressure (Pa)
  !!!          RLAYER   mean layer mixing ratio (kg/kg)
  !!!          VT       terminal velocity of stone (m/s)
  !!!          LDEPTH   depth of the layer (m)
  !!!  IN/OUTPUT: D     diameter (m)
  !!!          
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

      REAL*8 DIN, D
!f2py REAL*8, INTENT(IN) :: DIN, D
      REAL*8 DOUT
!f2py REAL*8, INTENT(OUT) :: DOUT
      REAL TLAYER, PLAYER, RLAYER, LDEPTH, VT
!f2py REAL, INTENT(IN) :: TLAYER, PLAYER, RLAYER, LDEPTH, VT
      REAL*8 eenv, delta, ewet, de, der, wetold, wetbulb, wetbulbk
      REAL*8 tdclayer, tclayer, eps, b, hplayer
      REAL*8 a
      REAL*8 ka, lf, lv, t0, dv, pi, rv, rhoice, &
           tres, delt, esenv, rhosenv, essfc, rhosfc, dsig
      REAL*8 dmdt, mass, massorg, newmass, gamma, r, rho, re
      INTEGER wcnt
      
      D = DIN
      
      
      !Convert temp to Celsius, calculate dewpoint in celsius
      eps = 0.622
      tclayer = TLAYER - 273.155  ! C
      a = 2.53E11  !Pa
      b = 5.42E3   !K
      tdclayer = b / LOG(a*eps / (rlayer*player))  !RY89 eq. 2.28
      hplayer = player / 100.
      
      !Calculate partial vapor pressure
      eenv = (player*rlayer) / (rlayer+eps)  !Pa
      eenv = eenv / 100.  !convert to mb
      
      !Estimate wet bulb temperature (C)
      gamma = 6.6E-4*player
      delta = (4098.0*eenv)/((tdclayer+237.7)*(tdclayer+237.7))
      wetbulb = ((gamma*tclayer)+(delta*tdclayer))/(gamma+delta)
      
      !Iterate to get exact wet bulb
      wcnt = 0
      DO WHILE (wcnt .lt. 11)
        ewet = 6.108*(exp((17.27*wetbulb)/(237.3 + wetbulb))) 
        de = (0.0006355*hplayer*(tclayer-wetbulb))-(ewet-eenv)
        der= (ewet*(.0091379024 - (6106.396/(273.155+wetbulb)**2))) &
             - (0.0006355*hplayer)
        wetold = wetbulb
        wetbulb = wetbulb - de/der
        wcnt = wcnt + 1
        IF ((abs(wetbulb-wetold)/wetbulb.gt.0.0001)) THEN
           EXIT
        ENDIF
      ENDDO
      
      wetbulbk = wetbulb + 273.155  !convert to K
      ka = .02 ! thermal conductivity of air. W m-1 K-1
      lf = 3.34e5 ! latent heat of melting/fusion J/kg
      lv = 2.5e6  ! latent heat of vaporization  J/Kg
      t0 = 273.155 ! temp of ice/water melting interface K
      dv = 0.25e-4 ! diffusivity of water vapor (m2/s)
      !RH87 has dv ~ 0.211 cm2/s ... *(1 m/100 cm)^2 .. 0.211E-4 m2/s
      pi = 3.1415927
      rv = 1004. - 287. ! gas constant for water vapor
      !! Not clear where this number comes from. RH89 pg. 12 gives 
      !! 461.5 J kg-1 K-1. Might be for wet bulb? Leaving alone here. RAS 20231218.
      rhoice = 917.0 ! density of ice (kg/m**3)
      r = D/2. ! radius of stone (m)
      
      !Compute residence time in warm layer
      tres = LDEPTH / VT
        
      !Calculate dmdt based on eqn (3) of Goyer et al. (1969)
      !Reynolds number...from pg 317 of Atmo Physics (Salby 1996)
      !Just use the density of air at 850 mb...close enough.
      !rho = 85000./(287.*TLAYER)
      rho = player/(287.*TLAYER)  !we have the layer pressure, just use it
      !re = rho*r*VT*.01/1.7e-5
      !units fix - r is now in meters, not mm. Plus we need D, not r. RAS 20230414
      re = rho*r*2*VT/1.7e-5

      !Temperature difference between environment and hailstone surface
      delt = wetbulb !- 0.0 !assume stone surface is at 0C
                            !wetbulb is in Celsius

      !Difference in vapor density of air stream and equil vapor
      !density at the sfc of the hailstone
      esenv = 610.8*(exp((17.27*wetbulb)/  &
               (237.3 + wetbulb))) ! es environment in Pa
      rhosenv = esenv/(rv*wetbulbk)
      essfc = 610.8*(exp((17.27*(t0-273.155))/  &
               (237.3 + (t0-273.155)))) ! es environment in Pa
      rhosfc = essfc/(rv*t0)
      dsig = rhosenv - rhosfc

      !Calculate new mass growth
      dmdt = (-1.7*pi*r*(re**0.5)/lf)*((ka*delt)+((lv-lf)*dv*dsig))
      IF (dmdt.gt.0.) dmdt = 0
      mass = dmdt*tres
      
      !Find the new hailstone diameter
      massorg = 1.33333333*pi*r*r*r*rhoice
      newmass = massorg + mass
      if (newmass.lt.0.0) newmass = 0.0
      D = 2.*(0.75*newmass/(pi*rhoice))**0.333333333
      
      
      DOUT = D
  END SUBROUTINE MELT_G1969orig
