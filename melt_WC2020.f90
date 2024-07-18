! File containing the melt subroutine from HAILCAST, modified
! to be variable by height

!f2py -c melt.f90 -m melt

SUBROUTINE MELT_WC2020short(DIN, TLAYER, PLAYER, RLAYER, VT, SEKDEL, DENSE, DENSA, DOUT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  This is a spherical hail melting estimate based on Wang and Chueh 
  !!!  2020 (AR). Assuming a short-lobed sphere with a lobe factor of 6.
  !!!
  !!!  INPUT:  DIN      diameter prior to melting 
  !!!          SEKDEL   timstep (s)
  !!!          TLAYER   environmental temperature (K)
  !!!          PLAYER   environmental pressure (Pa)
  !!!          RLAYER   environmental water vapor mixing ratio (kg/kg)
  !!!          VT       terminal velocity of stone (m/s)
  !!!          DENSE    density of hail (kg/m3)
  !!!          DENSA    density of air (kg/m3)
  !!!          LDEPTH   depth of the layer (m)
  !!!  IN/OUTPUT: DOUT     diameter after melting (m)
  !!!          
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

      REAL*8 DIN, SEKDEL, DENSA, DENSE
!f2py REAL*8, INTENT(IN) :: DIN, SEKDEL, DENSE, DENSA
      REAL*8 DOUT
!f2py REAL*8, INTENT(OUT) :: DOUT
      REAL*8 TLAYER, PLAYER, RLAYER, VT
!f2py REAL*8, INTENT(IN) :: TLAYER, PLAYER, RLAYER, VT

      !REAL*8 delta, ewet, de, der, wetold, wetbulb, wetbulbk
      REAL*8 tclayer
      REAL*8 D, DCM, r, eps, dmlt, GMASS, RE, massorg, eenv,  &
              pi, DI, g, AK, ANU, E, H, AH, AE, &
              rv, ALF, ALV, t0, newmass
      REAL*8 RATIO, ESAT, ALV_SI, ALS_SI, RHOKOR, eenv_si, &
              TS_SI, DELRW, DELRWC, RHOAIR  !suffix _SI indicates is a duplicate value, but in SI units
      
      ! Set constantsi
      ALF = 79.7 ! latent heat of melting (cal/g)
      ALV = 597.3  ! latent of evaporation (cal/g)
      eps = 0.622 ! molecular weight ration of water to dry air (unitless)
      g = 9.81 !m/s-2
      t0 = 273.155 ! temp of ice/water melting interface (K)
      pi = 3.1415927
      rv = 461.48 !water vapor gas constant (J K-1 kg-1)

      ! Convert diameter to cm
      D = DIN !m
      DCM = DIN * 100 ! cm
      r = D/2. ! radius of stone (m)

      !SI unit constants (for VAPORCLOSE)
      ALV_SI = 2500000. !J/kg
      ALS_SI = 2836050. !J/kg
      TS_SI = 273.155 !K

      !Mass of stone in kg
      GMASS = (DENSE * PI * (D**3.)) / 6. !kg m-3 * m3 = kg

      !Calculate partial vapor pressure
      eenv = (PLAYER*RLAYER) / (RLAYER+eps)
      eenv = eenv / 100.  !convert to mb = hPa

      !!!  CALCULATE THE DIFFUSIVITY DI (cm2/s)
      !D0=0.226*1.E-4  ! change to m2/s, not cm2/s
      !DI=D0*(tclayer/273.155)**1.81*(100000./PLAYER) !(cm2/s)
      !RAS 20230414 - changed to Pruppacher and Klett eq. (13-3) so 
      ! units remain entirely cm2/s
      !tclayer = TLAYER - 273.155 ! Kelvin to C
      !DI = 0.211*(tclayer/273.155)**1.94*(1013.25/PLAYER)
      !RAS 20231215 - surprisingly, wants the temperature still in K,
      ! but pressure in hPa, and produces results in cm2/s. Sigh.
      DI = 0.211*(TLAYER/273.155)**1.94*(1013.25/(PLAYER*1E-2))

      !Calculate thermal conductivity
      ! Per Table A1 of RH87, AK wants T in C.
      tclayer = TLAYER - 273.155 ! Kelvin to C
      AK=(5.8+0.0184*tclayer)*1.E-5  !cal/(cm*sec*C)  
        
      !Dynamic viscosity - wants T in K, NOT C. Rogers and Yau pg 102.
      !ANU = (0.00001718)*(273.155+120.)/(tclayer+120.)*(tclayer/273.155)**(1.5) !kg m-1 s-1
      ANU = (0.00001718)*(273.155+120.)/(TLAYER+120.)*(TLAYER/273.155)**(1.5) !output kg m-1 s-1
      

      !!!  CALCULATE THE VENTILATION COEFFICIENT - NEEDED FOR GROWTH FROM VAPOR
      !The coefficients in the ventilation coefficient equations have been
      !experimentally derived, and are expecting cal-C-g units.  Do some conversions.

      !!!!  CALCULATE THE REYNOLDS NUMBER - unitless
      !RAS bug fix 20230414
      RE=D*VT*DENSA/ANU  !D, VT, DENSA, ANU all in SI units; RE now unitless 

      ! Calculate ventilation coefficents
      H=(0.71)**(0.333333333)*(RE**0.50) !ventilation coefficient heat (fh)
      E=(0.60)**(0.333333333)*(RE**0.50) !ventilation coefficient vapor (fv)

      !!!   SELECT APPROPRIATE VALUES OF AH AND AE ACCORDING TO Reynolds number
      AH = 0.
      AE = 0.
      !Assumes we can treat all Re roughly equally
      AH=0.002953*H**2.-0.2583*H+11.13
      AE=0.002953*E**2.-0.2583*E+11.13

      !Calculate the difference in vapor density in environment and hailstone surface
      !Set some SI unit constants, since VAPORCLOSE can handle (and uses) SI units
      RATIO = 1./273.155
      !!! First, the hailstone.
      !Assume stone is always in "wet growth" (i.e., hailstone surface is 0C and we 
      ! should be worried about saturation wrt liquid water)
      !Also assume hailstone temperature is 273.155 (see contstants set above)
      ESAT=611.*EXP(ALV_SI/rv*(RATIO-1./TS_SI))  !From (2.11) of RY89. Units Pa
      RHOKOR=ESAT/(rv*TS_SI)  !Pa /(J K-1 kg-1 K) = Pa kg / J =  kg/m3
      
      !!!  NOW FOR THE AMBIENT/IN-CLOUD CONDITIONS 
      !It is not necessarily saturated!! 
      !Convert previously calculated eenv to SI units:
      eenv_si = eenv * 100. !Pa
      !Now, environmental water vapor density using RY89 (2.1)
      RHOAIR = eenv_si / (rv * TLAYER) !Pa / (J K-1 kg-1 K) --> !Pa / (J kg-1)


      !!!  CALC THE DIFFERENCE(KG/M3): <0 FOR CONDENSATION, 
      !!!  >0 FOR EVAPORATION
      !DELRW=(RHOKOR-RHOOMG)  !units kg/m3
      DELRW=(RHOKOR-RHOAIR)  !units kg/m3
      !! Bug fix - should actually be vapor density in air minus that in stone
      DELRW=(RHOAIR-RHOKOR)

      ! Convert to non-SI units
      DELRWC = DELRW * (1.E3) * (1.E-6)  !g/cm3


      !Calculate new mass growth.
      !tclayer will be equivalent to T_env (K) - T_sfc (K), where T_sfc = 273.155
      tclayer = TLAYER - 273.155 ! Kelvin to C

      dmlt = -SEKDEL  / ALF * (2.*PI*DCM*AH*AK*tclayer + 2.*PI*DCM*AE*ALV*DI*DELRWC) ! gram
      dmlt = dmlt/1000 ! convert g to kg
      IF (dmlt.gt.0.) dmlt = 0

      !Find the new hailstone diameter
      massorg = 1.33333333*pi*r*r*r*DENSE
      newmass = massorg + dmlt
      if (newmass.lt.0.0) newmass = 0.0
      D = 2.*(0.75*newmass/(pi*DENSE))**0.333333333
      
      DOUT = D

  END SUBROUTINE MELT_WC2020short

  SUBROUTINE MELT_WC2020long(DIN, TLAYER, PLAYER, RLAYER, VT, SEKDEL, DENSE, DENSA, DOUT)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  This is a spherical hail melting estimate based on Wang and Chueh 
    !!!  2020 (AR). Assuming a long-lobed sphere with a lobe factor of 6.
    !!!
    !!!  INPUT:  DIN      diameter prior to melting 
    !!!          SEKDEL   timstep (s)
    !!!          TLAYER   environmental temperature (K)
    !!!          PLAYER   environmental pressure (Pa)
    !!!          RLAYER   environmental water vapor mixing ratio (kg/kg)
    !!!          VT       terminal velocity of stone (m/s)
    !!!          DENSE    density of hail (kg/m3)
    !!!          DENSA    density of air (kg/m3)
    !!!          LDEPTH   depth of the layer (m)
    !!!  IN/OUTPUT: DOUT     diameter after melting (m)
    !!!          
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
  
        REAL*8 DIN, SEKDEL, DENSA, DENSE
  !f2py REAL*8, INTENT(IN) :: DIN, SEKDEL, DENSE, DENSA
        REAL*8 DOUT
  !f2py REAL*8, INTENT(OUT) :: DOUT
        REAL*8 TLAYER, PLAYER, RLAYER, VT
  !f2py REAL*8, INTENT(IN) :: TLAYER, PLAYER, RLAYER, VT
  
        !REAL*8 delta, ewet, de, der, wetold, wetbulb, wetbulbk
        REAL*8 tclayer
        REAL*8 D, DCM, r, eps, dmlt, GMASS, RE, massorg, eenv,  &
                pi, DI, g, AK, ANU, E, H, AH, AE, &
                rv, ALF, ALV, t0, newmass
        REAL*8 RATIO, ESAT, ALV_SI, ALS_SI, RHOKOR, eenv_si, &
                TS_SI, DELRW, DELRWC, RHOAIR  !suffix _SI indicates is a duplicate value, but in SI units
        
        ! Set constantsi
        ALF = 79.7 ! latent heat of melting (cal/g)
        ALV = 597.3  ! latent of evaporation (cal/g)
        eps = 0.622 ! molecular weight ration of water to dry air (unitless)
        g = 9.81 !m/s-2
        t0 = 273.155 ! temp of ice/water melting interface (K)
        pi = 3.1415927
        rv = 461.48 !water vapor gas constant (J K-1 kg-1)
  
        ! Convert diameter to cm
        D = DIN !m
        DCM = DIN * 100 ! cm
        r = D/2. ! radius of stone (m)
  
        !SI unit constants (for VAPORCLOSE)
        ALV_SI = 2500000. !J/kg
        ALS_SI = 2836050. !J/kg
        TS_SI = 273.155 !K
  
        !Mass of stone in kg
        GMASS = (DENSE * PI * (D**3.)) / 6. !kg m-3 * m3 = kg
  
        !Calculate partial vapor pressure
        eenv = (PLAYER*RLAYER) / (RLAYER+eps)
        eenv = eenv / 100.  !convert to mb = hPa
  
        !!!  CALCULATE THE DIFFUSIVITY DI (cm2/s)
        !D0=0.226*1.E-4  ! change to m2/s, not cm2/s
        !DI=D0*(tclayer/273.155)**1.81*(100000./PLAYER) !(cm2/s)
        !RAS 20230414 - changed to Pruppacher and Klett eq. (13-3) so 
        ! units remain entirely cm2/s
        !tclayer = TLAYER - 273.155 ! Kelvin to C
        !DI = 0.211*(tclayer/273.155)**1.94*(1013.25/PLAYER)
        !RAS 20231215 - surprisingly, wants the temperature still in K,
        ! but pressure in hPa, and produces results in cm2/s. Sigh.
        DI = 0.211*(TLAYER/273.155)**1.94*(1013.25/(PLAYER*1E-2))
  
        !Calculate thermal conductivity
        ! Per Table A1 of RH87, AK wants T in C.
        tclayer = TLAYER - 273.155 ! Kelvin to C
        AK=(5.8+0.0184*tclayer)*1.E-5  !cal/(cm*sec*C)  
          
        !Dynamic viscosity - wants T in K, NOT C. Rogers and Yau pg 102.
        !ANU = (0.00001718)*(273.155+120.)/(tclayer+120.)*(tclayer/273.155)**(1.5) !kg m-1 s-1
        ANU = (0.00001718)*(273.155+120.)/(TLAYER+120.)*(TLAYER/273.155)**(1.5) !output kg m-1 s-1
        
  
        !!!  CALCULATE THE VENTILATION COEFFICIENT - NEEDED FOR GROWTH FROM VAPOR
        !The coefficients in the ventilation coefficient equations have been
        !experimentally derived, and are expecting cal-C-g units.  Do some conversions.
  
        !!!!  CALCULATE THE REYNOLDS NUMBER - unitless
        !RAS bug fix 20230414
        RE=D*VT*DENSA/ANU  !D, VT, DENSA, ANU all in SI units; RE now unitless 
  
        ! Calculate ventilation coefficents
        H=(0.71)**(0.333333333)*(RE**0.50) !ventilation coefficient heat (fh)
        E=(0.60)**(0.333333333)*(RE**0.50) !ventilation coefficient vapor (fv)
  
        !!!   SELECT APPROPRIATE VALUES OF AH AND AE ACCORDING TO Reynolds number
        AH = 0.
        AE = 0.
        !Assumes we can treat all Re roughly equally
        AH = -4.589E-6*H**3.0 + 5.47E-3*H**2 - 0.2463*H + 9.533
        AE = -4.589E-6*E**3.0 + 5.47E-3*E**2 - 0.2463*E + 9.533
  
        !Calculate the difference in vapor density in environment and hailstone surface
        !Set some SI unit constants, since VAPORCLOSE can handle (and uses) SI units
        RATIO = 1./273.155
        !!! First, the hailstone.
        !Assume stone is always in "wet growth" (i.e., hailstone surface is 0C and we 
        ! should be worried about saturation wrt liquid water)
        !Also assume hailstone temperature is 273.155 (see contstants set above)
        ESAT=611.*EXP(ALV_SI/rv*(RATIO-1./TS_SI))  !From (2.11) of RY89. Units Pa
        RHOKOR=ESAT/(rv*TS_SI)  !Pa /(J K-1 kg-1 K) = Pa kg / J =  kg/m3
        
        !!!  NOW FOR THE AMBIENT/IN-CLOUD CONDITIONS 
        !It is not necessarily saturated!! 
        !Convert previously calculated eenv to SI units:
        eenv_si = eenv * 100. !Pa
        !Now, environmental water vapor density using RY89 (2.1)
        RHOAIR = eenv_si / (rv * TLAYER) !Pa / (J K-1 kg-1 K) --> !Pa / (J kg-1)
  
  
        !!!  CALC THE DIFFERENCE(KG/M3): <0 FOR CONDENSATION, 
        !!!  >0 FOR EVAPORATION
        !DELRW=(RHOKOR-RHOOMG)  !units kg/m3
        DELRW=(RHOKOR-RHOAIR)  !units kg/m3
        !! Bug fix - should actually be vapor density in air minus that in stone
        DELRW=(RHOAIR-RHOKOR)
  
        ! Convert to non-SI units
        DELRWC = DELRW * (1.E3) * (1.E-6)  !g/cm3
  
  
        !Calculate new mass growth.
        !tclayer will be equivalent to T_env (K) - T_sfc (K), where T_sfc = 273.155
        tclayer = TLAYER - 273.155 ! Kelvin to C
  
        dmlt = -SEKDEL  / ALF * (2.*PI*DCM*AH*AK*tclayer + 2.*PI*DCM*AE*ALV*DI*DELRWC) ! gram
        dmlt = dmlt/1000 ! convert g to kg
        IF (dmlt.gt.0.) dmlt = 0
  
        !Find the new hailstone diameter
        massorg = 1.33333333*pi*r*r*r*DENSE
        newmass = massorg + dmlt
        if (newmass.lt.0.0) newmass = 0.0
        D = 2.*(0.75*newmass/(pi*DENSE))**0.333333333
        
        DOUT = D
  
    END SUBROUTINE MELT_WC2020long