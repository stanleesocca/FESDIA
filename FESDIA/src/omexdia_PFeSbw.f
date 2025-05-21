!==========================================================================
! THE OMEXDIA model with P, Fe and S, implemented in FORTRAN
! version with dynamic water concentrations
!
! Karline Soetaert, nioz-yerseke
! Modules and transport functions are in file FEScommon.f
!==========================================================================

!==========================================================================
! initialisation = same as fesdiamod - see file fesdia_P.f
!==========================================================================

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the fesdia model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================

      SUBROUTINE fesdiamodbw (neq, t, Conc, dConc, yout, ip)
      use commonFESDIA
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i

      DOUBLE PRECISION  :: t,Conc(27*NP1),dConc(27*NP1),yout(*)
      DOUBLE PRECISION ::  FdetBW, SdetBW, O2BW, NO3BW, NO2bw, NH3bw,            &
     &                PO4bw,FePBW, CaPBW, DICbw, FeBW, FeOH3bw, H2Sbw,           &
     &                SO4bw, CH4bw, Padsbw, ALKbw, FeOH3BBW, Mnbw, 		    &
     &		        MnO2bw, MnO2Bbw, FeSbw, FeS2bw, S0bw, MnCO3bw, FeCO3BW,    &
     &                MnSBW
      DOUBLE PRECISION :: dFdetBW,dSdetBW,dO2BW,dNO3BW,dNO2BW,dNH3BW,            &
     &                dDICbw,dFeBW,dFeOH3bw,dH2Sbw,dSO4bw,dCH4bw,                &
     &                dPO4BW,dFePBW,dCaPBW,dPadsbw,dALKbw,dFeOH3BBW,		    &
     &		        dMnbw, dMnO2bw, dMnO2Bbw, dFeSbw, dFeS2bw, dS0bw,          &
     &                dMnCO3bw, dFeCO3BW, dMnSBW   

      DOUBLE PRECISION :: pF, FDETdepo, SDETdepo, FePdepo, CaPdepo,              &
     &                    FeOH3depo, FePflux, MnO2depo
      CHARACTER(len=80) msg

!............................ statements ..................................

!     check memory allocated to output variables
       IF (ip(1)< noutdia+nforcsdia) THEN
         CALL rexit("nout should be at least nout")
       ENDIF  

! from Conc to fdet, sdet, o2,...
       DO I = 1, N
        Fdet(I)   = Conc(       1+I)
        Sdet(I)   = Conc(   Np1+1+I)
        O2(I)     = Conc( 2*Np1+1+I)
        NO3(I)    = Conc( 3*Np1+1+I)
        NO2(I)    = Conc( 4*Np1+1+I)
        NH3(I)    = Conc( 5*Np1+1+I)
        DIC(I)    = Conc( 6*Np1+1+I)
        Fe(I)     = Conc( 7*Np1+1+I)
        FeOH3(I)  = Conc( 8*Np1+1+I)
        H2S(I)    = Conc( 9*Np1+1+I)
        SO4(I)    = Conc(10*Np1+1+I)
        CH4(I)    = Conc(11*Np1+1+I)
        PO4(I)    = Conc(12*Np1+1+I)
        FeP(I)    = Conc(13*Np1+1+I)
        CaP(I)    = Conc(14*Np1+1+I)
        Pads(I)   = Conc(15*Np1+1+I)
        ALK(I)    = Conc(16*Np1+1+I)
	 FeOH3B(I)  = Conc(17*Np1+1+I)
	 Mn(I)      = Conc(18*Np1+1+I)
	 MnO2(I)    = Conc(19*Np1+1+I)
	 MnO2B(I)   = Conc(20*Np1+1+I)
        FeS(I)     = Conc(21*Np1+1+I)
        FeS2(I)    = Conc(22*Np1+1+I)
        S0(I)      = Conc(23*Np1+1+I)
        MnCO3(I)   = Conc(24*Np1+1+I)
        FeCO3(I)   = Conc(25*Np1+1+I)
        MnS(I)   = Conc(26*Np1+1+I)
       ENDDO
       Fdetbw   = Conc(       1)
       Sdetbw   = Conc(   Np1+1)
       O2bw     = Conc( 2*Np1+1)
       NO3bw    = Conc( 3*Np1+1)
       NO2bw    = Conc( 4*Np1+1)
       NH3bw    = Conc( 5*Np1+1)
       DICbw    = Conc( 6*Np1+1)
       Febw     = Conc( 7*Np1+1)
       FeOH3bw  = Conc( 8*Np1+1)
       H2Sbw    = Conc( 9*Np1+1)
       SO4bw    = Conc(10*Np1+1)
       CH4bw    = Conc(11*Np1+1)
       PO4bw    = Conc(12*Np1+1)
       FePbw    = Conc(13*Np1+1)
       CaPbw    = Conc(14*Np1+1)
       Padsbw   = Conc(15*Np1+1)
       Alkbw    = Conc(16*Np1+1)
       FeOH3BBW = Conc(17*Np1+1)
       Mnbw     = Conc(18*Np1+1)
       MnO2bw   = Conc(19*Np1+1)
       MnO2Bbw  = Conc(20*Np1+1)
       FeSbw    = Conc(21*Np1+1)
       FeS2bw   = Conc(22*Np1+1)
       S0bw     = Conc(23*Np1+1)
       MnCO3bw  = Conc(24*Np1+1)
       FeCO3bw  = Conc(25*Np1+1)
       MnSbw    = Conc(26*Np1+1)
! --------------------------------------------------------------------------
       FDETFLux = Carbonflux*pFast
       SDETflux = Carbonflux - FDETflux

!       FeOH3surfflux = FeOH3depo*pFast
 !      FeOH3Bsurfflux = FeOH3depo - FeOH3surflux 

       FDETdepo  = FdetBW*Cfall
       SDETdepo  = SdetBW*Cfall
       FePdepo   = FePBW*FePfall
       FeOH3depo = FeOH3BW*FeOH3fall
       CaPdepo   = CaPBW*CaPfall
       FeOH3depo = FeOH3BBW*FeOH3fall    ! FeOH3depo (need to formulate correctly for FeoxB)
       MnO2depo  = MnO2bw*FeOH3fall

       CALL FESDIAtransolid(FDETdepo,SDETdepo,FePdepo,0.D0,                     &
     &                      FeOH3depo,CaPdepo,FeOH3depo, MnO2depo,              &
     &                      MnO2depo)

       if (Hwater > 0) THEN
       FePflux = 0.d0
        dFDETBW  = (FDETflux -   FDETdepo)/ Hwater
        dSDETBW  = (SDETflux -   SDETdepo)/ Hwater
        dFePBW   = (FePflux  -   FePdepo )/ Hwater
        dFeOH3BW = (FeOH3flux- FeOH3depo )/ Hwater
!        dFeOH3BW  = (FeOH3surfflux - FeOH3depo )/ Hwater
        dCaPBW   = (CaPflux  -   CaPdepo )/ Hwater
!	dFeOH3BBW = (FeOH3surfflux- FeOH3Bdepo )/ Hwater
	dFeOH3BBW = (FeOH3flux- FeOH3depo )/ Hwater
	dMnO2BW   = (MnO2flux - MnO2depo )/ Hwater
	dMnO2BBW  = (MnO2flux - MnO2depo )/ Hwater
       ! dFeSBW    = (FeSflux - MnO2depo )/ Hwater
       ! dFeS2BW   = (FeSflux - MnO2depo )/ Hwater
       ! dS0BW     = (FeSflux - MnO2depo )/ Hwater
       ! dMnCO3BW  = (FeSflux - MnO2depo )/ Hwater
       ! dFeCO3BW  = (FeSflux - MnO2depo )/ Hwater
       ELSE
        dFDETBW = 0.D0
        dSDETBW = 0.D0
        dFePBW  = 0.D0
        dFeOH3BW  = 0.D0
        dCaPBW  = 0.D0
	 dFeOH3BBW = 0.D0
	 dMnO2BW = 0.D0
	 dMnO2BBW = 0.D0
       ENDIF
       dPadsBW   = 0.D0
       dFeSBW    = 0.D0
       dFeS2BW   = 0.D0
       dS0BW     = 0.D0
       dMnCO3BW  = 0.D0
       dFeCO3BW  = 0.D0
       dMnSBW    = 0.D0
       
       CALL FESDIAtranliquid(O2bw, NO3bw, NO2bw, NH3bw, CH4bw, PO4bw,            &
     &                      Febw, H2Sbw, SO4bw, DICbw, ALKbw, Mnbw)

       if (Hwater > 0) THEN
        dO2BW   = (bwO2  -  O2BW) * relax - O2flux  /Hwater
        dNO3BW  = (bwNO3 - NO3BW) * relax - NO3flux /Hwater
        dNO2BW  = (bwNO2 - NO2BW) * relax - NO2flux /Hwater
        dNH3BW  = (bwNH3 - NH3BW) * relax - NH3flux /Hwater
        dCH4BW  = (bwCH4 - CH4BW) * relax - CH4flux /Hwater
        dPO4BW  = (bwPO4 - PO4BW) * relax - PO4flux /Hwater
        dFeBW   = (bwFe  - FeBW ) * relax -  Feflux /Hwater
        dH2SBW  = (bwH2S - H2SBW) * relax - H2Sflux /Hwater
        dSO4BW  = (bwSO4 - SO4BW) * relax - SO4flux /Hwater
        dDICBW  = (bwDIC - DICBW) * relax - DICflux /Hwater
        dALKBW  = (bwALK - ALKBW) * relax - ALKflux /Hwater
        dMnBW  =  (bwMn  - MnBW)  * relax - Mnflux  /Hwater
       ELSE
        dO2BW   = (bwO2  -  O2BW) * relax
        dNO3BW  = (bwNO3 - NO3BW) * relax
        dNO2BW  = (bwNO2 - NO2BW) * relax
        dNH3BW  = (bwNH3 - NH3BW) * relax
        dCH4BW  = (bwch4 - ch4BW) * relax
        dPO4BW  = (bwPO4 - PO4BW) * relax
        dFeBW   = (bwFe  - FeBW ) * relax
        dH2SBW  = (bwH2S - H2SBW) * relax
        dSO4BW  = (bwSO4 - SO4BW) * relax
        dDICBW  = (bwDIC - DICBW) * relax
        dALKBW  = (bwALK - ALKBW) * relax
        dMnBW   = (bwMn  - MnBW)  * relax
       ENDIF

       CALL FESDIAbiochem

       CALL MPBFESDIAsimple 

       CALL FESDIAout     (yout)

! from dfdet, dsdet, do2,... to dconc
       DO I = 1, N
         dConc(       1+I) =  dFdet(I)
         dConc(   NP1+1+I) =  dSdet(I)
         dConc( 2*NP1+1+I) =  dO2(I)
         dConc( 3*NP1+1+I) =  dNO3(I)
         dConc( 4*NP1+1+I) =  dNO2(I)
         dConc( 5*NP1+1+I) =  dNH3(I)
         dConc( 6*NP1+1+I) =  dDIC(I)
         dConc( 7*NP1+1+I) =  dFe(I)
         dConc( 8*NP1+1+I) =  dFeOH3(I)
         dConc( 9*NP1+1+I) =  dH2S(I)
         dConc(10*NP1+1+I) =  dSO4(I)
         dConc(11*NP1+1+I) =  dCH4(I)
         dConc(12*NP1+1+I) =  dPO4(I)
         dConc(13*NP1+1+I) =  dFeP(I)
         dConc(14*NP1+1+I) =  dCaP(I)
         dConc(15*NP1+1+I) =  dPads(I)
         dConc(16*NP1+1+I) =  dALK(I)
	 dConc(17*NP1+1+I) =  dFeOH3B(I)
	 dConc(18*NP1+1+I) =  dMn(I)
	 dConc(19*NP1+1+I) =  dMnO2(I)
	 dConc(20*NP1+1+I) =  dMnO2B(I)
        dConc(21*NP1+1+I) =  dFeS(I)
        dConc(22*NP1+1+I) =  dFeS2(I)
        dConc(23*NP1+1+I) =  dS0(I)
        dConc(24*NP1+1+I) =  dMnCO3(I)
        dConc(25*NP1+1+I) =  dFeCO3(I)
        dConc(26*NP1+1+I) =  dMnS(I)
       ENDDO 
       dConc(       1) =  dFdetBW
       dConc(   NP1+1) =  dSdetBW
       dConc( 2*NP1+1) =  dO2BW
       dConc( 3*NP1+1) =  dNO3BW
       dConc( 4*NP1+1) =  dNO3BW
       dConc( 5*NP1+1) =  dNH3BW
       dConc( 6*NP1+1) =  dDICBW
       dConc( 7*NP1+1) =  dFeBW
       dConc( 8*Np1+1) =  dFeOH3BW
       dConc( 9*NP1+1) =  dH2SBW
       dConc(10*NP1+1) =  dSO4BW
       dConc(11*NP1+1) =  dCH4BW
       dConc(12*NP1+1) =  dPO4BW
       dConc(13*NP1+1) =  dFePBW
       dConc(14*NP1+1) =  dCaPBW
       dConc(15*NP1+1) =  0.d0   ! Pads
       dConc(16*NP1+1) =  dALKbw
       dConc(17*NP1+1) =  dFeOH3BBW
       dConc(18*NP1+1) =  dMnBW
       dConc(19*NP1+1) =  dMnO2BW
       dConc(20*NP1+1) =  dMnO2BBW
       dConc(21*NP1+1) =  0.D0 
       dConc(22*NP1+1) =  0.D0 
       dConc(23*NP1+1) =  0.D0 
       dConc(24*NP1+1) =  0.D0 
       dConc(25*NP1+1) =  0.D0 
       dConc(26*NP1+1) =  0.D0 

      END SUBROUTINE fesdiamodBW

