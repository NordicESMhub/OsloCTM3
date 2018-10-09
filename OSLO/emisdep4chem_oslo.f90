!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Calculate emission and deposition as input to chemistry.
!//=========================================================================
module emisdep4chem_oslo
  !// ----------------------------------------------------------------------
  !// Routines for generating emissions for treatment in chemistry.
  !// Contains:
  !//   - subroutine emis4chem
  !//   - subroutine getEMISX
  !//   - subroutine getVDEP_oslo
  !//
  !// Ole Amund Sovde, September - October 2009
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'emisdep4chem_oslo.f90'
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine emis4chem(BEMIS,BTEM,UTTAU,MP)
    !// --------------------------------------------------------------------
    !// Calculates emissions for use in Oslo chemistry. Looping through the
    !// emission arrays ala p-srce.f, we sum up the emissions for each
    !// transported component. Final unit is still kg/s, it is converted
    !// to molec/(cm3*s) in chemistry sections.
    !//
    !// Operates on one IJ-block.
    !//
    !// This routine is called every NOPS (operator split), with best
    !// choice of NROPSM=3, i.e. each hour when NRMETD=8.
    !//
    !// It may be that if 3D emissions never have diurnal variations, they
    !// can be put in a seperate array and only be updated every NMET.
    !//
    !// Ole Amund Sovde, September - October 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, LPAR, NPAR, IDBLK, JDBLK, LWEPAR, &
         LEMISDEP_INCHEM
    use cmn_ctm, only: JMON, JDATE, NTM, GMTAU, XGRD, YGRD, IDAY, &
         MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, LMMAP,ETAA,ETAB
    use cmn_chem, only: E2DS, E3DSNEW, NEMLIT, LITSRC, LITFAC, NLIT, &
         NE2TBL, NY2TBL, NM2TBL, E2LTBL, E2STBL, &
         NE3TBL, NY3TBL, NM3TBL, E3LTBL, E3STBL, TNAME
    use cmn_met, only: MYEAR, SFU, SFV, LBLH, ZOFLE, P
    use cmn_parameters, only: R_UNIV
    use cmn_oslo, only: trsp_idx, FF_TYPE, NEFIR, &
         ECOMP_FIR, EMIS_FIR, EPAR_FIR_LM, &
         E2CTBL, E2LocHourTBL, E2LocHourSCALE, &
         E22dTBL, E22dSCALE, &
         E2vertTBL, E2vertSCALE, NE2vertLVS
    use emissions_aircraft, only: EPAR_AC, ECOMP_TRNR, EMIS_AC
    use emissions_megan, only: add_meganBiogenic
    use emissions_ocean, only: add_oceanOCemis
    use emissions_volcanoes, only: add_volcEMIS
    use sulphur_oslo, only: DMSseaconc
    use utilities_oslo, only: SZA_PN
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)  :: MP                     !// IJ-block index
    real(r8), intent(in) :: UTTAU                  !// Time UTC (GMTAU)
    real(r8), intent(in) :: BTEM(LPAR,IDBLK,JDBLK) !// Temperature

    !// Output: Grid box emissions for this IJ-block
    real(r8), intent(out) :: BEMIS(LPAR,NPAR,IDBLK,JDBLK)

    !// Locals
    real(r8) :: ADDTC, SCALTEMP, SCALLIGHT, SCALING
    real(r8) :: vv, rk600 !// DMS short term variations
    real(r8) :: SZA, U0, SOLF        !// To find night/day
    integer :: I,II, J,JJ, L, N, NN, M, MM, LOCHR, PN
    real(r8) :: local_hour(IDBLK), dPtot, P1, P2,ztot,Z1,Z2
    real(r8), parameter :: LOCDT=24._r8/IPAR
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis4chem'
    !// --------------------------------------------------------------------

    !// Just to make sure we really belong here
    if (.not. LEMISDEP_INCHEM) then
       write(6,'(a)') f90file//':'//subr// &
            ': Cannot handle LEMISDEP_INCHEM == false'
       stop
    end if

    !// Initialize this IJ-block
    BEMIS(:,:,:,:) = 0._r8

    !// Local hours for IDBLK indices
    do I = MPBLKIB(MP), MPBLKIE(MP)
       II   = I - MPBLKIB(MP) + 1
       !// Find local hour LOCHR in grid box i (+1 to get index).
       local_hour(II) = mod((UTTAU + LOCDT*real(i-1,r8)),24._r8) + 1._r8
    end do

    !// 2D emissions
    do M = 1, NE2TBL
      !// Skip emissions if year does not match year-table
      if (.not. (NY2TBL(M).eq.MYEAR .or. NY2TBL(M).eq.9999)) cycle

      if (NM2TBL(M).eq.JMON .or. NM2TBL(M).eq.99) then
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ   = J - MPBLKJB(MP) + 1
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II   = I - MPBLKIB(MP) + 1

          !// Calculate scaling factors that may be used below
          !// Temperature scaling (see explanation below)
          if (E22dTBL(M) .gt. 0) then
            SCALTEMP = exp( (95000._r8 * (BTEM(1,II,JJ) - 303._r8))/ &
                                 (R_UNIV * 303._r8 * BTEM(1,II,JJ)) ) &
                 / (1._r8 + exp( (230000._r8 * (BTEM(1,II,JJ) - 314._r8))/ &
                                  (R_UNIV * 303._r8 * BTEM(1,II,JJ)) ) )

            !// Find solar zenith angle (SZA) and if polar night (PN=1)
            call SZA_PN(UTTAU,IDAY,YGRD(J),XGRD(I), SZA,PN)
          else
            !// Should not need to set this, as the it is only used
            !// inside similar if-test (see below).
            SCALTEMP = 1._r8
          end if

          do N = 1, NPAR
           if (E2LTBL(N,M)) then

            !// Get emissions
            !// Scale with scaling factor for N & M (unit kg/s)
            ADDTC = E2STBL(N,M)*E2DS(I,J,1,M) !// kg/s

            if (ADDTC .eq. 0._r8) cycle

            !// Add scaling for diurnal variations
            if (E2LocHourTBL(M) .gt. 0) then
              !// Find local hour LOCHR in grid box i (+1 to get index).
              LOCHR = mod(int(UTTAU + LOCDT*(i-1)),24) + 1
              !// Category number for table is E2CTBL(M)
              !// Diurnal number for table is E2LocHourTBL(M)
              ADDTC = ADDTC * E2LocHourSCALE(LOCHR,E2CTBL(M),E2LocHourTBL(M))
            end if

            !// Add scaling for 2d variations
            !// We have checked that only one of E2LocHourTBL(M)
            !// and E22dTBL(M) is set.
            if (E22dTBL(M) .gt. 0) then
              !// Will use 2d scaling.
              !// 1. Temperature and light scaling of natural emissions.
              !//    Applies to monthly emissions. Method is based on
              !//    Guenther et al. (1995, JGR, doi:10.1029/94JD02950),
              !//                     exp((Ct1*(T-Ts))/(R*Ts*T))
              !//         E = Eadj * -----------------------------------
              !//                     1.d0 + exp((Ct2*(T-Tm))/(R*Ts*T))
              !//    Empirical constants: Tm=314K, Ct1=95000, Ct2=230000.
              !//    Eadj is emissions given at a standard temperature Ts.
              !//
              !//    - This treatment uses only daytime temperatures, thus
              !//      a day/night adjustment must be carried out also.
              !//
              !//    - Surface temperatures have previously been put out
              !//      for 1997-2010 for each meteorological time step.
              !//      For each hour we calculate the fraction
              !//               exp((Ct1*(T-Ts))/(R*Ts*T))
              !//         CT = -----------------------------------
              !//               1.d0 + exp((Ct2*(T-Tm))/(R*Ts*T))
              !//      and generate a monthly mean of this value for daylit
              !//      hours. In the polar night, we use all 24 hours of
              !//      the day in this mean (in case there are emissions).
              !//
              !//    - We count the numbers of sunlit hours in each grid
              !//      box for each month, to get a daylight scaling
              !//      (LHmean). If there is 50% daylight during one
              !//      month, the emissions should be doubled during
              !//      daytime and zero at night.
              !//
              !//    - We find the adjusted emission field by:
              !//         Eadj = E / (CTmean * LHmean)
              !//
              !//    - For each time step in the model, if there is
              !//      sunlight or polar night, we then calculate the CT
              !//      from the current temperature, and then find the
              !//      actual emissions (i.e. scaled up or down from Ts):
              !//                        exp((Ct1*(T-Ts))/(R*Ts*T))
              !//      Eactual = Eadj * -----------------------------------
              !//                        1.d0 + exp((Ct2*(T-Tm))/(R*Ts*T))
              !//
              !//      Hot days will increase emissions compared to the
              !//      mean, and we do not take into account that plants
              !//      may die due to draught. Colder temperatures have
              !//      smaller emissions than the mean.
              !//
              !// 2. Daylight scaling.
              !//    Assume the monthly emissions are to be distributed on
              !//    daylight only.
              !//    - For each month, we find fraction of daylight in each
              !//      grid box.
              !//    - Emissions are divided by this fraction; i.e. we
              !//      increase the emission rate (since we only will use
              !//      the emissions at daytime).
              !//    - We check for sunlight and use the adjusted emission
              !//      rates.
              !//
              if (E22dTBL(M).eq.1) then
                !// Scale with both sunlight and temperature
                !// Use scaling for polar night also.
                if (PN.eq.1 .or. SZA.lt.90._r8) then
                  SCALING = SCALTEMP * E22dSCALE(II,JJ,MP,1)
                else
                  SCALING = 0._r8 !// No emissions at night
                end if
              else if (E22dTBL(M).eq.2) then
                !// Scale with sunlight
                !// Use scaling for polar night also.
                if (PN.eq.1 .or. SZA.lt.90._r8) then
                  SCALING = E22dSCALE(II,JJ,MP,2)
                else
                  SCALING = 0._r8
                end if
              else if (E22dTBL(M).eq.3) then
                !// Scale with heating degree day (HDD), defined as
                !//   HDD = max(0, 288.15K - T)
                !// Weighting by the meteorolgical time step, we get
                !// HDDsum_month for each month, and HDDsum_year for the
                !// year.
                !// E22dSCALE is then HDDsum_month/HDDsum_year, and applies
                !// for specific years.
                !// It is further multiplied by 12 to create correct
                !// scalings to be applied on monthly fields of unit
                !// kg/s (done in emisutils.f90).
                SCALING = E22dSCALE(II,JJ,MP,3)
                !//
                !// It is also possible to use 3-hourly updating of the HDD
                !// (as in Stohl etal, ACPD, 2013), but CTM3 is not set up
                !// for it yet. The time step (DT) fraction of the year is
                !// given by
                !//   HDD * DT / HDDsum_year
                !// i.e. the emissions in kg/s must first be converted to
                !// kg/yr before this scaling:
                !//   E = Ein * 31536000s/yr * HDD * DT / HDDsum_year
                !// The multiplication with DT is done in chemistry, but
                !// the rest can be done here.
                !// E22dSCALE therefore needs to be calculated as
                !// 1/HDDsum_year.
              else
                SCALING = 1._r8 !// Will not happen due to test above.
              end if
              ADDTC = ADDTC * SCALING
            end if
            if (ADDTC.ne.ADDTC) then
               write(6,'(a)') f90file//':'//subr//': ADDTC NAN',SCALING,M
               stop
            end if

            !// Put into emission array, surface level
            !BEMIS(1,N,II,JJ) = BEMIS(1,N,II,JJ) + ADDTC

            !// Distribute vertically on EM_VRTLVS layers or not
            if (E2vertTBL(M) .eq. 0) then
               !// put all in surface level (L=1)
               BEMIS(1,N,II,JJ) = BEMIS(1,N,II,JJ) + ADDTC
            else
               !// Distribute in some near-surface levels
               !// Thickness (m) ca 0-16,16-41,41-77,77-128
               !// Use LMMAP in case of collapsed layers near surface.
               do L = 1, NE2vertLVS
                  BEMIS(LMMAP(L),N,II,JJ) = BEMIS(LMMAP(L),N,II,JJ) &
                       + ADDTC * E2vertSCALE(L,E2vertTBL(M))
               end do
            end if


           end if !// if (E2LTBL(N,M)) then
          end do !// do N = 1, NPAR
         end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
        end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
      end if
    end do !// do M = 1, NE2TBL


    !// 3D emissions
    do M = 1, NE3TBL
      !// Skip emissions if year does not match year-table
      if (.not. (NY3TBL(M).eq.MYEAR .or. NY3TBL(M).eq.9999)) cycle

      if (NM3TBL(M).eq.JMON .or. NM3TBL(M).eq.99) then
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ   = J - MPBLKJB(MP) + 1
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II   = I - MPBLKIB(MP) + 1
          do N = 1, NTM
            if (E3LTBL(N,M))  then
              do L=1,LPAR
                !// Scale with scaling factor for N & M
                !ADDTC = E3STBL(N,M)*E3DS(I,J,L,M)
                ADDTC = E3STBL(N,M)*E3DSNEW(L,I,J,M) !// kg/s
                !// Put into emission array, level L
                BEMIS(L,N,II,JJ) = BEMIS(L,N,II,JJ) + ADDTC
              end do
             end if
          end do !// do N = 1, NTM
         end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
        end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
      end if
    end do



    !// Forest fires
    if (FF_TYPE .gt. 0) then
      !// Forest fires are included.
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP),MPBLKIE(MP)
         II = I - MPBLKIB(MP) + 1

         !// Following Roberts et al. (2009, doi: 10.5194/bg-6-849-2009),
         !// we emit more between 7 and 19 local time.
         !// First estimate is 50% at night and 150% during day.
         !// Roberts et al. may suggest even more during day.
         if (local_hour(II) .ge. 7._r8 .and. local_hour(II) .lt. 19._r8) then
            !SCALLIGHT = 1.5_r8 ! Emit 50% more for 7-19
            SCALLIGHT = 1.9_r8 ! Emit 90% more for 7-19
            !SCALLIGHT = 1.0_r8 ! Emit 100% for 7-19
         else
            !SCALLIGHT = 0.5_r8 ! Emit 50% for 19-7
            SCALLIGHT = 0.1_r8 ! Emit 10% for 19-7
            !SCALLIGHT = 1.0_r8 ! Emit 100% for 19-7
         end if

        do MM = 1, NEFIR
          !// Get transport number of component emitted
          N = ECOMP_FIR(MM)
          !// Skip if not included
          if (N .le. 0) cycle

          !// Put into BEMIS

          if (FF_TYPE .eq. 6 .or. FF_TYPE .eq. 7) then
             !// Distribute emissions in PBL, according to air mass.
             !// Combined with putting all into L=1 of EMIS_FIR

             !// Distribute emissions according to air mass
             P2 = ETAA(1) + ETAB(1) * P(I,J) ! Will use this in loop below
             P1 = ETAA(LBLH(I,J) + 1) + ETAB(LBLH(I,J) + 1) * P(I,J) !// Top of LBLH
             dPtot = P2 - P1
             if (dPtot .le. 0._r8) then
                write(6,'(a,2es14.5,i5)') f90file//':'//subr// &
                     ': dPtot <= 0: ',P2,P1, LBLH(I,J)
                stop
             end if
             do L = 1, LBLH(I,J)
                P1 = P2 !// i.e. ETAA(L) + ETAB(L) * P(I,J), bottom of layer
                P2 = ETAA(L+1) + ETAB(L+1) * P(I,J) !// Top of layer
                SCALING = SCALLIGHT * (P1 - P2) / dPtot
                if (SCALING .lt. 0._r8) then
                   write(6,'(a,3es14.5,i5)') f90file//':'//subr// &
                        ': SCALING IS NEGATIVE',P2,P1,dPtot, LBLH(I,J)
                   stop
                end if
                !// Distribute the surface emission vertically
                ADDTC = EMIS_FIR(1,MM,II,JJ,MP) & !// kg/s
                     * SCALING
                !// Put into emission array, level L
                BEMIS(L,N,II,JJ) = BEMIS(L,N,II,JJ) + ADDTC
             end do !// do L = 1, LBLH(I,J)
          else if (FF_TYPE .eq. 8 .or. FF_TYPE .eq. 9 &
                   .or. FF_TYPE .eq. 51) then
             !// Distribute emissions in PBL, according to dZ.
             !// Combined with putting all into L=1 of EMIS_FIR

             !// Distribute emissions according to air mass
             Z2 = ZOFLE(1,I,J)    ! Will use this in loop below
             Z1 = ZOFLE(LBLH(I,J)+1,I,J) !// Top of LBLH
             ztot = Z1 - Z2
             do L = 1, LBLH(I,J)
                Z1 = Z2 !// i.e. ZOFLE(L,I,J) bottom of layer
                Z2 = ZOFLE(L+1,I,J) !// Top of layer
                SCALING = SCALLIGHT * (Z2 - Z1) / ztot
                !// Distribute the surface emission vertically
                ADDTC = EMIS_FIR(1,MM,II,JJ,MP) & !// kg/s
                     * SCALING
                !// Put into emission array, level L
                BEMIS(L,N,II,JJ) = BEMIS(L,N,II,JJ) + ADDTC
             end do !// do L = 1, LBLH(I,J)

          else
             !// Use vertical predefined distribution
             do L = 1, EPAR_FIR_LM
                !// Emissions are kg/s (scaled to molec/cm3/s in getEMISX)
                ADDTC = EMIS_FIR(L,MM,II,JJ,MP) * SCALLIGHT

                !// Put into emission array, level L
                BEMIS(L,N,II,JJ) = BEMIS(L,N,II,JJ) + ADDTC
             end do !// do L = 1,EPAR_FIR_LM
          end if

        end do !// do MM = 1, NEFIR
       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
      end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    end if !// if (FF_TYPE .gt. 0) then


    !// Lightning emissions:
    !// LITSRC is determined for each met field - called by main
    if (NEMLIT .gt. 0) then
      !// LITSRC is redefined for CTM3: It is given as fraction of total
      !// activity per seconds.
      !// LITFAC(M) scales the fraction/s to kg/s.
      do M = 1, NEMLIT
        N  = NLIT(M)
       do J = MPBLKJB(MP),MPBLKJE(MP)
         JJ = J - MPBLKJB(MP) + 1
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II = I - MPBLKIB(MP) + 1

          do L = 1, LWEPAR
            !// Emissions are kg/s (scaled to molec/cm3/s in getEMISX)
            ADDTC = LITSRC(L,I,J)*LITFAC(M) !// kg/s
            !// Put into emission array, level L
            BEMIS(L,N,II,JJ) = BEMIS(L,N,II,JJ) + ADDTC
          end do !// do L = 1, LWEPAR

        end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
       end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
      end do !// do M = 1, NEMLIT
    end if !// if (NEMLIT .gt. 0) then


    !// Aircraft emissions:
    do M = 1, EPAR_AC
      N = ECOMP_TRNR(M)
      if (N .gt. 0) then
        !// Component is included
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ = J - MPBLKJB(MP) + 1
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II = I - MPBLKIB(MP) + 1
           do L = 1, LPAR
             !// Emissions are kg/s (scaled to molec/cm3/s in getEMISX)
             ADDTC = EMIS_AC(L,M,II,JJ,MP) !// kg/s
             !// Put into emission array, level L
             BEMIS(L,N,II,JJ) = BEMIS(L,N,II,JJ) + ADDTC
           end do
         end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
        end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
      end if
    end do

    !// Short term variations
    !// DMS
    do N = 1, NPAR
      if (TNAME(N).eq.'DMS') then
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ = J - MPBLKJB(MP) + 1
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II = I - MPBLKIB(MP) + 1
           vv = SFU(I,J)**2 + SFV(I,J)**2 !// "10m" wind squared
           !// Velocity based on Nightingale et al., Global Biogeochemical
           !// Cycles vol.14, no.1, p 373-387, 2000
           !// Change from cm/hr to m/s
           rk600 = (0.222_r8*vv + 0.333_r8*vv**0.5_r8) * 1.e-2_r8 / 3600._r8
           !// Scale DMSseaconc from kg/m to kg/s and add to tracer array.
           !// Then scale to molec/cm3/s
           ADDTC = rk600 * DMSseaconc(I,J,JMON)  !// kg/s
           !// Put into emission array, surface level
           BEMIS(1,N,II,JJ) = BEMIS(1,N,II,JJ) + ADDTC
         end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
        end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
        exit
      end if !// if (TNAME(N).eq.'DMS') then
    end do !// do N = 1, NPAR


    !// Volcanoes SO2 from HTAP or other dataset.
    call add_volcEMIS(BEMIS,1._r8,JDATE,JMON,MP)

    !// Oceanic carbon emissions
    call add_oceanOCemis(BEMIS,1._r8,MP)

    !// Biogenic emissions (MEGANv2.10)
    call add_meganBiogenic(BEMIS, 1._r8, MP)


    !// --------------------------------------------------------------------
  end subroutine emis4chem
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getEMISX(OCEMIS,DVCOL,I,J,II,JJ,L_START,L_END, &
       LPAR, MP, NPAR, TRACER_ID_MAX, chem_idx, TMASS, EMISX)
    !// --------------------------------------------------------------------
    !// Get column emissions from layer L_START to L_END, retrieved from
    !// OCEMIS(LPAR,NPAR) and put into EMISX(TRACER_ID_MAX,LPAR).
    !//
    !// Units of [kg/s] are converted to [molec/cm3/s].
    !//
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: AVOGNR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LPAR, MP, NPAR, TRACER_ID_MAX
    real(r8), intent(in), dimension(LPAR,NPAR) :: OCEMIS
    real(r8), intent(in), dimension(LPAR) :: DVCOL
    integer, intent(in) :: I,J, II,JJ, L_START,L_END
    integer, intent(in), dimension(NPAR) :: chem_idx
    real(r8), intent(in), dimension(NPAR) :: TMASS
    !// Output
    real(r8), intent(out), dimension(TRACER_ID_MAX,LPAR) :: EMISX

    !// Locals
    integer :: TRACER_ID, N, L
    real(r8)  :: RTMP
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis4chem'
    !// --------------------------------------------------------------------

    !// Initialize (may want to remove it for efficiency)
    EMISX(:,:) = 0._r8

    !// Only emissions for transported species
    do N = 1,NPAR
       TRACER_ID = chem_idx(N)
       RTMP = 1.e-3_r8 * AVOGNR / TMASS(N)
       !// No need to check TRACER_ID, it is positive for all
       !// transported tracers
       if (TRACER_ID .le. 0) then
          write(6,'(a,7i5)') f90file//':'//subr// &
               ': VERY VERY WRONG (N,ID)',N,chem_idx(N),I,J,II,JJ,MP
          stop
       end if
       do L = L_START, L_END
          EMISX(TRACER_ID,L) = OCEMIS(L,N) &
               * RTMP / DVCOL(L) !// Convert [kg/s] to [molec/cm3/s]
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine getEMISX
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getVDEP_oslo(I,J,VDEP_OSLO,VDEP,DZ,TRACER_ID_MAX, &
       chem_idx,NPAR)
    !// --------------------------------------------------------------------
    !// Get the deposition frequency [1/s] for all species.
    !// VDEP is given in m/s, so we divide by the thickness of the layer.
    !//
    !// Keep in mind that VDEP_OSLO is NOT initialized for non-transported
    !// species are. This should not be a problem; non-transported species
    !// are not deposited.
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I,J, TRACER_ID_MAX, NPAR
    real(r8), intent(in)  :: DZ
    integer, dimension(NPAR), intent(in)          :: chem_idx
    real(r8), dimension(NPAR), intent(in) :: VDEP

    !// Output
    real(r8), dimension(TRACER_ID_MAX),intent(out)  :: VDEP_OSLO

    !// Local
    integer :: N
    !// --------------------------------------------------------------------

    !// Keep in mind that VDEP_OSLO is NOT initialized for non-transported
    !// species are. This should not be a problem; non-transported species
    !// are not deposited.

    !// Convert from m/s to 1/s
    !// Could also hardcode, but it is better to do it generally.
    do N = 1, NPAR
       !// Possible to include KEDDY-stuff here instead of in
       !// drydeposition_oslo.f
       VDEP_OSLO(chem_idx(N)) = VDEP(N) / DZ

    end do

    !// --------------------------------------------------------------------
  end subroutine getVDEP_oslo
  !// ----------------------------------------------------------------------



  !// --------------------------------------------------------------------
end module emisdep4chem_oslo
!//=========================================================================
