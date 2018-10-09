!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Source routine ala UCI (sources as separate process).
!//=========================================================================
module source_uci
  !//-----------------------------------------------------------------------
  !// MODULE: source_uci
  !// DESCRIPTION: Routine to emit species into tracer array, i.e.
  !//              as a separate process.
  !//
  !// Contains
  !//   subroutine SOURCE
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'source_uci.f90'
  !// ----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine SOURCE (BTT,BXT,BXX,BYT,BYY,BXY,BZT,DTSRCE,MP)
    !//---------------------------------------------------------------------
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: IPAR, LPAR, NPAR, IDBLK, JDBLK, LEMISDEP_INCHEM, &
         NE2DS
    use cmn_ctm, only: JMON, JDATE, NTM, GMTAU, MPBLKIB, MPBLKIE, &
         MPBLKJB, MPBLKJE, LMMAP, ETAA, ETAB
    use cmn_chem, only: NE2TBL, NY2TBL, NM2TBL, E2LTBL, E2STBL, E2DS, &
                        NE3TBL, NY3TBL, NM3TBL, E3LTBL, E3STBL, E3DSNEW, &
                        NEMLIT, NLIT, LITSRC, LITFAC, TNAME
    use cmn_met, only: MYEAR, SFU, SFV, LBLH, P, ZOFLE
    use cmn_oslo, only: trsp_idx, FF_TYPE, NEFIR, &
         ECOMP_FIR, EMIS_FIR, EPAR_FIR_LM, &
         E2CTBL, E2LocHourTBL, E2LocHourSCALE, &
         E2vertTBL, NE2vertLVS, E2vertSCALE
    use sulphur_oslo, only: DMSseaconc
    use emissions_aircraft, only: EPAR_AC, ECOMP_TRNR, EMIS_AC
    use emissions_megan, only: add_meganBiogenic
    use emissions_ocean, only: add_oceanOCemis
    use emissions_volcanoes, only: add_volcEMIS
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in) :: DTSRCE
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
         BXT, BXX, BYT, BYY, BXY, BZT

    !// Locals
    real(r8) :: ADDTC, ZNR, E2DS6(6)
    real(r8) :: DT
    real(r8) :: vv, rk600, zhr !// DMS short term variations
    integer :: I,II, J,JJ, L, N, NN, M, MM, LOCHR
    real(r8) :: local_hour(IDBLK), dPtot, P1, P2, ztot, Z1, Z2, &
         SCALING, SCALLIGHT
    real(r8), parameter :: LOCDT = 24.d0 / IPAR
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'SOURCE'
    !// --------------------------------------------------------------------


    !// Swich for Oslo chemistry treatment of emissions. If they are to 
    !// be treated as production terms in chemistry, they will be fetched
    !// in oc_main.f, so we just return from here.
    if (LEMISDEP_INCHEM) return

    !// All emissions are kg/s
    DT = DTSRCE

    !// Local hours for IDBLK indices
    do I = MPBLKIB(MP), MPBLKIE(MP)
       II   = I - MPBLKIB(MP) + 1
       !// Find local hour LOCHR in grid box i (+1 to get index).
       local_hour(II) = mod((GMTAU + LOCDT*real(i-1,r8)),24._r8) + 1._r8
    end do


    !// Changed the loop order; loop first over emission tables, then
    !// in the order J,I,N

    !// 2D emissions (monthly or yearly)
    !//---------------------------------------------------------------------
    if (NE2DS .eq. 6)  then
      !// 2D tables and tracers:  find non-zero emission factors
      do M = 1, NE2TBL
        !// Skip emissions if year does not match year-table
        if (.not. (NY2TBL(M) .eq. MYEAR .or. NY2TBL(M) .eq. 9999)) cycle

        if (NM2TBL(M) .eq. JMON .or. NM2TBL(M) .eq. 99) then

          do J = MPBLKJB(MP), MPBLKJE(MP)
            JJ   = J - MPBLKJB(MP) + 1
            do I = MPBLKIB(MP), MPBLKIE(MP)
              II   = I - MPBLKIB(MP) + 1
              do N = 1, NTM
                !// Check if emission set is used by tracer N
                if (E2LTBL(N,M)) then
                  !// Get emissions and moments
                  do NN = 1, NE2DS
                    E2DS6(NN) = E2DS(I,J,NN,M)
                  end do
                  !// Scale with scaling factor for N & M
                  ADDTC = DT * E2STBL(N,M) * E2DS6(1)
                  !// Add scaling for diurnal variations
                  if (E2LocHourTBL(M).gt. 0) then
                    !// Find local hour LOCHR in grid box i (+1 to get index).
                    LOCHR = mod(int(GMTAU + LOCDT*(i-1)),24) + 1
                    !// Category number for table is E2CTBL(M)
                    !// Diurnal number for table is E2LocHourTBL(M)
                    ADDTC = ADDTC * E2LocHourSCALE(LOCHR,E2CTBL(M),E2LocHourTBL(M))
                  end if

                  !// Distribute vertically on E2vertLVS layers or not
                  if (E2vertTBL(M) .eq. 0) then
                     !// Put all in surface level (L=1)
                     BTT(1,N,II,JJ) = BTT(1,N,II,JJ) + ADDTC
                     BXT(1,N,II,JJ) = BXT(1,N,II,JJ) + ADDTC * E2DS6(2)
                     BXX(1,N,II,JJ) = BXX(1,N,II,JJ) + ADDTC * E2DS6(3)
                     BYT(1,N,II,JJ) = BYT(1,N,II,JJ) + ADDTC * E2DS6(4)
                     BYY(1,N,II,JJ) = BYY(1,N,II,JJ) + ADDTC * E2DS6(5)
                     BXY(1,N,II,JJ) = BXY(1,N,II,JJ) + ADDTC * E2DS6(6)
                  else
                     !// Distribute in some near-surface levels
                     !// Thickness (m) ca 0-16,16-41,41-77,77-128
                     !// Use LMMAP in case of collapsed layers near surface.
                     do L = 1, NE2vertLVS
                        BTT(LMMAP(L),N,II,JJ) = BTT(LMMAP(L),N,II,JJ) &
                             + ADDTC * E2vertSCALE(L,E2vertTBL(M))
                        BXT(LMMAP(L),N,II,JJ) = BXT(LMMAP(L),N,II,JJ) &
                             + ADDTC * E2DS6(2) * E2vertSCALE(L,E2vertTBL(M))
                        BXX(LMMAP(L),N,II,JJ) = BXX(LMMAP(L),N,II,JJ) &
                             + ADDTC * E2DS6(3) * E2vertSCALE(L,E2vertTBL(M))
                        BYT(LMMAP(L),N,II,JJ) = BYT(LMMAP(L),N,II,JJ) &
                             + ADDTC * E2DS6(4) * E2vertSCALE(L,E2vertTBL(M))
                        BYY(LMMAP(L),N,II,JJ) = BYY(LMMAP(L),N,II,JJ) &
                             + ADDTC * E2DS6(5) * E2vertSCALE(L,E2vertTBL(M))
                        BXY(LMMAP(L),N,II,JJ) = BXY(LMMAP(L),N,II,JJ) &
                             + ADDTC * E2DS6(6) * E2vertSCALE(L,E2vertTBL(M))
                     end do
                  end if

                end if
              end do
            end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
          end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
        end if
      end do
    else  ! 2-D emission mass only
      do M = 1, NE2TBL
        !// Skip emissions if year does not match year-table
        if (.not. (NY2TBL(M) .eq. MYEAR .or. NY2TBL(M) .eq. 9999)) cycle

        if (NM2TBL(M) .eq. JMON .or. NM2TBL(M) .eq. 99) then
          do J = MPBLKJB(MP), MPBLKJE(MP)
            JJ   = J - MPBLKJB(MP) + 1
            do I = MPBLKIB(MP), MPBLKIE(MP)
              II   = I - MPBLKIB(MP) + 1
              do N = 1, NTM
                if (E2LTBL(N,M))  then
                  !// Scale with scaling factor for N & M
                  ADDTC = DT * E2STBL(N,M) * E2DS(I,J,1,M)
                  !// Add scaling for diurnal variations
                  if (E2LocHourTBL(M).gt. 0) then
                    !// Find local hour LOCHR in grid box i (+1 to get index).
                    LOCHR = mod(int(GMTAU + LOCDT*(i-1)),24) + 1
                    !// Category number for table is E2CTBL(M)
                    !// Diurnal number for table is E2LocHourTBL(M)
                    ADDTC = ADDTC * E2LocHourSCALE(LOCHR,E2CTBL(M),E2LocHourTBL(M))
                  end if

                  !// Distribute vertically on E2vertLVS layers or not
                  if (E2vertTBL(M) .eq. 0) then
                     !// Put all in surface level (L=1)
                     BTT(1,N,II,JJ) = BTT(1,N,II,JJ) + ADDTC
                  else
                     !// Distribute in some near-surface levels
                     !// Thickness (m) ca 0-16,16-41,41-77,77-128
                     !// Use LMMAP in case of collapsed layers near surface.
                     do L = 1, NE2vertLVS
                        BTT(LMMAP(L),N,II,JJ) = BTT(LMMAP(L),N,II,JJ) &
                             + ADDTC * E2vertSCALE(L,E2vertTBL(M))
                     end do
                  end if

                end if
              end do
            end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
          end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
        end if
      end do
    end if


    !// 3D emissions (monthly or yearly)
    !//---------------------------------------------------------------------
    !// 3D tables and tracers:  find non-zero emission factors
    do M = 1, NE3TBL
      !// Skip emissions if year does not match year-table
      if (.not. (NY3TBL(M) .eq. MYEAR .or. NY3TBL(M) .eq. 9999)) cycle

      if (NM3TBL(M) .eq. JMON .or. NM3TBL(M) .eq. 99) then
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ   = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II   = I - MPBLKIB(MP) + 1
            do N = 1, NTM
              if (E3LTBL(N,M))  then
                do L = 1, LPAR
                  !// Scale with scaling factor for N & M
                  ADDTC = DT * E3STBL(N,M) * E3DSNEW(L,I,J,M)
                  !// Put source in L=1, no moments:
                  BTT(L,N,II,JJ) = BTT(L,N,II,JJ) + ADDTC
                end do
              end if
            end do
          end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
        end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
      end if
    end do


    !// Forest fires (short term variations)
    !//---------------------------------------------------------------------
    if (FF_TYPE.gt. 0) then
      !// Forest fires are included.
      do J = MPBLKJB(MP), MPBLKJE(MP)
        JJ = J - MPBLKJB(MP) + 1
        do I = MPBLKIB(MP), MPBLKIE(MP)
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

            !// Put into BTT

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
                  !// Put into emission array, level L
                  BTT(L,N,II,JJ) = BTT(L,N,II,JJ) &
                       + EMIS_FIR(1,MM,II,JJ,MP) * SCALING * DT
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
                  !// Put into emission array, level L
                  BTT(L,N,II,JJ) = BTT(L,N,II,JJ) &
                       + EMIS_FIR(1,MM,II,JJ,MP) * SCALING * DT
               end do !// do L = 1, LBLH(I,J)

            else
               !// Use vertical predefined distribution
               do L = 1, EPAR_FIR_LM
                  BTT(L,N,II,JJ) = BTT(L,N,II,JJ) &
                       + EMIS_FIR(L,MM,II,JJ,MP) * SCALLIGHT * DT
               end do !// do L = 1,EPAR_FIR_LM
            end if

          end do
        end do
      end do
    end if


    !// Lightning emissions:
    !// LITSRC is determined for each met field - called by main
    !//---------------------------------------------------------------------
    if (NEMLIT.gt. 0) then
      !// LITSRC is fraction of annual emissions per second.
      !// LITFAC(M) scales the fraction/s to kg/s.
      do M = 1, NEMLIT
        N  = NLIT(M)
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II = I - MPBLKIB(MP) + 1
            do L = 1, LPAR
              !// CTM3: reverse order for LITSRC
              BTT(L,N,II,JJ) = BTT(L,N,II,JJ) &
                   + LITSRC(L,I,J) * LITFAC(M) & !// kg/s
                     * DT
            end do
          end do
        end do
      end do
    end if


    !// Aircraft emissions:
    !//---------------------------------------------------------------------
    do M = 1, EPAR_AC
      N = ECOMP_TRNR(M)
      if (N .gt. 0) then
        !// Component is included
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II = I - MPBLKIB(MP) + 1
            do L = 1, LPAR
              !// Emissions are kg/s
              ADDTC = EMIS_AC(L,M,II,JJ,MP) * DT
              !// Put into emission array, level L
              BTT(L,N,II,JJ) = BTT(L,N,II,JJ) + ADDTC
            end do
          end do
        end do
      end if
    end do



    !// Other short term variations (STV)
    !//---------------------------------------------------------------------
    !// DMS
    do N = 1, NPAR
      if (TNAME(N) .eq. 'DMS') then
        zhr = 1._r8 / 3600._r8
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II = I - MPBLKIB(MP) + 1
            vv = SFU(I,J)**2 + SFV(I,J)**2 !// "10m" wind squared
            !// Velocity based on Nightingale et al., Global Biogeochemical
            !// Cycles vol.14, no.1, p 373-387, 2000
            !// Change from cm/hr to m/s
            rk600 = (0.222_r8*vv + 0.333_r8*vv**0.5_r8) * 1.e-2_r8 * zhr
            !// Scale DMSseaconc from kg/m to kg/s and add to tracer array.
            ADDTC = DT * rk600 * DMSseaconc(I,J,JMON)
            BTT(1,N,II,JJ) = BTT(1,N,II,JJ) + ADDTC
          end do
        end do
        exit !// Done DMS
      end if
    end do

    !// Volcanoes SO2 from HTAP or other dataset.
    !//---------------------------------------------------------------------
    !// Added directly to BTT.
    call add_volcEMIS(BTT,DT,JDATE,JMON,MP)

    !// Oceanic carbon emissions
    call add_oceanOCemis(BTT,DT,MP)

    !// Biogenic emissions (MEGANv2.10)
    call add_meganBiogenic(BTT,DT,MP)



    !//---------------------------------------------------------------------
  end subroutine SOURCE
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
end module source_uci
!//=========================================================================
