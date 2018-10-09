!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Nitrate aerosols.
!//=========================================================================
module nitrate
  !// ----------------------------------------------------------------------
  !// MODULE: nitrate.
  !// DESCRIPTION: Contains routines for treating nitrate aerosols,
  !//              currently with equlibrium model eqsam.
  !//
  !// Contains
  !//   subroutine nitrate_init
  !//   subroutine nitrate_master
  !//
  !// Should add output diagnostics here.
  !// 
  !// To f90 and Oslo CTM3
  !// Ole Amund Sovde, November 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only:LPAR, IDBLK, JDBLK, MPBLK
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// The fine and coarse mode of H2O were only used in diagnostics (3-hourly
  !// output) and should not need to be transported. Treat them privately in
  !// this module until otherwise is needed.
  !// H2Ofine (component 78) and H2Ocoarse (component 79).
  real(r8), dimension(LPAR,IDBLK,JDBLK,MPBLK) :: H2Ofine, H2Ocoarse
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'nitrate.f90'
  !// ----------------------------------------------------------------------
  save !// All variables are to be saved.
  private
  public nitrate_init, nitrate_master
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine nitrate_init(trsp_idx,TRACER_ID_MAX,LSULFUR,LSALT)
    !// --------------------------------------------------------------------
    !// Check that all species needed to run NITRATE are included.
    !//
    !// Ole Amund Sovde, November 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: TRACER_ID_MAX, trsp_idx(TRACER_ID_MAX)
    logical, intent(in) :: LSULFUR, LSALT
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'nitrate_init'
    !// --------------------------------------------------------------------

    !// Check for needed packages to run NITRATE
    if ( (.not. LSULFUR) .or. (.not. LSALT) &
         .or. trsp_idx(73).lt.0 .or. trsp_idx(61).lt.0 &
         .or. trsp_idx(62).lt.0 .or. trsp_idx(63).lt.0 &
         .or. trsp_idx(64).lt.0 .or. trsp_idx(65).lt.0 ) then

       write(6,'(a)') f90file//':'//subr//': MISSING data for NITRATE run'
       if (.not.LSULFUR) write(6,'(a)') '  * Sulphur module must be turned on'
       if (.not.LSALT)   write(6,'(a)') '  * Salt module must be turned on'
       if (trsp_idx(73).lt.0 .or. trsp_idx(61).lt.0 .or. trsp_idx(62).lt.0 &
            .or. trsp_idx(63).lt.0 .or. trsp_idx(64).lt.0 &
            .or. trsp_idx(65).lt.0) then
          write(6,'(a)') 'Check that species 73, 61, 62, 63, 64, 65 are included'
          write(6,'(a)') ' i.e. SULPHUR and SALT packages.'
          write(6,'(a)') 'Species MUST be transported!'
          stop
       end if

    end if
    !// --------------------------------------------------------------------
  end subroutine nitrate_init
  !// ----------------------------------------------------------------------
    



  !// ----------------------------------------------------------------------
  subroutine nitrate_master(BTT, &   !I/O [kg] tracer in box
                            BTEM, &  !I [K] Temperature
                            QHUM, &  !I [kg/kg] absolute humidity
                            BDV, &   !I [m3] volume of grid box
                            AIRB, &  !I [kg] air mass in grid box
                            MP )
    !// --------------------------------------------------------------------
    !// Based on CTM2 routine metzgerequil.
    !//
    !// Purpose: Given a total conenctration of inorganic species
    !// Put this concentration in equilibrium between gas and aerosol phase

    !// Calls EQSAM
    !// Reference: http:/www.mpch-mainz.mpg.de/~metzger

    !// Note that the way this is written assumes that you have given the
    !// right molecular weight of sea salt in input file. Early model runs of
    !// seasalt used molar weight of seasalt= molar weight of air so that
    !// mixing ratio (mass) and (volume) are equal. That is just to simplify
    !// since CTM outputs mixing ratio (volume) and we really want to have
    !// seasalt mixing ratio (mass). For these runs, this has changed.
    !// Make sure you know what you do...
    !//
    !// To Oslo CTM3
    !// Ole Amund Sovde, November 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, NPAR, IDBLK, JDBLK, LOSLOCSTRAT
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TMASS
    use cmn_parameters, only: R_UNIV
    use cmn_oslo, only: trsp_idx, LMTROP, LVS2ADD2TC
    use seasalt, only: Nsaltbins, salt_trsp_idx
    use eqsam_v03d, only: eqsam_v03d_sub  !Use the eqsam module
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: &
         BTEM(LPAR,IDBLK,JDBLK), &
         AIRB(LPAR,IDBLK,JDBLK), &
         BDV(LPAR,IDBLK,JDBLK), &
         QHUM(IPAR,JPAR,LPAR)
    integer, intent(in) :: MP !// IJ-block
    !// In/Out [kg]
    real(r8), intent(inout) :: BTT(LPAR,NPAR,IDBLK,JDBLK)

    !// Parameters needed to get Relative humidity.
    !// Coefficients used to calculate saturation partial pressure of water
    !// at a certain temperature (Seinfeld & Pandis, atmospheric chemistry
    !// and physics pp 781, table 15.2)
    real(r8), parameter :: a0 = 6.1078_r8
    real(r8), parameter :: a1 = 4.4365e-1_r8
    real(r8), parameter :: a2 = 1.4289e-2_r8
    real(r8), parameter :: a3 = 2.6507e-4_r8
    real(r8), parameter :: a4 = 3.0312e-6_r8
    real(r8), parameter :: a5 = 2.0341e-8_r8
    real(r8), parameter :: a6 = 6.1368e-11_r8

    !// Parameters needed for calling Metzger's code
    integer, parameter  :: nca = 11         !Size of input array
    integer, parameter  :: nco = 35         !Size of output array

    !// LOCAL
    !// Dimensaions are really included in isrpia.inc, but I rename them here
    !// So that this code is really decoupled from ISORROPIA
    real(r8)  :: &
         TEMP(LPAR), &
         mol2umolm3(LPAR), &   !Conversion mol->umol/m3
         RH(LPAR)              !Relative humidity
    real(r8)  :: yi(LPAR,nca)    !Input, mostly in umol/m3
    real(r8)  :: yo(LPAR,nco)    !Output in umol/m3
    integer :: &
         idx_cmp, &      !Counter for salt component
         I,J,II,JJ,L, &  !Loop counters
         TNR, &          !Temprary number
         iopt, &         !Indication for solid/liquid/metastable aerosol
         imode, &        !Loop ixd for aerosol modes
         in(6)           !Needed by EQSAM, don't know what it is good for.
    real(r8) :: &
         P0, &           ![Pa] saturation water pressure
         TC, &           ![degree C] temperature
         weightvapor, &  ![kg] weight of vapor in grid cell
         pvapor          ![Pa] Partial pressure of vapor in grid cell
    real(r8)  :: LOSS, LOSSFRAC !// To remove in stratospheric 
    integer :: LMAX !// Integrate nitrate module up to this level

    !// Local parameter
    real(r8), parameter :: ZEROVAL = 1.e-20_r8
    !// Flag for writing a lot of things
    logical,parameter :: Ldebug_nit = .false.
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'nitrate_master'
    !// --------------------------------------------------------------------

    !// We want metastable aerosols
    iopt = 1

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// CTM2 does nitrate only up to max LMTROP in latitude band
       !// The same applies for sulphate.
       LMAX = maxval(LMTROP(:,J)) + LVS2ADD2TC

       !// IMPORTANT!
       !// Above LMAX, in the stratosphere, we may assume there is sulphate
       !// presenent, so NH4NO3 should be converted to HNO3. The SO4 tracer
       !// should _not_ be increased, because SO4 is never removed in the
       !// formation of (NH4)2SO4.
       !// This treatment should be included below I-loop or at the bottom,
       !// for NO3fine and NO3coarse.

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Tropopause height when using stratospheric chemistry
          if (LOSLOCSTRAT) LMAX = LMTROP(I,J)

          !// Set some variables for the vertical
          !// ------------------------------------------------------------------
          do L = 1, LMAX !// Do not calculate nitrate in stratosphere

             !// Conversion from mass to umol/m3
             !// DV is m3, mass is kg
             !// CONC = MASS * 1.d3 [g] / TMASS [g/mol] * 1.d6 [umol/mol]
             !//                / DV [m3]
             !// 1.d9 = 1.d3 (g/kg) * 1.d6 (umol/mol)
             mol2umolm3(L) = 1.e9_r8 / BDV(L,II,JJ)

             !// Temperature
             TEMP(L) = BTEM(L,II,JJ)

             !// Relative humidity
             TC = BTEM(L,II,JJ) - 273.15_r8
             if (Ldebug_nit) write(6,*) 'TC', TC

             !// Saturation water pressure [hPa]
             P0 = a0 + a1*TC + a2*TC*TC + a3*TC**3 + a4*TC**4 &
                     + a5*TC**5 + a6*TC**6
             P0 = P0 * 1.e2_r8 !// convert to Pa

             if (Ldebug_nit) then
                write(6,*)'P0',P0
                write(6,*)'Q ',L,QHUM(I,J,L)
                write(6,*)'AIR',L,AIRB(L,II,JJ)
             end if
         
             !// weight of H2O (in kg)
             weightvapor = QHUM(I,J,L) * AIRB(L,II,JJ)
             if (Ldebug_nit) write(6,*) 'weightvapor', weightvapor
         
             !// Partial pressure of H2O using ideal gas (pV=nRT)
             pvapor = (weightvapor * 1.e3_r8 / 18.02_r8) * R_UNIV &
                  * BTEM(L,II,JJ) / BDV(L,II,JJ)
             if (Ldebug_nit) write(6,*) 'Pvapor', pvapor

             !// Relative humidity
             RH(L) = Pvapor / P0        !RELATIVE HUMIDITY
             if (RH(L) .gt. 1._r8) then
                RH(L) = 1._r8
             else if (RH(L) .lt. 0._r8) then
                RH(L) = 0.01_r8
             end if
             if (Ldebug_nit) write(6,*) 'Relative humidity ', L, RH(L)
         
             !// Don't do calculation if T (degrees C) < -50 
             IF (TC .lt. -50._r8) RH(L) = 0._r8

          end do !// do L = 1, LMAX


          !// Set variables independent of mode
          do L = 1, LMAX !// Do not calculate nitrate in stratosphere

             !// Get the temperature (Kelvin)
             yi(L,1) = TEMP(L)

             !// Get the relative humidity
             yi(L,2) = RH(L)

             !// Get the total pressure ! NOT USED; send in a random number
             yi(L,11) = 1.e3_r8        !hPa
            
          end do !// do L = 1, LMAX


          !// Loop over fine and coarse mode
          !// --------------------------------------------------------------
          do imode = 1, 2

             select case(imode)


             case(1) 
                !// FINE MODE
                !// --------------------------------------------------------
                do L = 1, LMAX !// Do not calculate nitrate in stratosphere

                   !// Natrium content
                   yi(L,6) = 0._r8        !No natrium

                   !// Chloride content
                   yi(L,7) = 0._r8        !No chloride
            
                   !// Sulfate content (mass, must be converted to umolec/m3)
                   !// Total sulfate,  convert from mass to mole
                   yi(L,4) = BTT(L,trsp_idx(73),II,JJ)  / TMASS(trsp_idx(73)) &
                        * mol2umolm3(L)
            
                   !// Ammonia/ammonium content (total gas NH3 (tracer 61)
                   !// + fine NH4 (tracer 62)). Convert to correct mole
                   !// fraction.
                   yi(L,3) = ( BTT(L,trsp_idx(61),II,JJ)/TMASS(trsp_idx(61)) &
                             + BTT(L,trsp_idx(62),II,JJ)/TMASS(trsp_idx(62)) ) &
                             * mol2umolm3(L)
            
                   !// Nitrate content (Total HNO3 gas (4) + fine NO3 (64))
                   !// Convert to correct mole fraction.
                   yi(L,5) = ( BTT(L,trsp_idx(4),II,JJ)/TMASS(trsp_idx(4)) &
                             + BTT(L,trsp_idx(64),II,JJ)/TMASS(trsp_idx(64)) ) &
                             * mol2umolm3(L)

                   !// Set Kalium (K+), Calsium (Ca++) and Magnesium (Mg++)
                   !// to zero
                   yi(L,8) = 0._r8
                   yi(L,9) = 0._r8
                   yi(L,10)= 0._r8
                        
                end do !// do L = 1, LMAX

             case(2)                   
                !// COARSE MODE
                !// --------------------------------------------------------
                do L = 1, LMAX !// Do not calculate nitrate in stratosphere
            
                   !// Get Natrium
                   yi(L,6) = 0._r8
                   !// Fetch Natrium from salt species (mass; convert to
                   !// umolec/m3). Total Na (Assume equal ammounts of Na & Cl)
                   do idx_cmp = 1, Nsaltbins
                      !// Sum up all salt species
                      TNR = salt_trsp_idx(idx_cmp)
                      yi(L,6) = yi(L,6)  &
                           + BTT(L,TNR,II,JJ)/TMASS(TNR) * mol2umolm3(L)
                   end do
            
                   !// Chloride equals natrium (assume molar ratio of 1:1
                   !// in seasalt for Na and Cl)
                   yi(L,7) = yi(L,6)
            
                   !// Sulfate (we assume none in coarse mode)
                   yi(L,4) = 0._r8
            
                   !// Gas NH3 and coarse NH4 (Total NH3gas (61)
                   !// + coarse NH4 (63))
                   yi(L,3) = ( BTT(L,trsp_idx(61),II,JJ)/TMASS(trsp_idx(61)) &
                             + BTT(L,trsp_idx(63),II,JJ)/TMASS(trsp_idx(63)) ) &
                             * mol2umolm3(L)
 
                   !// Gas and coarse nitrate (Total HNO3 gas (4)
                   !// + coarse NO3 (65))
                   yi(L,5) = ( BTT(L,trsp_idx(4),II,JJ)/TMASS(trsp_idx(4)) &
                             + BTT(L,trsp_idx(65),II,JJ)/TMASS(trsp_idx(65)) ) &
                             * mol2umolm3(L)

                   !// Set Kalium (K+), Calsium (Ca++) and Magnesium (Mg++)
                   !// to zero
                   yi(L,8) = 0._r8
                   yi(L,9) = 0._r8
                   yi(L,10)= 0._r8

                end do !// do L = 1, LMAX

             case default
                write(6,'(a,i7)') f90file//':'//subr// &
                     ': error in IMODE parameter',IMODE
                stop
             end select

             !// write(6,*)'eqs calling eqsam for ICOARSE=',ICOARSE
             !// write(6,'(a,(32e10.3))')'so4',yi(:,4)
             !// write(6,'(a,(32e10.3))')'HNO3',yi(:,5)
             !// write(6,'(a,(32e10.3))')'NH4',yi(:,3)
             !// write(6,'(a,(32e10.3))')'salt',yi(:,6)


             !// CALL EQSAM (NOTE: DOUBLE PRECISION FOR ALL INPUTS)
             !// -----------------------------------------------------------
             call eqsam_v03d_sub( &
                  yi, &       !I input concentrations (umol/m3)
                  yo, &       !O output concentrations (umol/m3)
                  nca, &      !I Number of input parameters
                  nco, &      !I Number of output parameters
                  iopt, &     !I Option (normally 1 for metastable)
                  LMAX, &     !I Number of steps to perform
                  LPAR, &     !I Size of input/output blocks
                  6, &        !I Unit to write to (screen)
                  in )        !I Don't know what this is good for
      
             !// write(6,'(a,4e15.3)')  'eqs hno3(a)/nh4(a)/hno3(g)/nh4(g)', &
             !//     yo(:,20),yo(:,19),yo(:,9),yo(:,10)
 
             !// write(6,'(a,4e15.3)') 'eqs checking output, NH3',&
             !//     yo(:,3),yo(:,19),yo(:,16),yo(:,10)

             !// write(6,'(a,4e15.3)') 'eqs checking output, HNO3', &
             !//     yo(:,5),yo(:,20),yo(:,17),yo(:,9)
      

             !// Put back into BTT (mass, must convert from umol/m3)
             !// -----------------------------------------------------------
             select case(imode)

             case(1) 
                !// FINE MODE
                !// --------------------------------------------------------
                do L = 1, LMAX !// Do not calculate nitrate in stratosphere

                   !// Fine mode ammonium (NH4+ and NH3(aq)). Back to mass.
                   BTT(L,trsp_idx(62),II,JJ) = max(ZEROVAL, &
                        yo(L,19) * TMASS(trsp_idx(62))/mol2umolm3(L) )
            
                   !// Fine mode nitrate (aerosol HNO3). Back to mass.
                   BTT(L,trsp_idx(64),II,JJ) = max(ZEROVAL, &
                        yo(L,20) * TMASS(trsp_idx(64))/mol2umolm3(L) )
            
                   !// Gas HNO3 left after equil. with fine mode. Back to mass.
                   BTT(L,trsp_idx(4),II,JJ) = max(ZEROVAL, &
                        yo(L,9) * TMASS(trsp_idx(4))/mol2umolm3(L) )

                   !// Gas NH3 left after equil. with fine mode. Back to mass.
                   BTT(L,trsp_idx(61),II,JJ) = max(ZEROVAL, &
                        yo(L,10) * TMASS(trsp_idx(61))/mol2umolm3(L) )

                   !// Water; (NOTE: H2O is in ugram in Metzger's code)
                   !// Not transported, only for diagnostics
                   H2Ofine(L,II,JJ,MP) = max(ZEROVAL, yo(L,12)/mol2umolm3(L))
            
                end do !// do L = 1, LMAX

             case(2)                   
                !// COARSE MODE
                !// --------------------------------------------------------
                do L = 1, LMAX !// Do not calculate nitrate in stratosphere

                   !// Coarse ammonium (NH4+ and NH3(aq))
                   BTT(L,trsp_idx(63),II,JJ) = max(ZEROVAL, &
                        yo(L,19) * TMASS(trsp_idx(63))/mol2umolm3(L) )
            
                   !// Coarse nitrate (NO3- & HNO3(aq))
                   BTT(L,trsp_idx(65),II,JJ) = max(ZEROVAL, &
                        yo(L,20) * TMASS(trsp_idx(65))/mol2umolm3(L) )
            
                   !// Gas HNO3 left after equilibrium with coarse
                   BTT(L,trsp_idx(4),II,JJ) = max(ZEROVAL, &
                        yo(L,9) * TMASS(trsp_idx(4))/mol2umolm3(L) )
            
                   !// Gas NH3 left after equilibrium with coarse mode
                   BTT(L,trsp_idx(61),II,JJ) = max(ZEROVAL, &
                        yo(L,10) * TMASS(trsp_idx(61))/mol2umolm3(L) )

                   !// Water; (NOTE: H2O is in ugram/m3 in Metzger's code)
                   !// Not transported, only for diagnostics
                   H2Ocoarse(L,II,JJ,MP) = max(ZEROVAL, yo(L,12)/mol2umolm3(L))

                end do !// do L = 1, LMAX
         
             case default
                write(6,'(a,i7)') f90file//':'//subr// &
                     ': error in coarse/fine indicator',IMODE
             end select
     
     
             !// write(6,*)'done equilibrium indicator',imode

          end do !// do imode = 1, 2

          !// Stratospheric treatment
          do L = LMAX+1, LPAR

             !// Assume 65% to be removed (roughly 1 hour life time)
             LOSSFRAC = 0.65_r8

             !// Remove NO3fine (64); assume it produces HNO3
             LOSS = BTT(L,trsp_idx(64),II,JJ) * LOSSFRAC
             BTT(L,trsp_idx(64),II,JJ) = max(ZEROVAL, &
                  BTT(L,trsp_idx(64),II,JJ) * (1._r8 - LOSSFRAC))

             !// Move to HNO3 & NOy
             if (LOSS .gt. ZEROVAL) then
                !// Convert from mass of NO3fine to mass of HNO3
                BTT(L,trsp_idx(4),II,JJ) = BTT(L,trsp_idx(4),II,JJ) &
                     + LOSS * TMASS(trsp_idx(4)) / TMASS(trsp_idx(64))
                !// Need to add to NOy also
                if (LOSLOCSTRAT) &
                     BTT(L,trsp_idx(147),II,JJ) = BTT(L,trsp_idx(147),II,JJ) &
                     + LOSS * TMASS(trsp_idx(147)) / TMASS(trsp_idx(64))
             end if

             !// Remove NO3coarse (65); assume it produces HNO3
             LOSS = BTT(L,trsp_idx(65),II,JJ) * LOSSFRAC
             BTT(L,trsp_idx(65),II,JJ) = max(ZEROVAL, &
                  BTT(L,trsp_idx(65),II,JJ) * (1._r8 - LOSSFRAC))

             !// Move to HNO3 & NOy
             if (LOSS .gt. ZEROVAL) then
                !// Convert from mass of NO3coarse to mass of HNO3
                BTT(L,trsp_idx(4),II,JJ) = BTT(L,trsp_idx(4),II,JJ) &
                     + LOSS * TMASS(trsp_idx(4)) / TMASS(trsp_idx(65))
                !// Need to add to NOy also
                if (LOSLOCSTRAT) &
                     BTT(L,trsp_idx(147),II,JJ) = BTT(L,trsp_idx(147),II,JJ) &
                     + LOSS * TMASS(trsp_idx(147)) / TMASS(trsp_idx(65))
             end if


             !// Remove NH4fine and NH4coarse
             BTT(L,trsp_idx(62),II,JJ) = &
                  max(ZEROVAL, BTT(L,trsp_idx(62),II,JJ) * (1._r8 - LOSSFRAC))
             BTT(L,trsp_idx(63),II,JJ) = &
                  max(ZEROVAL, BTT(L,trsp_idx(63),II,JJ) * (1._r8 - LOSSFRAC))

          end do

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine nitrate_master
  !// ----------------------------------------------------------------------

      
  !// ----------------------------------------------------------------------
end module nitrate
!//=========================================================================
