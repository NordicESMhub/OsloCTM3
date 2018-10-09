!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Routines for handling aircraft emissions.
!//=========================================================================
module emissions_aircraft
  !// ----------------------------------------------------------------------
  !// Treat aircraft emissions.
  !// The default is to read 1x1 data and interpolate to current resolution,
  !// but by hard-coding you can use the old horizontal-resolution specific
  !// data from CTM2.
  !//
  !// Emissions are put into EMIS_AC_IN, on flight levels and CTM horizontal
  !// resolution, and each meteorological time step it is interpolated
  !// vertically into EMIS_AC.
  !//
  !// Emissions are not used if the species listed in ECOMP_AC are not
  !// transported.
  !//
  !// How to treat emissions close to surface is controlled by LADJUST_SFC.
  !//
  !// Ole Amund Sovde, November 2010
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8, r4
  use cmn_size, only: IPAR, JPAR, LPAR, IDBLK, JDBLK, MPBLK
  use cmn_ctm, only: JMON, GMTAU, ETAA, ETAB, XDEDG, YDEDG, &
       MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
  use cmn_chem, only: TNAME
  use cmn_met, only: P
  use cmn_parameters, only: AVOGNR
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Max levels of aircraft emissions
  integer, parameter :: maxlevels = 66

  !// Number of species emitted (used to be only 1: NO)
  integer, parameter :: EPAR_AC = 19
  !// List of all these component species
  character(len=10), dimension(EPAR_AC), parameter :: &
       ECOMP_NAMES = (/'CO', 'C2H4', 'C2H6','C3H6', 'C6H14',&
       'C6HXR', 'C6HXR_SOA', &
       'CH2O','CH3CHO','NO','C3H8','ACETONE','SO2',&
       'Tolmatic','Benzene','omFF1fob','omFF1fil','bcFF1fob','bcFF1fil'/)

  !// Chemical idx and transport idx
  integer, dimension(EPAR_AC) :: ECOMP_AC, ECOMP_TRNR

  !// Diurnal variation (local hour start 00 - 23)
  !// IMPORTANT!
  !// This variation should NOT be used if the aircraft data are already
  !// given with temporal variation. This variation is superimposed
  !// in the emission routines, and should be used with care!
  !// See oc_emis4chem (oc_emisdep4chem.f) and SOURCE (p-srce_oc2.f)
  real(r8), dimension(24) :: ac_diurnal = (/0.17_r8, 0.17_r8, 0.17_r8, &
       0.17_r8, 0.17_r8, 0.17_r8, 0.39_r8, 0.6_r8, 1.5_r8, 1.5_r8, &
       1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, &
       1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 0.6_r8, 0.39_r8 /)

  !// Emissions on flight levels
  !// Data can be means of each hour.
  real(r8) :: EMIS_AC_IN(maxlevels,EPAR_AC,IPAR,JPAR,24)

  !// Final emissions array (MP-block structure)
  real(r8) :: EMIS_AC(LPAR,EPAR_AC,IDBLK,JDBLK,MPBLK)

  !// Read from Ltracer_emis.inp by emis_input in emissions_oslo.f90
  integer :: AirScenYear
  character(len=12) :: AirScen
  character(len=160) :: AirEmisPath

  !// Number of levels in dataset
  integer :: actual_levels

  !// ----------------------------------------------------------------------

  !// Vertical variables
  real(r8) :: PAC_CNTR(maxlevels), &       !// Center pressure of level
            PAC_LEDG_STD(maxlevels+1), & !// Lower edge
            ZAC_LEDG(maxlevels+1)        !// Height of box bottoms

  !// Parameters; conversions and switches
  real(r8),parameter :: ft2km=0.3048e-3_r8

  !// Surface treatment: LADJUST_SFC=.true. will adjust aircraft emissions
  !// between aircraft surface and 500hPa to models surface and 500hPa.
  !// If false, emissions below CTM surface are put into surface layer.
  logical, parameter :: LADJUST_SFC=.true.
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'emissions_aircraft.f90'
  !// ----------------------------------------------------------------------

  save
  private
  public aircraft_emis_master, EPAR_AC, ECOMP_AC, ECOMP_TRNR, EMIS_AC, &
       AirScen, AirScenYear,AirEmisPath, ac_diurnal, aircraft_h2o_zero
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine aircraft_emis_master(NMET,NOPS,LNEW_MONTH,DT)
    !// --------------------------------------------------------------------
    !// Master routine for reading and interpolating aircraft emissions.
    !//
    !// Called from outside parallel region.
    !//
    !// Ole Amund Sovde, November 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NMET,NOPS
    logical, intent(in) :: LNEW_MONTH
    real(r8), intent(in)  :: DT
    !// --------------------------------------------------------------------

    !// Zero aircraft H2O tracer below 400hPa
! trop chem will remove chemical H2O at surface 
    !call aircraft_zero_400hpa(DT)

    !// Update aircraft emisions
    if (LNEW_MONTH) then
       call aircraft_set_species()
       call aircraft_emis_update()
    end if

    !// Short term variations (in vertical) every meteorological time step
    if (NOPS.eq.1) call aircraft_emis_vinterp()


    !// --------------------------------------------------------------------
  end subroutine aircraft_emis_master
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine aircraft_set_species()
    !// --------------------------------------------------------------------
    !// Finds chemi_idx and trsp_idx for emitted species.
    !//
    !// Called from aircraft_emis_master.
    !//
    !// Amund Sovde Haslerud, January 2017
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR
    use cmn_oslo, only: trsp_idx, chem_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: NAC, N
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'aircraft_set_species'
    !// --------------------------------------------------------------------

    !// From ECOMP_NAMES, get component ids in ECOMP_AC and transport number
    !// in ECOMP_TRNR.
    ECOMP_AC(:)   = -99
    ECOMP_TRNR(:) = -99

    do NAC = 1, EPAR_AC
       do N = 1, NPAR
          if (trim(ECOMP_NAMES(NAC)) .eq. trim(TNAME(N))) then
             ECOMP_AC(NAC) = chem_idx(N)
             ECOMP_TRNR(NAC) = N
             exit
          end if
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine aircraft_set_species
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine aircraft_emis_update()
    !// --------------------------------------------------------------------
    !// Controls read-in of emission data, called by aircraft_emis_master.
    !//
    !// Ole Amund Sovde, November 2010
    !// --------------------------------------------------------------------
    use utilities_oslo, only: US76_Atmosphere
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Locals
    Integer :: I,J,L,H,LAC,NAC,IOS

    !// For converting flight altitude to pressure
    real(r8) :: dft_ac, dz_ac, z_ac, pfree
    real(r8) :: sigma, delta, theta

    !// Field to read in
    real(r4) :: in_src(IPAR,JPAR,maxlevels)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'aircraft_emis_update'
    !// --------------------------------------------------------------------

    !// Initialize
    EMIS_AC_IN(:,:,:,:,:) = 0._r8
 
    !// --------------------------------------------------------------------
    !// Define vertical variables
    !// - NASA aircraft emission arrays H2O_src and NOXsrc contain data
    !//   for layers 0-1 km, 1-2 km, ..., 22-23 km above mean sea level.
    !// - Lee aircraft emission arrays H2Osrc and NOXsrc contain data
    !//   for layers 0-610 m, 610-1220 m, ..., 13420-14030 m above mean sea
    !//   level.
    !// --------------------------------------------------------------------
    if (trim(AirScen).eq.'NONE') then
       !// No aircraft emissions
       EMIS_AC(:,:,:,:,:) = 0._r8
       return
    else if (trim(AirScen).eq.'TradeOff_1e') then
       !// Lee/QinetiQ scenarios: QinetiQ 2000 revised
       dft_ac = 2000._r8
       dz_ac  = dft_ac*ft2km
       actual_levels = 23
    else if (AirScen(1:7) .eq. 'REACT4C') then
       !// React4C
       dft_ac = 2000._r8
       dz_ac  = dft_ac*ft2km
       actual_levels = 23
    else if (trim(AirScen).eq.'Quantify_MA8') then
       !// Quantify emissions May 2008
       dft_ac = 2000._r8
       dz_ac  = dft_ac*ft2km
       actual_levels = 23
    else if (trim(AirScen).eq.'CEDS') then
       !// CEDS/CMIP6 emissions
       dft_ac = 2000._r8
       dz_ac  = dft_ac*ft2km
       actual_levels = 25
    else if (trim(AirScen).eq.'CEDS2017') then
       !// CEDS/CMIP6 emissions
       dft_ac = 2000._r8
       dz_ac  = dft_ac*ft2km
       actual_levels = 25
    else if (    trim(AirScen).eq.'z1' &     ! 2025 subsonic
             .or.trim(AirScen).eq.'z2' &     ! 2025 sub- & supersonic
             .or.trim(AirScen).eq.'z3' &     ! 2050 subsonic
             .or.trim(AirScen).eq.'z4' &     ! 2050 sub- & supersonic
             .or.trim(AirScen).eq.'z5' &     ! 2050 EI(nox) increase
             .or.trim(AirScen).eq.'z6' &     ! 2050 sup. fleet size increase
             .or.trim(AirScen).eq.'z7' &     ! 2050 mach nr decrease (-4000ft)
             .or.trim(AirScen).eq.'z8' &     ! 2050 sup. range increase
             .or.trim(AirScen).eq.'z9') Then ! 2050 sup. cruise alt. reduction
       !// Airbus (SCENIC) scenarios
       dft_ac = 1000._r8
       dz_ac  = dft_ac*ft2km
       actual_levels = 66
    else
       !// Undefined scenario
       write(6,'(a)') f90file//':'//subr// &
            ': Aircraft scenario not defined: '//trim(AirScen)
       stop
    end if

    write(6,'(a)') '-------------------------------------------------' &
         //'----------------------'

    !// Vertical spacing for aircraft data starts at 0m, so we treat the
    !// levels as edge values.
    !//
    !// Edges given in altitude
    ZAC_LEDG(1) = 0._r8
    do LAC = 2, actual_levels+1
       ZAC_LEDG(LAC) = ZAC_LEDG(LAC-1) + dz_ac
    end do

    !// Will convert this spacing to pressure levels:
    !//   Starts at sea surface; elevated airports will be put into their
    !//   respective levels. Doing a pressure interpolation using actual
    !//   atmosphere may put surface emissions under or above surface.
    PAC_LEDG_STD(:) = 0._r8
    z_ac = 0._r8
    write(6,'(a)') f90file//':'//subr//': Aircraft emissions'
    write(6,'(a)') 'Setting up standard atmosphere pressure levels [edges]'
    do LAC = 1, actual_levels + 1
       call US76_atmosphere(z_ac, sigma, delta, theta)
       !// Save pressure in hPa
       PAC_LEDG_STD(LAC) = delta*1013.25_r8
       write(6,'(a,i3,f10.1,f10.2)') '  L / Z / P: ',LAC,z_ac,PAC_LEDG_STD(LAC)
       z_ac = z_ac + dz_ac
    end do

    !// Center of pressure levels
    do LAC = 1, actual_levels
       PAC_CNTR(LAC) = 0.5_r8*(PAC_LEDG_STD(LAC) + PAC_LEDG_STD(LAC+1))
    end do


    !// Aircraft emissions info
    write(6,'(a,a)')  '* Aircraft scenario: ',AirScen
    write(6,'(a,i3)') '  Number of layers: ',actual_levels


    !// Interpolate from original data
    call read_original_res()

    !// Print total emissions
    !// Convert from NO to N
    !write(6,'(a,f16.7)') '    Total Tg(N)/year:   ',&
    !     sum(EMIS_AC_IN(:,1,:,:,:))*3600._r8*365.e-9_r8 * 14._r8/30._r8
    !write(6,'(a,f16.7)') '    Total Tg(H2O)/year: ',&
    !     sum(EMIS_AC_IN(:,2,:,:,:))*3600._r8*365.e-9_r8
    !// Convert from NO to N
    write(6,'(a)') '-------------------------------------------------' &
         //'----------------------'

    !// --------------------------------------------------------------------
  end subroutine aircraft_emis_update
  !// ----------------------------------------------------------------------




  !// Short term variations
  !// ----------------------------------------------------------------------
  subroutine aircraft_emis_vinterp()
    !// --------------------------------------------------------------------
    !// Interpolate field on pressure levels. Horizontal resolution of
    !// input data must match CTM resolution.
    !//
    !// Called from outside parallel region; interpolation is done in parallel.
    !//
    !// Ole Amund Sovde, November 2010
    !// --------------------------------------------------------------------
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Temporary arrays used in vertical interpolation
    real(r8) :: emis_tmp1D(maxlevels), emis_ctm_tmp1D(LPAR)
    real(r8) :: CTMP_CNTR(LPAR), CTMP_LEDG(LPAR+1)

    !// Other variables
    integer :: I,J,L,N,LAC,NAC,hour,H, II,JJ, MP
    real(r8) :: ESUM_AC(EPAR_AC,2), RTMP

    !// File variables
    logical :: fnr_ok
    integer :: ifnr

    real(r8) :: tmpdata(ipar,jpar,lpar)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'aircraft_emis_vinterp'
    !// --------------------------------------------------------------------

    !// Initialize
    EMIS_AC(:,:,:,:,:) = 0._r8


!$omp parallel private (I,J,L,N,NAC,hour,II,JJ,MP, CTMP_CNTR,CTMP_LEDG,&
!$omp                   emis_tmp1D, emis_ctm_tmp1D) &
!$omp          shared (ETAA,ETAB,P,EMIS_AC_IN,PAC_LEDG_STD,actual_levels, &
!$omp                  EMIS_AC,MPBLKJB,MPBLKJE,MPBLKIB,MPBLKIE,GMTAU, &
!$omp                  ECOMP_TRNR) &
!$omp          default(NONE)
!$omp do
    do MP = 1, MPBLK
      !// Loop over latitude (J is global, JJ is block)
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1

        !// Loop over longitude (I is global, II is block)
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1


          !// Time of day (UTC)
          hour = nint(GMTAU) + 1
          if (hour .eq. 25) hour = 1

          !// CTM center pressure
          do L = 1, LPAR
             CTMP_CNTR(L) = 0.5_r8*( ETAA(L) + ETAB(L+1) + &
                P(I,J)*(ETAB(L) + ETAB(L+1)) )
          end do
          !// CTM lower edge pressure
          do L = 1, LPAR+1
             CTMP_LEDG(L) = ETAA(L) + P(I,J)*ETAB(L)
          end do

          !// Initialize the whole column; only data up to actual_levels
          !// will be used
          do NAC = 1, EPAR_AC

            N = ECOMP_TRNR(NAC)
            !// For CTM3, transported species have 0<N<=NPAR.
            if (N .lt. 0) cycle

            !// Initialize emissions to be interpolated
            emis_tmp1D(:) = 0._r8

            !// Fetch column emissions for this hour
            emis_tmp1D(1:actual_levels) = EMIS_AC_IN(1:actual_levels,NAC,I,J,hour)

            !// Interpolate component if there are emissions in column
            if (sum(emis_tmp1D) .gt. 0._r8) then

              !// Initialize CTM column emissions
              emis_ctm_tmp1D(:) = 0._r8

              !// Do the interpolation
              call ac_interp(emis_tmp1D, PAC_LEDG_STD, actual_levels, &
                  emis_ctm_tmp1D, CTMP_LEDG, LPAR, LADJUST_SFC,i,j)

              !// Put into final array
              do L = 1, LPAR
                 EMIS_AC(L,NAC,II,JJ,MP) = emis_ctm_tmp1D(L)
              end do

            end if

          end do !// do NAC = 1, EPAR_AC
        end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
      end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
    end do !// do MP = 1, MPBLK
!$omp end do
!$omp end parallel


    !// Check interpolation
    ESUM_AC(:,:) = 0._r8
    do H = 1, 24
      do J = 1, JPAR
        do I = 1, IPAR
          do NAC = 1, EPAR_AC
            do LAC = 1, actual_levels
              !// Data are kg/s, sum up
              ESUM_AC(NAC,1) = ESUM_AC(NAC,1) + EMIS_AC_IN(LAC,NAC,I,J,H)
            end do
          end do
        end do
      end do
    end do
    !// Convert from kg/s to Tg/year (already summed up over 24 hours)
    ESUM_AC(:,1) = ESUM_AC(:,1)*3600._r8*365._r8*1.e-9_r8

    !// Check values after interpolation
    do MP = 1, MPBLK
      !// Loop over latitude (J is global, JJ is block)
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        !// Loop over longitude (I is global, II is block)
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// Loop over emitted components
          do NAC = 1, EPAR_AC
            !// Loop vertically
            do L = 1, LPAR
              !// Data are kg/s, sum up
              ESUM_AC(NAC,2) = ESUM_AC(NAC,2) + EMIS_AC(L,NAC,II,JJ,MP)
            end do
          end do
        end do
      end do
    end do
    !// Convert from kg/s to Tg/year (one hour data only)
    ESUM_AC(:,2) = ESUM_AC(:,2)*3600._r8*24._r8*365._r8*1.e-9_r8

    !// Compare values before and after interpolation
    do NAC = 1, EPAR_AC
       N = ECOMP_TRNR(NAC)

       if (N .gt. 0) then
         RTMP = abs((ESUM_AC(NAC,2) - ESUM_AC(NAC,1))/ESUM_AC(NAC,1))
         if (RTMP .gt. 1.e-3_r8) then
           write(6,'(a)') f90file//':'//subr// &
                ': Problems in vertical interpolation of aircraft emissions!'
           write(6,'(a)') '    Component: '//trim(TNAME(N))
           write(6,'(a,f15.7)') ' [Tg/yr] before interpolation:',ESUM_AC(NAC,1)
           write(6,'(a,f15.7)') ' [Tg/yr] after  interpolation:',ESUM_AC(NAC,2)
           stop
         end if
       end if
    end do




    !// Code to print out the original and the interpolated data for NO
    if (.false.) then
       !// Find non-used file number for input file
       fnr_ok = .true.
       ifnr = 8
       do while (fnr_ok)
          ifnr = ifnr + 1
          inquire(ifnr,opened=fnr_ok)
       end do

       open(ifnr,file='acdata2.dta',form='unformatted')
       tmpdata(:,:,:) = 0._r8
       do LAC = 1, actual_levels
          do J = 1, JPAR
             do I = 1, IPAR
                tmpdata(i,j,LAC) = EMIS_AC_IN(LAC,1,I,J,1) !// Use hour 1
             end do
          end do
       end do
       write(ifnr) real(tmpdata,r4)
       tmpdata(:,:,:)=0._r8
       do MP = 1, MPBLK
         !// Loop over latitude (J is global, JJ is block)
         do J = MPBLKJB(MP),MPBLKJE(MP)
           JJ    = J - MPBLKJB(MP) + 1

           !// Loop over longitude (I is global, II is block)
           do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, LPAR
                tmpdata(i,j,l) = EMIS_AC(L,1,II,JJ,MP)
             end do
           end do
         end do
       end do
       write(ifnr) real(tmpdata,r4)
       close(ifnr)
       write(6,'(a)') f90file//':'//subr// &
            ': STOP after writing binary data!'
       stop
    end if



    !// --------------------------------------------------------------------
  end subroutine aircraft_emis_vinterp
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ac_interp(indata,inp_ledg_org,inlm, outdata,outp_ledg,outlm, &
       LADJ_SFC,LON_IND,LAT_IND)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Interpolate e.g. aircraft emissions columnwise.
    !// In/out must be mass equivalents.
    !//
    !// Assumes that aircraft emissions are given on pressure levels.
    !// These are based on geographical height converted by using the
    !// standard atmosphere, and the surface may therefore be lower
    !// (higher pressure) than the model surface.
    !//
    !// LADJ_SFC=T will modify the pressure levels close to the surface to
    !// match the model surface pressure, so that emissions in level 1 of
    !// aircraft data always will be at surface.
    !//
    !// LADJ_SFC=F will put data below surface into the model surface layer,
    !// and if the model pressure is higher than in level 1 of aircraft data,
    !// emissions supposed to be at the surface will end up above surface in
    !// the model.
    !//
    !// Ole Amund Sovde, February-November 2010
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: inlm, outlm  !// Resolutions
    real(r8), intent(in) :: indata(inlm)  !// Emissions to be interpolated
    real(r8), intent(in) :: inp_ledg_org(inlm+1),& !// Pressure for indata
                          outp_ledg(outlm+1)     !// Pressure for outdata (CTM)
    logical, intent(in) :: LADJ_SFC          !// Switch for sfc treatment
    integer, intent(in) :: LON_IND, LAT_IND !// Lon and lat indices

    !// Input/output
    real(r8), intent(inout) :: outdata(outlm)

    !// Locals
    real(r8)  :: dp_in, dp_out, wt, pmin
    real(r8)  :: inp_ledg(inlm+1) !// inp_ledg_org can be modified by LADJ_SFC
    integer :: L, LAC,LAC500

    !// Switch to print out debug messages
    logical,parameter :: LDBG_AC=.false.
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'ac_interp'
    !// --------------------------------------------------------------------

    !// Initialize output
    outdata(:) = 0._r8

    
    if (LADJ_SFC) then
       !// Will re-fit aircraft standard pressures between sfc and 500hPa to
       !// actual (model) sfc pressure and 500hPa.

       !// Check surface pressure
       pmin = 500._r8
       do LAC = 1, 10
          if ((outp_ledg(1) - pmin) .gt. 100._r8) then
             !// Allow surface pressures down to 600hPa, otherwise adjust
             !// minimum pressure by steps of 100hPa. Not a big problems, since
             !// over mountains there is not much emissions close to surface.
             exit
          else
             if (LDBG_AC) write(6,'(a,i3,3(1x,f5.1),2(1x,i3))') &
                  'Low surface pressure:',&
                  LAC,outp_ledg(1),pmin,sum(indata(1:10)),lon_ind,lat_ind
             !// Adjust minimum pressure
             pmin = pmin - 100._r8
          end if
       end do

       !// Find level closest to pmin
       do LAC = 1, inlm
          if (inp_ledg_org(LAC) .le. pmin) then
             LAC500 = LAC
             exit
          end if
       end do

       !// Do the re-fitting from sfc to LAC500
       inp_ledg(1) = outp_ledg(1)
       do LAC = 2, LAC500
          inp_ledg(LAC) = inp_ledg(LAC-1) &
               - (inp_ledg_org(LAC-1) - inp_ledg_org(LAC)) &
               /(inp_ledg_org(1) - inp_ledg_org(LAC500)) &
               *(inp_ledg(1) - inp_ledg_org(LAC500))
       end do
       !// Keep original above LAC500
       do LAC = LAC500+1,inlm+1
          inp_ledg(LAC) = inp_ledg_org(LAC)
       end do

    else    
       !// Keep aircraft pressure and put emissions below surface into the
       !// model surface layer.
       inp_ledg(:) = inp_ledg_org(:)
    end if



    !// Loop over emission layers:
    do LAC = 1, inlm

       !// No interpolation if field is zero
       if (indata(LAC).eq.0._r8) cycle

       !// dp for indata
       dp_in = inp_ledg(LAC) - inp_ledg(LAC+1)

       !// Check if indata exist below lowest level of outdata?
       if (.not. LADJ_SFC) then

          !// Put levels below CTM surface into surface
          if ( inp_ledg(LAC  ).gt.outp_ledg(1) .and. &
               inp_ledg(LAC+1).gt.outp_ledg(1) ) then
             !// Below model surface; put into surface
             wt = 1._r8
             outdata(1) = outdata(1) + indata(LAC)*wt
             if (LDBG_AC) write(6,'(2(i2,1x),a,2f8.2,f12.6)') 1,LAC, &
                  'AC below surface a',inp_ledg(LAC),outp_ledg(1),outdata(1)
          else if ( inp_ledg(LAC  ).gt.outp_ledg(1) .and. &
                   inp_ledg(LAC+1).lt.outp_ledg(1) ) then
             wt = (inp_ledg(LAC) - outp_ledg(1))/dp_in
             outdata(1) = outdata(1) + indata(LAC)*wt
             if (LDBG_AC) write(6,'(2(i2,1x),a,2f8.2,f12.6)') 1,LAC, &
                  'AC below surface b',inp_ledg(LAC),outp_ledg(1),outdata(1)
          end if

       end if !// if (.not.LADJ_SFC) then



       !// Loop through CTM levels and put emissions into correct levels
       do L = 1, outlm
          dp_out = outp_ledg(L) - outp_ledg(L+1)


          if ( outp_ledg(L  ).le.inp_ledg(LAC  ) .and. &
               outp_ledg(L  ).gt.inp_ledg(LAC+1) ) then

             if ( outp_ledg(L+1).gt.inp_ledg(LAC+1) ) then
                !// A: CTM completely inside EMIS layer
                wt = dp_out/dp_in
                outdata(L) = outdata(L) + indata(LAC)*wt
                if (LDBG_AC) write(6,'(2(i2,1x),a,4f8.2,2f12.6)') L,LAC, &
                   'CTM inside (a) ',inp_ledg(LAC),outp_ledg(L), &
                   outp_ledg(L+1),inp_ledg(LAC+1),outdata(L),sum(outdata(1:L))
             else
                !// B: Lower part of CTM inside EMIS
                wt = (outp_ledg(L) - inp_ledg(LAC+1))/dp_in
                outdata(L) = outdata(L) + indata(LAC)*wt
                if (LDBG_AC) write(6,'(2(i2,1x),a,4f8.2,2f12.6)') L,LAC, &
                    'CTM ovlaps (b) ',inp_ledg(LAC),outp_ledg(L), &
                    inp_ledg(LAC+1),outp_ledg(L+1),outdata(L),sum(outdata(1:L))
                
             end if


          else if ( OUTP_LEDG(L  ).ge.INP_LEDG(LAC  ) .and. &
                   OUTP_LEDG(L+1).lt.INP_LEDG(LAC  ) ) then

             if ( OUTP_LEDG(L+1).lt.INP_LEDG(LAC+1) ) then
                !// C: EMIS layer completely within CTM layer
                wt = 1._r8
                outdata(L) = outdata(L) + indata(LAC)*wt
                if (LDBG_AC) write(6,'(2(i2,1x),a,4f8.2,2f12.6)')L,LAC, &
                   'AC within  (c) ',outp_ledg(L),inp_ledg(LAC), &
                   inp_ledg(LAC+1),outp_ledg(L+1),outdata(L),sum(outdata(1:L))

             else
                !// D: CTM layer starts below EMIS layer and ends
                !//    before EMIS ends
                wt = (INP_LEDG(LAC) - OUTP_LEDG(L+1))/dp_in
                outdata(L) = outdata(L) + indata(LAC)*wt
                if (LDBG_AC) write(6,'(2(i2,1x),a,4f8.2,2f12.6)')L,LAC, &
                   'AC upper c (d) ',outp_ledg(L),inp_ledg(LAC), &
                   outp_ledg(L+1),inp_ledg(LAC+1),outdata(L),sum(outdata(1:L))
             end if

          end if



       end do !// do L = 1, outlm

    end do !// do LAC = 1, inlm

    !// --------------------------------------------------------------------
  end subroutine ac_interp
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine read_original_res()
    !// --------------------------------------------------------------------
    !// Master routine for reading original res resolution data and
    !// interpolate horizontally.
    !//
    !// Ole Amund Sovde, January 2011
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    if (AirScen(1:8) .eq. 'TradeOff') then
       call read_tradeoff_original()
    else if (AirScen(1:7) .eq. 'REACT4C') then
       call read_react4c_original()
    else if (AirScen(1:8) .eq. 'Quantify') then
       call read_quantify_original()
    else if (AirScen(1:8) .eq. 'CEDS2017') then
       call read_ceds_2017()
    else if (AirScen(1:4) .eq. 'CEDS') then
       call read_ceds_original()
    end if

    !// --------------------------------------------------------------------
  end subroutine read_original_res
  !// ----------------------------------------------------------------------





  !// Read original resolution data and interpolate horizontally
  !// ----------------------------------------------------------------------
  subroutine read_tradeoff_original()
    !// --------------------------------------------------------------------
    !// Read original resolution data from TradeOff and interpolate
    !// horizontally.
    !//
    !// Ole Amund Sovde, January 2011
    !// --------------------------------------------------------------------
    use regridding, only: E_GRID
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Locals
    integer :: I,J,H,LAC,NAC,IOS, N, ICOMP
    integer :: IM_AC, JM_AC, LM_AC, MAX_ITER, grid_type

    !// For converting flight altitude to pressure
    real(r8) :: dft_ac, dz_ac, z_ac, pfree
    real(r8) :: sigma, delta, theta

    !// For getting file name
    character(len=160) :: FILENAME
    character(len=2) :: CMON
    !// File variables
    logical :: fnr_ok
    integer :: ifnr

    !// Interpolated field
    real(r8) :: in_src(IPAR,JPAR,maxlevels), R8CTM(IPAR,JPAR)
    real(r8) :: ESUM_AC(EPAR_AC)
    !// Lee grid: 1x1 deg horizontal, 610 m vertical.
    !// Horizontally: First grid box between 180W and 179W and
    !//               between 90S and 89N.
    !// Vertically: L=1 means interval between 0 and 610m 'above ground'
    real(r8)  :: I_in, J_in, NOx, Fuel, Dist, RFACT, EI_H2O
    integer :: ILee, JLee, LLee, L_in
    !// Travelled distance in Lee grid: km/(month*Leebox)
    !// Fuel consumption in Lee grid:   kg(Fuel)/(month*Leebox)
    !// NOx emissions in Lee grid:      kg(NO2)/(month*Leebox)
    real(r8), dimension(:),allocatable :: XBEDGE, YBEDGE
    real(r8), dimension(:,:,:),allocatable :: FuelLee, NOxLee
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_tradeoff_original'
    !// --------------------------------------------------------------------

    if (AirScen(1:8) .eq. 'TradeOff') then


       !// Tradeoff/QinetiQ 2000 revised

       !// Seasonal data for months 01-03, 04-06, 07-09 and 10-12.
       if (JMON.lt.4) then
          CMON='01'
       else if (JMON.ge.4 .and. JMON.le.6) then
          CMON='04'
       else if (JMON.ge.7 .and. JMON.le.9) then
          CMON='07'
       else if (JMON.ge.10 .and. JMON.le.12) then
          CMON='10'
       end if

       if (trim(AirScen) .eq. 'TradeOff_1e') then
          !// File name
          FILENAME = trim(AirEmisPath)//'TRADEOFF/'// &
               'Scenario_1_2000e_610_m_'//CMON//'_civmil_3D.dat'
       else
          write(6,'(a)') f90file//':'//subr// &
               ': No such scenario defined: '//trim(AirScen)
          stop
       end if

       !// Find non-used file number for input file
       fnr_ok = .true.
       ifnr = 8
       do while (fnr_ok)
          ifnr = ifnr + 1
          inquire(ifnr,opened=fnr_ok)
       end do

       !// Open file
       Open(ifnr,File=FILENAME,Form='formatted',Status='Old',iostat=IOS)
       if (IOS .ne. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Could not find aircraft emission file: '//trim(FILENAME)
          stop
       end if
       write(6,'(a,i2)') '  * Reading aircraft emissions for month: ',JMON
       write(6,'(4x,a)') trim(FILENAME)

       !// Original resolution
       IM_AC = 360
       JM_AC = 180
       LM_AC = actual_levels
       MAX_ITER = IM_AC*JM_AC*LM_AC
       grid_type = 2

       !// Allocate arrays
       allocate (XBEDGE(IM_AC+1),YBEDGE(JM_AC+1), &
            FuelLee(IM_AC,JM_AC,LM_AC), NOxLee(IM_AC,JM_AC,LM_AC))
       call get_xyedges(IM_AC,JM_AC,XBEDGE,YBEDGE,grid_type)

       !// Initialize
       FuelLee(:,:,:) = 0._r8
       NOxLee(:,:,:) = 0._r8

       Do N = 1, MAX_ITER
          if (trim(AirScen) .eq. 'TradeOff_1e') then
             !// For 1e scenario, dist is not included
             read(ifnr,'(f6.1,1x,f5.1,1x,i2.2,3(1x,e25.18))',iostat=IOS) &
                  I_in,J_in,L_in,Fuel,NOx
             Dist = 0._r8
          else
             read(ifnr,'(f6.1,1x,f5.1,1x,i2.2,3(1x,e25.18))',iostat=IOS) &
                  I_in,J_in,L_in,Dist,Fuel,NOx ! all other scenarios
          end if
          !// Check reading
          if (ios .gt. 0) then
             write(6,'(a,i3)') '*** Failing reading file, with iostat: ',ios
             stop
          else if (ios .eq. -1) then
             !// End of file
             exit
          else if (ios .eq. -2) then
             !// End of record
             write(6,'(a)') '*** End of record, check file/read-in!'
             stop
          end if
          !// Get indices 1-360, starting at -180W
          ILee = Nint(I_in - 0.5) + 181
          JLee = Nint(J_in - 0.5) + 91
          LLee = L_in + 1 !// Lee starts at 0
          if (ILee.lt.1.or.ILee.gt.360) print*,'ILee=',ILee,'!!!'
          if (JLee.lt.1.or.JLee.gt.180) print*,'JLee=',JLee,'!!!'
          !DistLee(ILee,JLee,LLee) = Dist
          FuelLee(ILee,JLee,LLee) = Fuel
          NOxLee(ILee,JLee,LLee)  = NOx
       End Do

       !// Close file
       close(ifnr)


       !// Interpolate horizontally over aircraft levels for each of the
       !// components.
       write(6,'(a)') '  * Interpolating horizontally to CTM grid'

       !// Sum of emissions before interpolation
       ESUM_AC(:) = 0._r8

       !// Converting factor from kg/month to kg/s
       RFACT = 12._r8/(3600._r8*24._r8*365._r8)

       !// NO emissions
       !// ------------------------------------------------------
       ICOMP = 43
       !// Find NAC
       NAC = -1
       do N = 1, EPAR_AC
          if (ECOMP_AC(N) .eq. ICOMP) then
             NAC = N
             exit
          end if
       end do

       !// Is component in ECOMP_AC?
       if (NAC .gt. 0) then

         !// Sum of original data as kg[NO]/year
         ESUM_AC(NAC) = sum(NOxLee)*12._r8*30._r8/46._r8

!$omp parallel private (LAC,R8CTM) &
!$omp          shared (NOxLee,IM_AC,JM_AC,XBEDGE,YBEDGE,in_src,LM_AC) &
!$omp          shared (XDEDG,YDEDG,RFACT) &
!$omp          default(NONE)
!$omp do
         do LAC = 1, LM_AC

           !// No need to interpolate zero field
           if (maxval(NOxLee(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
           end if
           !// Interpolate into R8CTM (no moments, only mean field
           call E_GRID(NOxLee(:,:,LAC),XBEDGE,YBEDGE,IM_AC,JM_AC, &
                  R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)

           !// Put into 3D array
           !// Convert from kg[NO2]/month to kg[NO]/s
           in_src(:,:,LAC) = R8CTM(:,:)/46._r8*30._r8*RFACT

         end do
!$omp end do
!$omp end parallel

         do J = 1, JPAR
           do I = 1, IPAR
             do LAC = 1, actual_levels
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
           end do
         end do
         !// Assume no hour-to-hour variation
         do H=2,24
           do J=1,JPAR
             do I=1,IPAR
                do LAC=1,actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
           end do
         end do

       end if !// if (NAC .gt. 0) then

       !// H2O emissions
       !// ------------------------------------------------------
       ICOMP = 114
       !// Find NAC
       NAC = -1
       do N = 1, EPAR_AC
          if (ECOMP_AC(N) .eq. ICOMP) then
             NAC = N
             exit
          end if
       end do

       !// Is component in ECOMP_AC?
       if (NAC .gt. 0) then

         !// Emission index of H2O from fuel [kg/kg]
         EI_H2O = 1.230_r8
         !// Sum of original data as kg[H2O]/year
         ESUM_AC(NAC) = sum(FuelLee)*12._r8*EI_H2O

!$omp parallel private (LAC,R8CTM) &
!$omp          shared (FuelLee,IM_AC,JM_AC,XBEDGE,YBEDGE,in_src,LM_AC) &
!$omp          shared (XDEDG,YDEDG,RFACT,EI_H2O) &
!$omp          default(NONE)
!$omp do
         do LAC = 1, LM_AC

           !// No need to interpolate zero field
           if (maxval(FuelLee(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
           end if
           !// Interpolate into R8CTM (no moments, only mean field
           call E_GRID(FuelLee(:,:,LAC),XBEDGE,YBEDGE,IM_AC,JM_AC, &
                  R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)

           !// Put into 3D array
           !// Convert to kg[H2O]/s
           in_src(:,:,LAC) = R8CTM(:,:)*RFACT*EI_H2O

         end do
!$omp end do
!$omp end parallel

         !// Put into emission array, first hour
         do J = 1, JPAR
           do I = 1, IPAR
             do LAC = 1, actual_levels
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
           end do
         end do
         !// Assume no hour-to-hour variation
         do H=2,24
           do J=1,JPAR
             do I=1,IPAR
                do LAC=1,actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
           end do
         end do

         !// H2O Emissions for 114 also applies to 148
         ICOMP = 148
         !// Find NAC for 148, call it ILEE
         ILEE = -1
         do N = 1, EPAR_AC
            if (ECOMP_AC(N) .eq. ICOMP) then
               ILEE = N
               exit
            end if
         end do
         if (ILEE .gt. 0) then
            ESUM_AC(ILEE) = ESUM_AC(NAC)
            EMIS_AC_IN(:,ILEE,:,:,:) = EMIS_AC_IN(:,NAC,:,:,:)
         end if


       end if !// if (NAC .gt. 0) then


       !// De-allocate
       deallocate (XBEDGE,YBEDGE, FuelLee, NOxLee)


    end if

    !// Print total emissions
    write(6,'(a)') '    Total emissions [Tg/year] for emitted species:'
    do NAC = 1, EPAR_AC
      N = ECOMP_TRNR(NAC)
      if (N .gt. 0) then
         write(6,'(a,2f16.7)') '    '//TNAME(N)// &
              ' before/after hor. interp.:',&
              ESUM_AC(NAC)*1.e-9_r8, &
              sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
      else
         write(6,'(a,i3,a,2f16.7)') '    CompID ',ECOMP_AC(NAC), &
              ' before/after hor. interp.:',&
              ESUM_AC(NAC)*1.e-9_r8, &
              sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
         write(6,'(15x,a)')'*** Not transported - will not be used! ***'
      end if
    end do

    !// ------------------------------------------------------------------
  end subroutine read_tradeoff_original
  !// ------------------------------------------------------------------






  !// ------------------------------------------------------------------
  subroutine read_react4c_original()
    !// ------------------------------------------------------------------
    !// Read REACT4C WP5 aircraft data and interpolate horizontally.
    !// Annual emissions.
    !//
    !// Ole Amund Sovde, February 2011
    !// ------------------------------------------------------------------
    use regridding, only: E_GRID
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------

    !// Locals
    integer :: I,J,H,LAC,NAC,IOS, N, ICOMP
    integer :: IM_AC, JM_AC, LM_AC, MAX_ITER, grid_type

    !// For converting flight altitude to pressure
    real(r8) :: dft_ac, dz_ac, z_ac, pfree
    real(r8) :: sigma, delta, theta

    !// For getting file name
    character(len=160) :: FILENAME
    !character(len=2) :: CMON
    !// File variables
    logical :: fnr_ok
    integer :: ifnr

    !// Interpolated field
    real(r8) :: in_src(IPAR,JPAR,maxlevels), R8CTM(IPAR,JPAR)
    real(r8) :: ESUM_AC(EPAR_AC)
    !// Lee grid: 1x1 deg horizontal, 610 m vertical.
    !// Horizontally: First grid box between 180W and 179W and
    !//               between 90S and 89N.
    !// Vertically: L=1 means interval between 0 and 610m 'above ground'
    real(r8)  :: lat, lon, NOx, Fuel, Dist,CO2,soot,particles
    real(r8)  :: RFACT_NOX, RFACT_FUEL, EI_H2O
    integer :: ILee, JLee, LLee, ilev,idist
    real(r8), dimension(:),allocatable :: XBEDGE, YBEDGE
    real(r8), dimension(:,:,:),allocatable :: R4C_Fuel, R4C_NO
    !// ------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_react4c_original'
    !// --------------------------------------------------------------------


    !// Annual average data.

    !// File name
    if (trim(AirScen) .eq. 'REACT4C_5a') then
       !// 2006 base case
       FILENAME = trim(AirEmisPath)//'REACT4C/BC_270810_v1_0.txt'
    else if (trim(AirScen) .eq. 'REACT4C_5b') then
       !// 2006 plus 2000ft
       FILENAME = trim(AirEmisPath)//'REACT4C/Plus_270810_v1_0.txt'
    else if (trim(AirScen) .eq. 'REACT4C_5c') then
       !// 2006 minus 2000ft
       FILENAME = trim(AirEmisPath)//'REACT4C/Minus_270810_v1_0.txt'
    else
       write(6,'(a)') f90file//':'//subr// &
            ': No such scenario defined: '//trim(AirScen)
       stop
    end if

    !// Find non-used file number for input file
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do

    !// Open file
    Open(ifnr,File=FILENAME,Form='formatted',Status='Old',iostat=IOS)
    if (IOS .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Could not find aircraft emission file: '//trim(FILENAME)
       stop
    end if
    write(6,'(a,i2)') '  * Reading aircraft emissions for month: ',JMON
    write(6,'(a)')    '    '//trim(FILENAME)
    !// Read header
    read(ifnr,*)

    !// Original resolution
    IM_AC = 360
    JM_AC = 180
    LM_AC = actual_levels
    MAX_ITER = IM_AC*JM_AC*LM_AC
    grid_type = 2

    !// Allocate arrays
    allocate (XBEDGE(IM_AC+1),YBEDGE(JM_AC+1), &
         R4C_FUEL(IM_AC,JM_AC,LM_AC), R4C_NO(IM_AC,JM_AC,LM_AC))
    call get_xyedges(IM_AC,JM_AC,XBEDGE,YBEDGE,grid_type)

    !// Initialize
    R4C_FUEL(:,:,:) = 0._r8
    R4C_NO(:,:,:) = 0._r8


    Do N = 1, MAX_ITER
       !// "LONG","LAT","FL","DIST","FUEL","NOx","CO2","Soot","Particles"
       !// -179.50,-86.50,17,312,1683.90,21.50,5304.28,0.02,2753623.90
       read(ifnr,*,iostat=IOS) &
            lon,lat,ilev,idist,Fuel,NOx,CO2,soot,particles

       !// Check reading
       if (ios .gt. 0) then
          write(6,'(a,i3)') f90file//':'//subr// &
               ': Failing reading file, with iostat: ',ios
          stop
       else if (ios .eq. -1) then
          !// End of file
          exit
       else if (ios .eq. -2) then
          !// End of record
          write(6,'(a,i3)') f90file//':'//subr// &
               ': End of record, check file/read-in!'
          stop
       end if
       !// Get indices 1-360, starting at -180W
       ILee = Nint(lon - 0.5) + 181
       JLee = Nint(lat - 0.5) + 91
       LLee = ilev + 1 !// Lee starts at 0
       if (ILee.lt.1.or.ILee.gt.360) print*,'ILee=',ILee,'!!!'
       if (JLee.lt.1.or.JLee.gt.180) print*,'JLee=',JLee,'!!!'
       !DistLee(ILee,JLee,LLee) = Dist
       R4C_FUEL(ILee,JLee,LLee) = Fuel
       R4C_NO(ILee,JLee,LLee)   = NOx/46._r8*30._r8 !// From NO2 to NO
    End Do

    !// Close file
    close(ifnr)


    !// Interpolate horizontally over aircraft levels for each of the
    !// components.
    write(6,'(a)') f90file//':'//subr// &
         ': Interpolating horizontally to CTM grid'

    !// Sum of emissions before interpolation
    ESUM_AC(:) = 0._r8


    !// Convert input units to useful units
    !// Dist            Distance travelled in grid cell (km/year)
    !// Fuel            Fuel burnt in grid cell (kg/year)
    !// NOx             Emissions per grid cell (kg of NOx as NO2/year)
    !// CO2             Emissions per grid cell (kg/year)
    !// Soot            Emissions per grid cell(kg/year)
    !// Particles       Emissions per grid cell(number of particles per 
    !// NOx: From kg/year to kg/s
    RFACT_NOX = 1._r8/(3600._r8*24._r8*365._r8)
    !// Fuel: From kg/year to kg/s
    RFACT_FUEL = 1._r8/(3600._r8*24._r8*365._r8)
    
    !// NO emissions
    !// ------------------------------------------------------
    ICOMP = 43
    !// Find NAC
    NAC = -1
    do N = 1, EPAR_AC
       if (ECOMP_AC(N) .eq. ICOMP) then
          NAC = N
          exit
       end if
    end do

    !// Is component in ECOMP_AC?
    if (NAC .gt. 0) then

       !// Sum of original data as kg[NO]/year (already converted to NO)
       ESUM_AC(NAC) = sum(R4C_NO)

!$omp parallel private (LAC,R8CTM) &
!$omp          shared (R4C_NO,IM_AC,JM_AC,XBEDGE,YBEDGE,in_src,LM_AC) &
!$omp          shared (XDEDG,YDEDG,RFACT_NOX) &
!$omp          default(NONE)
!$omp do
       do LAC = 1, LM_AC

          !// No need to interpolate zero field
          if (maxval(R4C_NO(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
          end if
          !// Interpolate into R8CTM (no moments, only mean field
          call E_GRID(R4C_NO(:,:,LAC),XBEDGE,YBEDGE,IM_AC,JM_AC, &
               R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)

          !// Put into 3D array
          !// Convert from kg[NO]/month to kg[NO]/s
          in_src(:,:,LAC) = R8CTM(:,:)*RFACT_NOX

       end do
!$omp end do
!$omp end parallel

       do J = 1, JPAR
          do I = 1, IPAR
             do LAC = 1, actual_levels
                !// Convert from kg/month to kg/s
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
          end do
       end do
       !// Assume no hour-to-hour variation
       do H=2,24
          do J=1,JPAR
             do I=1,IPAR
                do LAC=1,actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
          end do
       end do

    end if !// if (NAC .gt. 0) then


    !// H2O emissions
    !// ------------------------------------------------------
    ICOMP = 114
    !// Find NAC
    NAC = -1
    do N = 1, EPAR_AC
       if (ECOMP_AC(N) .eq. ICOMP) then
          NAC = N
          exit
       end if
    end do

    !// Is component in ECOMP_AC?
    if (NAC .gt. 0) then

       !// Emission index of H2O from fuel [kg/kg]
       EI_H2O = 1.230_r8
       !// Sum of original data as kg[H2O]/year
       ESUM_AC(NAC) = sum(R4C_FUEL)*EI_H2O

!$omp parallel private (LAC,R8CTM) &
!$omp          shared (R4C_FUEL,IM_AC,JM_AC,XBEDGE,YBEDGE,in_src,LM_AC) &
!$omp          shared (XDEDG,YDEDG,RFACT_FUEL,EI_H2O) &
!$omp          default(NONE)
!$omp do
       do LAC = 1, LM_AC

          !// No need to interpolate zero field
          if (maxval(R4C_FUEL(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
          end if
          !// Interpolate into R8CTM (no moments, only mean field
          call E_GRID(R4C_FUEL(:,:,LAC),XBEDGE,YBEDGE,IM_AC,JM_AC, &
                  R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)

          !// Put into 3D array
          !// Convert to kg[H2O]/s
          in_src(:,:,LAC) = R8CTM(:,:)*RFACT_FUEL*EI_H2O

       end do
!$omp end do
!$omp end parallel

       !// Put into emission array, first hour
       do J = 1, JPAR
          do I = 1, IPAR
             do LAC = 1, actual_levels
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
          end do
       end do
       !// Assume no hour-to-hour variation
       do H=2,24
          do J=1,JPAR
             do I=1,IPAR
                do LAC=1,actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
          end do
       end do


       !// H2O Emissions for 114 also applies to 148
       ICOMP = 148
       !// Find NAC for 148, call it ILEE
       ILEE = -1
       do N = 1, EPAR_AC
          if (ECOMP_AC(N) .eq. ICOMP) then
             ILEE = N
             exit
          end if
       end do
       if (ILEE .gt. 0) then
          EMIS_AC_IN(:,ILEE,:,:,:) = EMIS_AC_IN(:,NAC,:,:,:)
          ESUM_AC(ILEE) = ESUM_AC(NAC)
       end if

    end if !// if (NAC .gt. 0) then


    !// De-allocate
    deallocate (XBEDGE,YBEDGE, R4C_FUEL, R4C_NO)




    !// Print total emissions
    write(6,'(a)') '    Total emissions [Tg/year] for emitted species:'
    do NAC = 1, EPAR_AC
      N = ECOMP_TRNR(NAC)
      if (N .gt. 0) then
         write(6,'(a,2f16.7)') '    '//TNAME(N)// &
              ' before/after hor. interp.:',&
              ESUM_AC(NAC)*1.e-9_r8, &
              sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
      else
         write(6,'(a,i3,a,2f16.7,a)') '    Comp ',ECOMP_AC(NAC), &
              ' before/after hor. interp.:',&
              ESUM_AC(NAC)*1.e-9_r8, &
              sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
         write(6,'(15x,a)')'*** Not transported - will not be used! ***'
      end if
    end do

    !// --------------------------------------------------------------------
  end subroutine read_react4c_original
  !// ----------------------------------------------------------------------




  !// Read original resolution data and interpolate horizontally
  !// ----------------------------------------------------------------------
  subroutine read_quantify_original()
    !// --------------------------------------------------------------------
    !// Read original resolution data from TradeOff and interpolate
    !// horizontally.
    !//
    !// Ole Amund Sovde, January 2011
    !// --------------------------------------------------------------------
    use regridding, only: E_GRID
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Locals
    integer :: I,J,H,LAC,NAC,IOS, N, ICOMP
    integer :: IM_AC, JM_AC, LM_AC, MAX_ITER, grid_type

    !// For getting file name
    character(len=160) :: FILENAME
    character(len=2) :: CMON
    !// File variables
    logical :: fnr_ok
    integer :: ifnr

    !// Interpolated field
    real(r8) :: in_src(IPAR,JPAR,maxlevels), R8CTM(IPAR,JPAR)
    real(r8) :: ESUM_AC(EPAR_AC)
    !// Lee grid: 1x1 deg horizontal, 610 m vertical.
    !// Horizontally: First grid box between 180W and 179W and
    !//               between 90S and 89N.
    !// Vertically: L=1 means interval between 0 and 610m 'above ground'
    real(r8)  :: I_in, J_in, NOx, Fuel, Dist, RFACT, EI_H2O
    integer :: ILee, JLee, LLee, L_in
    !// Travelled distance in Lee grid: km/(month*Leebox)
    !// Fuel consumption in Lee grid:   kg(Fuel)/(month*Leebox)
    !// NOx emissions in Lee grid:      kg(NO2)/(month*Leebox)
    real(r8), dimension(:),allocatable :: XBEDGE, YBEDGE
    real(r8), dimension(:,:,:),allocatable :: FuelLee, NOxLee
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_quantify_original'
    !// --------------------------------------------------------------------

    !// File name
    if (trim(AirScen) .eq. 'Quantify_MA8') then
       !// Quantify May 2008
       FILENAME = trim(AirEmisPath)//'QUANTIFY/QTF_2000_Scaled.nc'
    else
       write(6,'(a)') f90file//':'//subr// &
            ': No such scenario defined: '//trim(AirScen)
       stop
    end if


    !// Original resolution
    IM_AC = 360
    JM_AC = 180
    LM_AC = actual_levels
    MAX_ITER = IM_AC*JM_AC*LM_AC
    !// Grid starts at 0E (converted to this in read-in routine)
    grid_type = 1

    !// Allocate arrays
    allocate (XBEDGE(IM_AC+1),YBEDGE(JM_AC+1), &
         FuelLee(IM_AC,JM_AC,LM_AC), NOxLee(IM_AC,JM_AC,LM_AC))
    call get_xyedges(IM_AC,JM_AC,XBEDGE,YBEDGE,grid_type)

    !// Initialize
    FuelLee(:,:,:) = 0._r8
    NOxLee(:,:,:) = 0._r8


    !// Interpolate horizontally over aircraft levels for each of the
    !// components.
    write(6,'(a)') f90file//':'//subr// &
         ': Interpolating horizontally to CTM grid'

    !// Sum of emissions before interpolation
    ESUM_AC(:) = 0._r8

    !// Converting factor from kg/month to kg/s
    RFACT = 12._r8/(3600._r8*24._r8*365._r8)

    !// NO emissions
    !// ------------------------------------------------------
    ICOMP = 43
    !// Find NAC
    NAC = -1
    do N = 1, EPAR_AC
       if (ECOMP_AC(N) .eq. ICOMP) then
          NAC = N
          exit
       end if
    end do

    !// Is component in ECOMP_AC?
    if (NAC .gt. 0) then

       !// Read data for this month
       call readmass_qfyAIR_month(JMON,NOxLee,FILENAME,'NOx',IM_AC,JM_AC,LM_AC)
       write(6,'(a)') '    Read NOx data'

       !// Sum of original data as kg[NO]/year
       ESUM_AC(NAC) = sum(NOxLee)*12._r8*30._r8/46._r8

!$omp parallel private (LAC,R8CTM) &
!$omp          shared (NOxLee,IM_AC,JM_AC,XBEDGE,YBEDGE,in_src,LM_AC) &
!$omp          shared (XDEDG,YDEDG,RFACT) &
!$omp          default(NONE)
!$omp do
       do LAC = 1, LM_AC

          !// No need to interpolate zero field
          if (maxval(NOxLee(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
          end if
          !// Interpolate into R8CTM (no moments, only mean field
          call E_GRID(NOxLee(:,:,LAC),XBEDGE,YBEDGE,IM_AC,JM_AC, &
               R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)

          !// Put into 3D array
          !// Convert from kg[NO2]/month to kg[NO]/s
          in_src(:,:,LAC) = R8CTM(:,:)/46._r8*30._r8*RFACT

       end do
!$omp end do
!$omp end parallel

       do J = 1, JPAR
          do I = 1, IPAR
             do LAC = 1, actual_levels
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
          end do
       end do
       !// Assume no hour-to-hour variation
       do H=2,24
          do J=1,JPAR
             do I=1,IPAR
                do LAC=1,actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
          end do
       end do

    end if !// if (NAC .gt. 0) then




    !// H2O emissions
    !// ------------------------------------------------------
    ICOMP = 114
    !// Find NAC
    NAC = -1
    do N = 1, EPAR_AC
       if (ECOMP_AC(N) .eq. ICOMP) then
          NAC = N
          exit
       end if
    end do

    !// Is component in ECOMP_AC?
    if (NAC .gt. 0) then

       !// Read data for this month
       call readmass_qfyAIR_month(JMON,FuelLee,FILENAME,'Fuel',IM_AC,JM_AC,LM_AC)
       write(6,'(a)') '    Read Fuel data'

       !// Emission index of H2O from fuel [kg/kg]
       EI_H2O = 1.230_r8
       !// Sum of original data as kg[H2O]/year
       ESUM_AC(NAC) = sum(FuelLee)*12._r8*EI_H2O

!$omp parallel private (LAC,R8CTM) &
!$omp          shared (FuelLee,IM_AC,JM_AC,XBEDGE,YBEDGE,in_src,LM_AC) &
!$omp          shared (XDEDG,YDEDG,RFACT,EI_H2O) &
!$omp          default(NONE)
!$omp do
       do LAC = 1, LM_AC

          !// No need to interpolate zero field
          if (maxval(FuelLee(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
          end if
          !// Interpolate into R8CTM (no moments, only mean field
          call E_GRID(FuelLee(:,:,LAC),XBEDGE,YBEDGE,IM_AC,JM_AC, &
                  R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)

          !// Put into 3D array
          !// Convert to kg[H2O]/s
          in_src(:,:,LAC) = R8CTM(:,:)*RFACT*EI_H2O

       end do
!$omp end do
!$omp end parallel

       !// Put into emission array, first hour
       do J = 1, JPAR
          do I = 1, IPAR
             do LAC = 1, actual_levels
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
          end do
       end do
       !// Assume no hour-to-hour variation
       do H=2,24
          do J=1,JPAR
             do I=1,IPAR
                do LAC=1,actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
          end do
       end do

       !// H2O Emissions for 114 also applies to 148
       ICOMP = 148
       !// Find NAC for 148, call it ILEE
       ILEE = -1
       do N = 1, EPAR_AC
          if (ECOMP_AC(N) .eq. ICOMP) then
             ILEE = N
             exit
          end if
       end do
       if (ILEE .gt. 0) then
          ESUM_AC(ILEE) = ESUM_AC(NAC)
          EMIS_AC_IN(:,ILEE,:,:,:) = EMIS_AC_IN(:,NAC,:,:,:)
       end if


    end if !// if (NAC .gt. 0) then


    !// De-allocate
    deallocate (XBEDGE,YBEDGE, FuelLee, NOxLee)



    !// Print total emissions
    write(6,'(a)') '    Total emissions [Tg/year] for emitted species:'
    do NAC = 1, EPAR_AC
      N = ECOMP_TRNR(NAC)
      if (N .gt. 0) then
         write(6,'(a,2f16.7)') '    '//TNAME(N)// &
              ' before/after hor. interp.:',&
              ESUM_AC(NAC)*1.e-9_r8, &
              sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
      else
         write(6,'(a,i3,a,2f16.7)') '    CompID ',ECOMP_AC(NAC), &
              ' before/after hor. interp.:',&
              ESUM_AC(NAC)*1.e-9_r8, &
              sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
         write(6,'(15x,a)')'*** Not transported - will not be used! ***'
      end if
    end do



    !// --------------------------------------------------------------------
  end subroutine read_quantify_original
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  Subroutine readmass_qfyAIR_month(month,qfyAIR,filename,varchar,&
       IEMIS,JEMIS,numlevels)
    !// --------------------------------------------------------------------
    !// Reads QUANTIFY netCDF-files.
    !// --------------------------------------------------------------------
    use netcdf
    use ncutils, only: handle_err
    Implicit None


    !INPUT
    character(*), intent(in) :: filename   ! Name of netCDFfile 
    character(*), intent(in) :: varchar    ! Name of variable to fetch
    integer,intent(in)     :: IEMIS        ! Number of longitudes
    integer,intent(in)     :: JEMIS        ! Number of latitudes
    integer,intent(in)     :: numlevels    ! Number of levels
    integer,intent(in)     :: month        ! The timestep to be read

    !OUTPUT
    real(r8),intent(out) :: qfyAIR(IEMIS,JEMIS,numlevels)   !Three dimensional field

    !LOCAL
    !LOCAL NETCDF DIMENSION IDs ETCETERA
    real(r4) :: field(IEMIS,JEMIS,numlevels) !Four dimensional field
    integer                  :: lon_dim_id      !Id for longitude dimension
    integer                  :: lon_id          !Id for variable longitude
    real(r8)                 :: lon(IEMIS)      !variable lon (in file)
    integer                  :: lat_dim_id      !Id for latitude dimension
    integer                  :: lat_id          !Id for latitude
    real(r8)                 :: lat(JEMIS)      !Variable for latitude
    integer                  :: time_dim_id     !Id for time dimension
    integer                  :: time_id         !Id for time
    integer                  :: lev_dim_id      !Id for time dimension
    integer                  :: lev_id          !Id for time
    real(r8)                 :: lev(numlevels)  !Variable for latitude
    integer                  :: field_dim_id(4) !Dimension id for field
    integer                  :: field_id        !Variable id for field
    integer                  :: srt_lon_lat_lev_time(4) !Start point 
    integer                  :: cnt_lon_lat_lev_time(4) !Count indexes
    integer                  :: nlons           !Longitudes in file
    integer                  :: nlats           !Latitudes in file
    integer                  :: nsteps          !Timesteps avaiable in file
    integer                  :: nlevs           !levels avaiable in file
    integer                  :: status          !status of process (0=OK)
    integer                  :: ncid            !file id 

    !// Values declarations
    integer :: igrid,jgrid
    !// Div.
    Integer :: I,J,L
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'readmass_qfyAIR_month'
    !// --------------------------------------------------------------------

    !Array which tells you where to start picking your 3D field
    srt_lon_lat_lev_time= (/ 1 , 1 , 1 , month /)    !Start array
    !Array which tells you how far to count when picking it
    cnt_lon_lat_lev_time= (/ IEMIS , JEMIS , numlevels , 1 /)

    status=nf90_noerr  !Status is 0 and should be kept that way !!

    !// Open the netCDF file for reading
    write(6,'(a)') '    Reading: '//trim(filename)
    status=nf90_open(filename, nf90_nowrite, ncid)
    if(status/=nf90_noerr)call handle_err(status)

    !Inquire dimension ids
    status = nf90_inq_dimid(ncid,"time",time_dim_id)
    if(status/=nf90_noerr)call handle_err(status)

    status = nf90_inq_dimid(ncid,"lon",lon_dim_id)
    if(status/=nf90_noerr)call handle_err(status)

    status = nf90_inq_dimid(ncid,"lat",lat_dim_id)
    if(status/=nf90_noerr)call handle_err(status)

    status = nf90_inq_dimid(ncid,"lev",lev_dim_id)
    if(status/=nf90_noerr)call handle_err(status)

    !Dimension id for 3D field /lon/lat/lev/time
    field_dim_id(1)=lon_dim_id
    field_dim_id(2)=lat_dim_id
    field_dim_id(3)=lev_dim_id
    field_dim_id(4)=time_dim_id

    !Inquire dimensions
    status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)

    if(status/=nf90_noerr)call handle_err(status)
    if(nlats/=JEMIS)then
       write(6,'(a)') f90file//':'//subr//': ERROR'
       write(6,*)'file '//trim(filename)//' reports lat = ',nlats
       write(6,*)'your array has dimension ',JEMIS
       stop
    end if
    status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
    if(status/=nf90_noerr)call handle_err(status)
    if(nlons/=IEMIS)then
       write(6,'(a)') f90file//':'//subr//': ERROR'
       write(6,*)'file '//trim(filename)//' reports lon = ',nlons
       write(6,*)'your array has dimension',IEMIS
       stop
    end if
    status = nf90_Inquire_Dimension(ncid,lev_dim_id,len=nlevs)
    if(status/=nf90_noerr)call handle_err(status)
    if(numlevels.gt.nlevs.or.nlevs.le.0)then
       write(6,'(a)') f90file//':'//subr//': ERROR'
       write(6,*)'file '//trim(filename)//' reports nlevs = ',nlevs
       write(6,*)'you try to read levels',numlevels
       stop
    end if
    status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
    if(status/=nf90_noerr)call handle_err(status)
    if(month.gt.nsteps.or.nsteps.le.0)then
       write(6,'(a)') f90file//':'//subr//': ERROR'
       write(6,*)'file '//trim(filename)//' reports nsteps = ',nsteps
       write(6,*)'you try to read timestep',month
       stop
    end if

    !Get variable ID
    status=nf90_inq_varid(ncid,varchar,field_id) 
    if(status/=nf90_noerr)call handle_err(status)

    !Finally after all this, you can get the variable you want !!
    !and put it in the threeDfield array
    status=nf90_get_var(ncid,field_id,field, &
         start=srt_lon_lat_lev_time,   &
         count=cnt_lon_lat_lev_time )
    if(status/=nf90_noerr)call handle_err(status)

    Do L=1,numlevels
       Do J=1,JEMIS
          Do I=181,360
             igrid=I-180
             !qfyAIR(igrid,J,L)=field(i,J,L,month)
             qfyAIR(igrid,J,L) = field(i,J,L)
          End Do
          Do I=1,180
             igrid=I+180
             !qfyAIR(igrid,J,L)=field(i,J,L,month)
             qfyAIR(igrid,J,L)=field(i,J,L)
          End Do
       End Do
    End Do

    !Closing file
    status=nf90_close(ncid)
    if(status/=nf90_noerr)call handle_err(status)

    !// --------------------------------------------------------------------
  end subroutine readmass_qfyAIR_month
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine read_ceds_original()
    !// --------------------------------------------------------------------
    !// Read original resolution data from CEDS (2016) and interpolate
    !// horizontally.
    !//
    !// MTL 07/01/17
    !// --------------------------------------------------------------------
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use ncutils, only: get_netcdf_var_1d, readnc_3d_from4d
    use regridding, only: E_GRID
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Locals
    integer :: I,J,H,LAC,NAC,N, ICOMP
    integer :: MAX_ITER, grid_type

    integer :: sTime, M, getY, start_year, daynr
    integer :: nLon, nLat, nLev, nTime

    !// For getting file name
    character(len=200) :: infile
    real(r8), dimension(EPAR_AC):: scalefac
    character(len=5), dimension(EPAR_AC) :: FILECOMP
    character(len=5) :: version
    character(len=13) :: infileyear

    !// File variables
    logical :: fnr_ok
    integer :: ifnr

    !// Interpolated field
    real(r8) :: in_src(IPAR,JPAR,maxlevels), R8CTM(IPAR,JPAR)
   
    real(r8) :: ESUM_AC(EPAR_AC), sumB, sumA

    real(r8), dimension(:), allocatable :: &
         XBEDGE, YBEDGE, XYBOX,inTime,inLon,inLat,inLev
    real(r8), dimension(:,:,:), allocatable :: HLFDUM, RDUM

    integer, dimension(12), parameter :: midmonth = &
         (/15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349 /)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_ceds_original'
    !// --------------------------------------------------------------------

    !// Generate file name
    getY = AirScenYear  
    start_year = 1750

    if (getY .ge. 1750 .and. getY .le. 1799) then 
       version = '06-18'
       infileyear = '175001-179912'
    else if (getY .ge. 1800 .and. getY .le. 1849) then
       version = '06-18'
       infileyear = '180001-184912'
    else if (getY .eq. 1850) then
       version = '06-18'
       infileyear = '185001-185012'
    else if (getY .ge. 1851 .and. getY .le. 1899) then
       version = '07-26'
       infileyear = '185101-189912'
    else if (getY .ge. 1900 .and. getY .le. 1949) then
       version = '07-26'
       infileyear = '190001-194912'
    else if (getY .ge. 1950 .and. getY .le. 1999) then
       version = '07-26'
       infileyear = '195001-199912'
    else if (getY .ge. 2000 .and. getY .le. 2014) then
       version = '07-26'
       infileyear = '200001-201412'
    endif
    
    if (trim(AirScen) .eq. 'CEDS') then
       write(6,'(a)') f90file//':'//subr//': Reading CEDS aircraft emissions'
    else
       write(6,'(a)') f90file//':'//subr//': No such CEDS scenario defined.'
       stop
    endif


    !// Set file component name and scale factor
    do NAC = 1, EPAR_AC
       if (trim(ECOMP_NAMES(NAC)) .eq. 'CO') then
          FILECOMP(NAC) = 'CO'
          scalefac(NAC) = 1._r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C2H4') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.1546_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C2H6') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0052_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C3H6') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0453_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C6H14') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.1812_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C6HXR_SOA' .or. &
                trim(ECOMP_NAMES(NAC)) .eq. 'C6HXR') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0115_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'CH2O') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.123_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'CH3CHO') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.2195_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'NO') then
          FILECOMP(NAC) = 'NOx'
          scalefac(NAC) = 0.652_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C3H8') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.000078_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'ACETONE') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.00369_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'SO2') then
          FILECOMP(NAC) = 'SO2'
          scalefac(NAC) = 1._r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'Tolmatic') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.04769_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'Benzene') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0168_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'omFF1fob') then
          FILECOMP(NAC) = 'OC'
          scalefac(NAC) = 0.5_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'omFF1fil') then
          FILECOMP(NAC) = 'OC'
          scalefac(NAC) = 0.5_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'bcFF1fob') then
          FILECOMP(NAC) = 'BC'
          scalefac(NAC) = 0.8_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'bcFF1fil') then
          FILECOMP(NAC) = 'BC'
          scalefac(NAC) = 0.2_r8
       else
          write(6,'(a)') f90file//':'//subr// &
               ': Unknown species: '//trim(ECOMP_NAMES(NAC))
          stop
       end if
    end do


    !//Get resolution (latitude/longitude/time)
    infile = trim(AirEmisPath)//trim(FILECOMP(1))// &
         '-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-v2016-'// &
         trim(version)//'_gr_'//trim(infileyear)//'.nc'

    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )
    call get_netcdf_var_1d( infile, 'level', inLev  )
    call get_netcdf_var_1d( infile, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nLev  = SIZE( inLev  )
    nTime = SIZE( inTime )

    if (nLev .ne. actual_levels) then
       write(6,'(a,2i5)') f90file//':'//subr// &
            ': nLev != actual_levels ',nLev,actual_levels
       stop
    end if

    !// Time step to fetch
    !// Day of mid-month, since 1750:
    daynr = midmonth(JMON) + (getY - start_year)*365
    sTime = -1
    do M = 1, nTime
       if (inTime(M) .eq. daynr) then
          sTime = M !// Found the time step
          exit
       end if
    end do

    if (sTime .eq. -1) then
       write(6,'(a)') f90file//':'//subr//': Wrong inTime index'
       stop
    else
       write(6,'(a,4i5,f12.1)') f90file//':'//subr// &
            ': mon/year/sTime/daynr:',JMON,getY,sTime,daynr,inTime(sTime)
    end if

    !// Allocate arrays
    allocate (XBEDGE(nLon+1),YBEDGE(nLat+1), &
              HLFDUM(nLon,nLat,nLev)&
              ,RDUM(nLon,nLat,nLev),XYBOX(nLat))

    !// Data starts at (180W,90S), i.e. grid_type 2
    grid_type = 2
    call get_xyedges(nLon,nLat,XBEDGE,YBEDGE,grid_type)

    !// Grid box areas
    do J = 1, nLat
       XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
            * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do

    !// Sum of emissions before interpolation
    ESUM_AC(:) = 0._r8

    do NAC = 1, EPAR_AC

       !// Cycle if species is not included
       if (ECOMP_AC(NAC) .le. 0) cycle

       !// Initialize
       RDUM(:,:,:) = 0._r8
       HLFDUM(:,:,:) = 0._r8
       in_src(:,:,:) = 0._r8

       !// Emissions
       !// ------------------------------------------------------

       infile = trim(AirEmisPath)//trim(FILECOMP(NAC))// &
            '-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-v2016-'// &
            trim(version)//'_gr_'//trim(infileyear)//'.nc'

       !// Read data for this month
       call readnc_3d_from4d(INFILE, 'lon', nLon, 'lat', nLat, &
            'level', nLev, 'time', sTime, 'AIR', RDUM)
       write(*,'(a)') '    Read '//FILECOMP(NAC)//' data'


       !// Scale with corresponding fraction
       HLFDUM(:,:,:) = RDUM(:,:,:) * scalefac(NAC)

       !// Scale with area (xybox is m2)--> kg/s. Done before interpolation,
       !// since interpolation routine needs field to be per grid box.
       do LAC = 1, nLev
          do J = 1, nLat
             HLFDUM(:,J,LAC) = HLFDUM(:,J,LAC) * XYBOX(J)
          end do
       end do
      
         
          
!$omp parallel private (LAC,R8CTM) &
!$omp          shared (HLFDUM,nLon,nLat,XBEDGE,YBEDGE,in_src,nLev) &
!$omp          shared (XDEDG,YDEDG) &
!$omp          shared (XYBOX) &
!$omp          default(NONE)
!$omp do
       do LAC = 1, nLev

          !// No need to interpolate zero field
          if (maxval(HLFDUM(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
          end if

          !// Interpolate into R8CTM (no moments, only mean field
          call E_GRID(HLFDUM(:,:,LAC),XBEDGE,YBEDGE,nLon,nLat, &
               R8CTM,XDEDG,YDEDG,IPAR,JPAR,2)

          !// Put into 3D array
          in_src(:,:,LAC) = R8CTM(:,:)

       end do
!$omp end do
!$omp end parallel

       !write(*,'(a,es20.12,es20.12,es20.12)'), 'Max/min in_src kg/s:',maxval(in_src),minval(in_src),sum(in_src)

       !// Put into emission array, first hour
       do J = 1, JPAR
          do I = 1, IPAR
             do LAC = 1, actual_levels
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
          end do
       end do
       !// Assume no hour-to-hour variation
       do H = 2, 24
          do J = 1, JPAR
             do I = 1, IPAR
                do LAC = 1, actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
          end do
       end do
          
       !// Sum of original data as kg/year
       ESUM_AC(NAC) = sum(HLFDUM) * 86400._r8 * 365._r8
    end do !// do NAC = 1, EPAR_AC



    !// De-allocate
    deallocate (XBEDGE,YBEDGE,XYBOX,HLFDUM,RDUM)

    !// Print total emissions
    write(6,'(a)') '    Total emissions [Tg/year] for emitted species:'
    do NAC = 1, EPAR_AC
       N = ECOMP_TRNR(NAC)
       if (N .gt. 0) then
          sumB = ESUM_AC(NAC)*1.e-9_r8
          sumA = sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
          if (abs(sumA-sumB)/sumB .gt. 1.e-5_r8) then
             write(6,'(a,es20.12,es20.12)') f90file//':'//subr// &
                  ':  Wrong sum before/after interp. ', sumB, sumA
             stop
          else
             write(6,'(a,2f16.7)') '    '//TNAME(N)//' emitted (Tg/yr): ',sumB
          end if
       else
          write(6,'(a)') '    '//ECOMP_NAMES(NAC)//' not included'
       end if
    end do


    !// --------------------------------------------------------------------
  end subroutine read_ceds_original
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_ceds_2017()
    !// --------------------------------------------------------------------
    !// Read original resolution data from CEDS and interpolate
    !// horizontally.
    !//
    !// Amund Sovde Haslerud, October 2017
    !// MTL 07/01/17
    !// --------------------------------------------------------------------
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use ncutils, only: get_netcdf_var_1d, readnc_3d_from4d
    use regridding, only: E_GRID
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Locals
    integer :: I,J,H,LAC,NAC,N, ICOMP
    integer :: MAX_ITER, grid_type

    integer :: sTime, M, getY, start_year, daynr
    integer :: nLon, nLat, nLev, nTime

    !// For getting file name
    character(len=200) :: infile
    real(r8), dimension(EPAR_AC):: scalefac
    character(len=5), dimension(EPAR_AC) :: FILECOMP
    character(len=10) :: version
    character(len=4) :: ctag
    character(len=13) :: infileyear

    !// File variables
    logical :: fnr_ok
    integer :: ifnr

    !// Interpolated field
    real(r8) :: in_src(IPAR,JPAR,maxlevels), R8CTM(IPAR,JPAR)
   
    real(r8) :: ESUM_AC(EPAR_AC), sumB, sumA

    real(r8), dimension(:), allocatable :: &
         XBEDGE, YBEDGE, XYBOX,inTime,inLon,inLat,inLev
    real(r8), dimension(:,:,:), allocatable :: HLFDUM, RDUM

    integer, dimension(12), parameter :: midmonth = &
         (/15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349 /)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_ceds_2017'
    !// --------------------------------------------------------------------

    !// Generate file name
    getY = AirScenYear  
    start_year = 1750

    version = '2017-05-18'
    ctag = '_gn_'
    if (getY .ge. 1750 .and. getY .le. 1799) then 
       infileyear = '175001-179912'
    else if (getY .ge. 1800 .and. getY .le. 1849) then
       infileyear = '180001-184912'
    else if (getY .eq. 1850) then
       infileyear = '185001-185012'
    else if (getY .ge. 1851 .and. getY .le. 1899) then
       infileyear = '185101-189912'
    else if (getY .ge. 1900 .and. getY .le. 1949) then
       infileyear = '190001-194912'
    else if (getY .ge. 1950 .and. getY .le. 1999) then
       infileyear = '195001-199912'
    else if (getY .ge. 2000 .and. getY .le. 2014) then
       infileyear = '200001-201412'
    endif
    

    !// Set file component name and scale factor
    do NAC = 1, EPAR_AC
       if (trim(ECOMP_NAMES(NAC)) .eq. 'CO') then
          FILECOMP(NAC) = 'CO'
          scalefac(NAC) = 1._r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C2H4') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.1546_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C2H6') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0052_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C3H6') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0453_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C6H14') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.1812_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C6HXR_SOA' .or. &
                trim(ECOMP_NAMES(NAC)) .eq. 'C6HXR') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0115_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'CH2O') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.123_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'CH3CHO') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.2195_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'NO') then
          FILECOMP(NAC) = 'NOx'
          scalefac(NAC) = 0.652_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'C3H8') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.000078_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'ACETONE') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.00369_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'SO2') then
          FILECOMP(NAC) = 'SO2'
          scalefac(NAC) = 1._r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'Tolmatic') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.04769_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'Benzene') then
          FILECOMP(NAC) = 'NMVOC'
          scalefac(NAC) = 0.0168_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'omFF1fob') then
          FILECOMP(NAC) = 'OC'
          scalefac(NAC) = 0.5_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'omFF1fil') then
          FILECOMP(NAC) = 'OC'
          scalefac(NAC) = 0.5_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'bcFF1fob') then
          FILECOMP(NAC) = 'BC'
          scalefac(NAC) = 0.8_r8
       else if (trim(ECOMP_NAMES(NAC)) .eq. 'bcFF1fil') then
          FILECOMP(NAC) = 'BC'
          scalefac(NAC) = 0.2_r8
       else
          write(6,'(a)') f90file//':'//subr// &
               ': Unknown species: '//trim(ECOMP_NAMES(NAC))
          stop
       end if
    end do


    !//Get resolution (latitude/longitude/time)
    infile = trim(AirEmisPath)//trim(FILECOMP(1))// &
         '-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-'// &
         trim(version)//ctag//trim(infileyear)//'.nc'

    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )
    call get_netcdf_var_1d( infile, 'level', inLev  )
    call get_netcdf_var_1d( infile, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nLev  = SIZE( inLev  )
    nTime = SIZE( inTime )

    if (nLev .ne. actual_levels) then
       write(6,'(a,2i5)') f90file//':'//subr// &
            ': nLev != actual_levels ',nLev,actual_levels
       stop
    end if

    !// Time step to fetch
    !// Day of mid-month, since 1750:
    daynr = midmonth(JMON) + (getY - start_year)*365
    sTime = -1
    do M = 1, nTime
       if (inTime(M) .eq. daynr) then
          sTime = M !// Found the time step
          exit
       end if
    end do

    if (sTime .eq. -1) then
       write(6,'(a)') f90file//':'//subr//': Wrong inTime index'
       stop
    else
       write(6,'(a,4i5,f12.1)') f90file//':'//subr// &
            ': mon/year/sTime/daynr:',JMON,getY,sTime,daynr,inTime(sTime)
    end if

    !// Allocate arrays
    allocate (XBEDGE(nLon+1),YBEDGE(nLat+1), &
              HLFDUM(nLon,nLat,nLev)&
              ,RDUM(nLon,nLat,nLev),XYBOX(nLat))

    !// Data starts at (180W,90S), i.e. grid_type 2
    grid_type = 2
    call get_xyedges(nLon,nLat,XBEDGE,YBEDGE,grid_type)

    !// Grid box areas
    do J = 1, nLat
       XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
            * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do

    !// Sum of emissions before interpolation
    ESUM_AC(:) = 0._r8

    do NAC = 1, EPAR_AC

       !// Cycle if species is not included
       if (ECOMP_AC(NAC) .le. 0) cycle

       !// Initialize
       RDUM(:,:,:) = 0._r8
       HLFDUM(:,:,:) = 0._r8
       in_src(:,:,:) = 0._r8

       !// Emissions
       !// ------------------------------------------------------

       infile = trim(AirEmisPath)//trim(FILECOMP(NAC))// &
            '-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-'// &
            trim(version)//ctag//trim(infileyear)//'.nc'

       !// Read data for this month
       call readnc_3d_from4d(INFILE, 'lon', nLon, 'lat', nLat, &
            'level', nLev, 'time', sTime, &
            trim(FILECOMP(NAC))//'_em_AIR_anthro', RDUM)
       write(*,'(a)') '    Read '//FILECOMP(NAC)//' data -> '//trim(ECOMP_NAMES(NAC))


       !// Scale with corresponding fraction
       HLFDUM(:,:,:) = RDUM(:,:,:) * scalefac(NAC)

       !// Scale with area (xybox is m2)--> kg/s. Done before interpolation,
       !// since interpolation routine needs field to be per grid box.
       do LAC = 1, nLev
          do J = 1, nLat
             HLFDUM(:,J,LAC) = HLFDUM(:,J,LAC) * XYBOX(J)
          end do
       end do
      
         
          
!$omp parallel private (LAC,R8CTM) &
!$omp          shared (HLFDUM,nLon,nLat,XBEDGE,YBEDGE,in_src,nLev) &
!$omp          shared (XDEDG,YDEDG) &
!$omp          shared (XYBOX) &
!$omp          default(NONE)
!$omp do
       do LAC = 1, nLev

          !// No need to interpolate zero field
          if (maxval(HLFDUM(:,:,LAC)) .eq. 0._r8) then
             in_src(:,:,LAC) = 0._r8
             cycle
          end if

          !// Interpolate into R8CTM (no moments, only mean field
          call E_GRID(HLFDUM(:,:,LAC),XBEDGE,YBEDGE,nLon,nLat, &
               R8CTM,XDEDG,YDEDG,IPAR,JPAR,2)

          !// Put into 3D array
          in_src(:,:,LAC) = R8CTM(:,:)

       end do
!$omp end do
!$omp end parallel

       !write(*,'(a,es20.12,es20.12,es20.12)'), 'Max/min in_src kg/s:',maxval(in_src),minval(in_src),sum(in_src)

       !// Put into emission array, first hour
       do J = 1, JPAR
          do I = 1, IPAR
             do LAC = 1, actual_levels
                EMIS_AC_IN(LAC,NAC,I,J,1) = in_src(I,J,LAC)
             end do
          end do
       end do
       !// Assume no hour-to-hour variation
       do H = 2, 24
          do J = 1, JPAR
             do I = 1, IPAR
                do LAC = 1, actual_levels
                   EMIS_AC_IN(LAC,NAC,I,J,H) = EMIS_AC_IN(LAC,NAC,I,J,1)
                end do
             end do
          end do
       end do
          
       !// Sum of original data as kg/year
       ESUM_AC(NAC) = sum(HLFDUM) * 86400._r8 * 365._r8
    end do !// do NAC = 1, EPAR_AC



    !// De-allocate
    deallocate (XBEDGE,YBEDGE,XYBOX,HLFDUM,RDUM)

    !// Print total emissions
    write(6,'(a)') '    Total emissions [Tg/year] for emitted species:'
    do NAC = 1, EPAR_AC
       N = ECOMP_TRNR(NAC)
       if (N .gt. 0) then
          sumB = ESUM_AC(NAC)*1.e-9_r8
          sumA = sum(EMIS_AC_IN(:,NAC,:,:,:))*3600._r8*365.e-9_r8
          if (abs(sumA-sumB)/sumB .gt. 1.e-5_r8) then
             write(6,'(a,es20.12,es20.12)') f90file//':'//subr// &
                  ':  Wrong sum before/after interp. ', sumB, sumA
             stop
          else
             write(6,'(a,2f16.7)') '    '//TNAME(N)//' emitted (Tg/yr): ',sumB
          end if
       else
          write(6,'(a)') '    '//ECOMP_NAMES(NAC)//' not included'
       end if
    end do


    !// --------------------------------------------------------------------
  end subroutine read_ceds_2017
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_xyedges(ISIZE,JSIZE,XBEDGE,YBEDGE,START)
    !// --------------------------------------------------------------------
    !// Set up edges for a certain resolution.
    !// START: 1 = Lower-left edge at (0,90S)
    !//        2= lower-left edges at (180W, 90S)
    !//
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: ISIZE,JSIZE
    integer, intent(in) :: START
    !// Output
    real(r8), intent(out) :: XBEDGE(ISIZE+1), YBEDGE(JSIZE+1)

    !// Locals
    real(r8) :: ioffset, joffset, factor
    integer :: I, J
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'get_xyedges'
    !// --------------------------------------------------------------------

    !// Check resolution
    if (ISIZE.eq.360) then
       
       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 1._r8
          joffset = 91._r8
          factor = 1._r8
       else if (START .eq. 2) then
          !// Standard 1x1 grid with (1,1)-box
          !// with lower-left edges at (180W, 90S)
          ioffset = 181._r8
          joffset = 91._r8
          factor = 1._r8
       else
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Cannot calculate XBEDGE and YBEDGE: Unknown map',START
          stop
       end if

    else if (ISIZE.eq.720) then

       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 0.5_r8
          joffset = 90.5_r8
          factor = 0.5_r8
       else if (START .eq. 2) then
          !// Standard 0.5x0.5 grid with (1,1)-box
          !// with lower-left edges at (180W, 90S)
          ioffset = 180.5_r8
          joffset = 90.5_r8
          factor = 0.5_r8
       else
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Cannot calculate XBEDGE and YBEDGE: Unknown map',START
          stop
       end if

    else if (ISIZE.eq.1440) then

       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 0.25_r8
          joffset = 90.25_r8
          factor = 0.25_r8
       else if (START .eq. 2) then
          !// Standard 0.25x0.25 grid with (1,1)-box
          !// with lower-left edges at (180W, 90S)
          ioffset = 180.25_r8
          joffset = 90.25_r8
          factor = 0.25_r8
       else
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Cannot calculate XBEDGE and YBEDGE: Unknown map',START
          stop
       end if

    else
       write(6,'(a,2i7)') f90file//':'//subr// &
            ': Wrong resolution in emissions; get_xyedges',ISIZE, JSIZE
       stop
    end if

    do I = 1, ISIZE+1
       XBEDGE(I) = (real(I,r8)*factor - ioffset)
    end do
    do J = 1, JSIZE+1
       YBEDGE(J) = (real(J,r8)*factor - joffset)
    end do

    !// --------------------------------------------------------------------
  end subroutine get_xyedges
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine aircraft_zero_400hpa(DT)
    !// --------------------------------------------------------------------
    !// Assume aircraft H2O to be removed/washed out below 400hPa,
    !// according to Danilin et al, 1998. Instead of setting to zero,
    !// use 1hr lifetime. This will hinder creating weird gradients that
    !// may give negatives after transport (due to moments).
    !//
    !// This is carried out outside parallel region, before updating
    !// aircraft emissions.
    !//
    !// Ole Amund Sovde, February 2011
    !// --------------------------------------------------------------------
    use cmn_ctm, only: STT, SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW
    use cmn_diag, only: STTBCK
    use cmn_oslo, only: LMTROP, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: DT

    !// Locals
    integer :: I,J,L,LMAX, TNR
    real(r8) :: ptop, fact
    !// Inverse lifetime H2O below 400hPa: assume 1 hour
    real(r8), parameter :: ZTAU = 1._r8/3600._r8
    !// --------------------------------------------------------------------
    if (trsp_idx(148) .le. 0) return

    fact = exp(-DT*ZTAU)
    LMAX = maxval(LMTROP)
    TNR = trsp_idx(148)
    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LMAX

             !// Pressure at grid box top
             ptop = etaa(L+1) + etab(L+1) * P(I,J)
             if (ptop .ge. 400._r8) then
                STT(I,J,L,TNR) = max(1.e-20_r8,STT(I,J,L,TNR)*fact)
                STTBCK(I,J,L,TNR) = STT(I,J,L,TNR)
                !// Also adjust moments
                SUT(I,J,L,TNR) = SUT(I,J,L,TNR)*fact
                SVT(I,J,L,TNR) = SVT(I,J,L,TNR)*fact
                SWT(I,J,L,TNR) = SWT(I,J,L,TNR)*fact
                SUU(I,J,L,TNR) = SUU(I,J,L,TNR)*fact
                SVV(I,J,L,TNR) = SVV(I,J,L,TNR)*fact
                SWW(I,J,L,TNR) = SWW(I,J,L,TNR)*fact
                SUV(I,J,L,TNR) = SUV(I,J,L,TNR)*fact
                SUW(I,J,L,TNR) = SUW(I,J,L,TNR)*fact
                SVW(I,J,L,TNR) = SVW(I,J,L,TNR)*fact
             else
                !// Exit L-loop
                exit
             end if

          end do
       end do
    end do
    !// --------------------------------------------------------------------
  end subroutine aircraft_zero_400hpa
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine aircraft_h2o_zero()
    !// --------------------------------------------------------------------
    !// Set aircraft H2O to zero.
    !//
    !// Ole Amund Sovde, March 2013
    !// --------------------------------------------------------------------
    use cmn_ctm, only: STT, SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW
    use cmn_diag, only: STTBCK
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'aircraft_h2o_zero'
    !// --------------------------------------------------------------------
    !// Set aircraft H2O to sero
    if (trsp_idx(148) .gt. 0) then
       write(6,'(a)') f90file//':'//subr//': Setting Aircraft H2O to zero'
       STT(:,:,:,trsp_idx(148)) = 0._r8
       SUT(:,:,:,trsp_idx(148)) = 0._r8
       SVT(:,:,:,trsp_idx(148)) = 0._r8
       SWT(:,:,:,trsp_idx(148)) = 0._r8
       SUU(:,:,:,trsp_idx(148)) = 0._r8
       SVV(:,:,:,trsp_idx(148)) = 0._r8
       SWW(:,:,:,trsp_idx(148)) = 0._r8
       SUV(:,:,:,trsp_idx(148)) = 0._r8
       SUW(:,:,:,trsp_idx(148)) = 0._r8
       SVW(:,:,:,trsp_idx(148)) = 0._r8
    end if
    !// --------------------------------------------------------------------
  end subroutine aircraft_h2o_zero
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
end module emissions_aircraft
!//=========================================================================
