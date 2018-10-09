!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Mineral dust.
!//=========================================================================
module dust_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: dust_oslo
  !// DESCRIPTION: Oslo CTM3 - Module for DUST package.
  !//
  !// Contains DUST variables and the following routines:
  !//   - subroutine dust_init
  !//   - subroutine dustinput_met
  !//   - subroutine SWradbdg_get
  !//   - subroutine psadjust
  !//   - subroutine dust_globalupdate
  !//   - subroutine oro_set
  !//   - subroutine plevs0
  !//   - subroutine plevs00
  !//   - subroutine dust_master
  !//   - subroutine dust_column
  !//   - subroutine dust_set_ssrd
  !//   - subroutine dust_set_strd
  !//   - subroutine dust_set_SWVL1
  !//   - subroutine dust_set_mbl_name_ff
  !//   - subroutine dustbdg2file
  !//
  !// Ole Amund Sovde, March 2016, October 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, LWEPAR, IDBLK, JDBLK, MPBLK, NPAR_DUST
  use cmn_parameters, only: G0, CPI, A0
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Number of DUST components
  integer, parameter :: Ndustbins = NPAR_DUST

  !// TIMESTEP STUFF
  logical :: Lfirstdust  !Logical (true if dust has not been called before)
  real(r8) :: calday

  !// List of all transport numbers used for DUST components
  integer :: dust_trsp_idx(Ndustbins)
  integer :: dust_idx_tot !// Total numbers of components in use

  !// BUDGET STUFF
  !// 2D fields for burden and budget processes; all dust tracers
  real(r8), dimension(Ndustbins,ipar,jpar) :: &
       dust_production, &
       dust_drydep, &
       dust_lscv, &
       dust_cscv, &
       dust_burden
  !// Flag to turn on 3D diagnostics. If false, only total budgets
  !// will be put out to budget file.
  logical, parameter :: LDUSTDIAG3D = .false.

  !// Time span added to budget average.
  real(r8) :: dustbudgetacctime

  !// File for averages
  character(len=80), parameter :: dustbdgfile='ctm3dstbudget.nc'
  !// Number of averages written to file
  integer :: dustaverages_written

  !// Rain last 24h in a grid square (mm)
  real(r8) :: rainlast24h(8,IDBLK,JDBLK,MPBLK)
  integer :: time_since_rain(IDBLK,JDBLK,MPBLK) !Time since significant rain
                                                !in a grid square (seconds)
  integer :: dryuptime(IDBLK,JDBLK,MPBLK)       !Time needed to dry up a give
                                                !grid squard (seconds)
  real(r8) :: newtime, &   !Time needed to try up rain in last 24h (seconds)
             timeleft, &   !Time left before soil is considered dry (seconds)
             rainingrid    !Sum of Rain in grid (mm) last 24h
  real(r8), dimension(JPAR) :: latwts
  !// Limit of rain (mm) in 24h before emission of dust is stopped
  real(r8), parameter :: rainlimit = 0.5_r8

  !// Mobilisation map name
  character(len=80) :: mbl_name = 'NOT_SET'
  !// Mobilisation fudge factor matching map name
  real(r8) :: mbl_fudgefactor = 0._r8

  !// Meteorological fields used in dust calculations
  real(r8), dimension(IDBLK,JDBLK,MPBLK) :: &
             SSRD, &      !Surface solar radiation downwards (W/m2)
             !SSRU, &      !Surface solar radiation upwards (W/m2)
             !             !Skip; it is calculated as SSRD*Albedo
             SFP, &       !Surface pressure (Pa)
             STRD, &      !Surface termal (longwave) radiation downwards
             SWVL1, &     !Volumetric water content (m3/m3)
             SNOW         !Snow depth
  !// Flags set by subroutines if metdata is used for some variables
  logical :: &
       LSSRD = .false., & !It is set to T if metdata is used for SSRD
       LSTRD = .false., & !It is set to T if metdata is used for STRD
       LSWVL1 = .false.   !It is set to T if metdata is used for SWVL1

  !// To diagnose prod/loss/burden
  real(r8) :: lastburden, lastlscv, lastcscv, lastprod, lastddep, lastchem
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'dust_oslo.f90'
  !// ----------------------------------------------------------------------
  !// All variables are to be saved.
  save
  !// Keep most of the module private
  private
  !// And specify what will be public
  public Ndustbins, dust_trsp_idx, &
       dust_init, dust_globalupdate, dust_master, dustbdg2file, &
       dustInstBdg, dust_set_ssrd, dust_set_strd, dust_set_SWVL1, &
       dust_set_mbl_name_ff
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine dust_init(NDAY,DT_DUST)
    !// --------------------------------------------------------------------
    !// Initialize DUST simulations. Only figures out the transport numbers
    !// and indices for DUST components.
    !//
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR
    use cmn_ctm, only: YGRD
    use cmn_chem, only: TNAME
    use cmn_parameters, only: R_AIR
    use grid, only: GAUSST2
    use pmgrid, only: PLON, PLAT, PLEV, PCNST, dst_nbr
    use dead_inirun, only: inirun
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY
    real(r8), intent(in) :: DT_DUST
    !// Local variables
    integer :: N
    character(len=11) :: tracername !// Name of tracer

    real(r8) :: logbin1, &    !// Logarithm of lowest bin length (log(m))
              logbinmax, &  !// Logarithm of highest bin length (log(m))
              totlogstep, & !// Logarithm of whole bin length
              logstep, &    !// Logarithm of one bin length
              logbin        !// Logarithm of arbitrary bin length
    real(r8), dimension(JPAR) :: plat_rdn
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'dust_init'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Initialising DUST'

    !// Check number of DUST species; CTM3 has one parameter for it (NPAR_DUST)
    !// DEAD has another (DST_NBR)
    if (NPAR_DUST .ne. DST_NBR) then
       print*,'*** DUST: Inconsistencies between CTM3 and DEAD'
       print*,'    NPAR_DUST (CTM3 params.h):', NPAR_DUST
       print*,'    DST_NBR (DUST params.h):  ', DST_NBR
       print*,'    Stopping!'
       stop
    end if
    if (PCNST .ne. DST_NBR) then
       print*,'*** DUST: Inconsistencies in DEAD'
       print*,'    PCNST   (DUST params.h): ', PCNST
       print*,'    DST_NBR (DUST params.h): ', DST_NBR
       print*,'    Stopping!'
       stop
    end if

    !// Check resolution
    if (IPAR.ne.PLON .or. JPAR.ne.PLAT .or.LPAR.ne.PLEV) then
       print*,'*** DUST: Inconsistencies between CTM3 and DEAD'
       print*,'    IPAR/PLON (CTM3 params.h/DUST params.h):', IPAR,PLON
       print*,'    JPAR/PLAT (CTM3 params.h/DUST params.h):', JPAR,PLAT
       print*,'    LPAR/PLEV (CTM3 params.h/DUST params.h):', LPAR,PLEV
       print*,'    Stopping!'
       stop
    end if


    !// Initialize number of DUST components
    dust_idx_tot = 0

    !// Get transport numbers for DUST tracers
    do N = 1, NPAR
       tracername = TNAME(N)
       if (tracername(1:4) .eq. 'DUST') then
          dust_idx_tot = dust_idx_tot + 1
          dust_trsp_idx(dust_idx_tot) = N
       end if
    end do


    !// Check if you read any DUST tracers at all
    if (dust_idx_tot .eq. 0) then
       write(6,'(a)') f90file//':'//subr//': '// &
            'YOU HAVE NO COMPONENTS STARTING WITH DUST'
       stop
    else if (dust_idx_tot .gt. Ndustbins) then
       write(6,'(a)') f90file//':'//subr//': DUST tracer problem:'
       write(6,'(a,i3)') '* Too many DUST tracers found ',dust_idx_tot
       write(6,'(a,i3)') '  Max DUST tracers in dust_oslo ',dust_idx_tot
       stop
    else if (dust_idx_tot .eq. Ndustbins) then
       write(6,'(a,i5)') f90file//':'//subr//': '// &
            'Number of DUST tracers/bins found ',dust_idx_tot
    else
       write(6,'(a,i5)') f90file//':'//subr//': '// &
            'DUST not set up to use dst_idx_tot less than Ndustbins ',dust_idx_tot
       stop
    end if



    !// Plus more diagnostic stuff
    dustbudgetacctime = 0._r8
    dustaverages_written = 0

    dust_production(:,:,:) = 0._r8
    dust_drydep(:,:,:)     = 0._r8
    dust_lscv(:,:,:)       = 0._r8
    dust_cscv(:,:,:)       = 0._r8
    dust_burden(:,:,:)     = 0._r8


    lastburden = 0._r8
    lastlscv = 0._r8
    lastcscv = 0._r8
    lastprod = 0._r8
    lastddep = 0._r8
    lastchem = 0._r8

    !// Initialize
    Lfirstdust = .TRUE.  !Dust has never been called before

    !// Get gaussian weights and latitudes in radians (used in dust model)
    call GAUSST2(JPAR,plat_rdn,latwts)

    !// Initialize DEAD stuff
    !// Parameters in cmn_parameters.f90 are used by inirun also
    call inirun(dt_dust,latwts,real(NDAY,r8),mbl_name,mbl_fudgefactor)

    write(6,'(a)') '* DUST initialized'

    !// --------------------------------------------------------------------
  end subroutine dust_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine dustinput_met(NOPS,MP)   !I 1/timestep [1/s]
    !// --------------------------------------------------------------------
    !// Purpose: Make special input for dust calculations
    !// based on standard input in the pwindXX.F file
    !//
    !// Author: Alf Grini
    !//
    !// To CTM3:
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    use cmn_ctm, only: AREAXY, NRMETD, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_met, only: SFT, P, PRECLS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)   :: NOPS, MP

    !// Local
    integer :: I,J,II,JJ,L
    real(r8) :: DTMET
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'dustinput_met'
    !// --------------------------------------------------------------------

    if (NOPS .ne. 1) then
       write(6,'(a,i5)') f90file//':'//subr//': '// &
            'Wrong NOPS',NOPS
       stop
    end if
    DTMET = 24._r8 / real(NRMETD,r8) * 3600._r8 !// met. time step [s]

    !// Set Surface thermal radiation downwards (Agreed with Gunnar Myhre)
    if (.not.LSTRD) STRD(:,:,MP) = 350._r8  !W/m2

    !// Possibly initialize Solar flux (recalculated in SWradbdg)
    if (.not.LSSRD) SSRD(:,:,MP) = 700._r8  !W/m2 downwards
    !SSRU(:,:,MP) = 0.e5_r8 * SSRD(:,:,MP) !W/m2 upwards


    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1
        !// Set surface pressure in Pa
        SFP(II,JJ,MP) = 100._r8 * P(I,J)

        !// Set Snow dependent on Surface temperature
        !// Assume 10 cm snow if surface T < 0 C
        if (SFT(I,J) .lt. 273._r8) then
           SNOW(II,JJ,MP) = 0.1_r8   !10 cm snow = 0.01m water equivalent
        else
           SNOW(II,JJ,MP) = 0._r8   !Now snow
        end if
      end do
    end do

    if (.not. LSWVL1) then
    !// From here and out, I calculate SWVL1
    !// Liberate first place in rainlast24h
    !// Remove data in place #8
    do L = 8,2,-1               !Go from 8 to 2 decreasing by 1
       rainlast24h(L,:,:,MP) = rainlast24h(L-1,:,:,MP)
    end do

    !// Get rain (mm) for this timestep and put in first place of rainlast24h
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        rainlast24h(1,II,JJ,MP) = & !in mm
             PRECLS(I,J,1) &   !kg/sek
             /1000._r8 &       !kg/m3
             *DTMET &         !sek (3h)
             /AREAXY(i,j) &  !m2
             *1000._r8        !mm/m
      end do
    end do

    !// Use SWVL1 as dummy for rain (in mm) last 24h
    !// Initialize it first
    SWVL1(:,:,MP) = 0._r8
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1
            
          !// Initialize rain in this grid
          rainingrid = 0._r8 
          do L=1,8
             rainingrid = rainingrid + rainlast24h(L,II,JJ,MP) !sum up in mm
          end do
            
          if (rainingrid .gt. rainlimit) then
             !// Rainlimit is in mm/h and is defined in cmndust.F

             !// Get time left to dry up last big rainfall (real)
             timeleft = real(dryuptime(II,JJ,MP) &
                             - time_since_rain(II,JJ,MP), r8)

             !// Get time needed to dry up rainfall in last 24h
             !// The assumption is that it takes 0.5 days to dry up
             !// 0.5mm-->5days to dry up 5 mm
             newtime = min(5._r8, rainingrid) & !mm/day
                       *(86400._r8)             !seconds per day
               
             !// Check if newtime > time left already
             if (newtime .gt. timeleft) then !We have to adjust the dryuptime
                dryuptime(II,JJ,MP) = nint(newtime)
                time_since_rain(II,JJ,MP) = 0
             else
                time_since_rain(II,JJ,MP) = time_since_rain(II,JJ,MP) &
                                            + nint(DTMET)
             end if
               
          else !Rain is smaller than rainlimit (less than 0.5 mm last 24h)
               
             !// Count time since a day with rainfall
             time_since_rain(II,JJ,MP) = time_since_rain(II,JJ,MP) &
                                         + nint(DTMET)
          end if
            
      end do
    end do

    !// Set the value of swvl1dust based on time_since_rain array
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        if (time_since_rain(II,JJ,MP) .le. dryuptime(II,JJ,MP)) then
           SWVL1(II,JJ,MP) = 0.5_r8  !(m3/m3) A very high soil water content
        else
           SWVL1(II,JJ,MP) = 0.01_r8 !(m3/m3) A very low soil water content
           dryuptime(II,JJ,MP) = 0       !Soil is dry, takes 0 sec to dry up
        end if
      end do
    end do
    end if !//if (.not. LSWVL1) then
    !// --------------------------------------------------------------------
  end subroutine dustinput_met
  !// ----------------------------------------------------------------------
      




  !// ----------------------------------------------------------------------
  subroutine SWradbdg_get(tau,lat,lon,soldec,CLDFRC,albedo,MP)
    !// --------------------------------------------------------------------
    !// Purpose: Given latitude, solar declination angle, local time,  
    !// latitudes, albedo and cloud cover, we should be able to calculate
    !// approximatel solar radiation in and out at the surface
    !//
    !// Theory: 
    !// An introduction to Three-Dimensional Climate modeling,
    !// W.M. Washington and C.L. Parkinson, Oxford University press, 1986
    !// pp 102-104 (eq. 3.137 to 3.139)
    !// For UiO subroutine users: Gunnar Myhre has a copy of this book
    !//
    !// Author: Alf Grini, 2002
    !//
    !// To CTM3
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input:
    integer, intent(in)   :: MP       !JI-block number
    real(r8), intent(in)  :: tau      !Model hour
    real(r8), intent(in)  :: lat(JPAR)  !Model latitudes in radians
    real(r8), intent(in)  :: lon(IPAR)  !Model longitude in radians
    real(r8), intent(in)  :: soldec   !Solar declination angle for this month
    real(r8), intent(in)  :: albedo(IPAR,JPAR) !Surface albedo (0-1)
    real(r8), intent(in)  :: CLDFRC(IPAR,JPAR,LWEPAR) !Cloud fraction


    !// Local:
    !// PHYSICAL CONSTANTS
    real(r8), parameter :: &
         S0 = 1367._r8,&       !Solar constant (W/m2)
         Ta = 0.76_r8,&        !Atmospheric transmisivity
         pi = CPI,&            !pi
         pi180 = pi / 180._r8,&
         cloudreducefactor = 0.75_r8 !Factor by which we reduce incoming
                                     !solar radiation in case of 100 %
                                     !cloud cover

    !// VARIABLES
    integer  :: I,J,II,JJ     !Counters for lon/lat
    real(r8) :: &
         H,&           !Hour angle
         cos_ssa,&     !Cosine of solar zenith angle
         phi,&         !latitude in radians
         delta,&       !Solar declination angle in radians
         maxcloud      !Maximum cloud cover in column
    !// --------------------------------------------------------------------

    !//To be consistant with the reference: 
    !// soldec is delta
    !// phi is lat(J)
    !// H is H

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// phi id lat(j), latitude in radians
          phi = lat(j)

          !// solar declination angle (soldec) is delta which is degrees
          !// in input
          delta = soldec * pi180

          !// Hour angle (0 at solar noon)
          !// Thus here at tau=0, we have H=0 for lon=0
          !// Formula fetched from pphot.f
          H = (tau*15._r8 - 180._r8)*pi180 + lon(I)
     
          !// Get max cloud fraction in colum
          maxcloud = maxval(cldfrc(I,J,:))


          !// Calculate cosinus SSA
          !// Eq. 3.137
          COS_SSA = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(H)
          !// Limit this to zero (don't want negative radiation in)
          COS_SSA = max(0.0,COS_SSA)

          !// Calculate Solar radiation downwards (W/m2)
          !// Eq. 3.139, mulitply with atmospheric transmissivity

          SSRD(II,JJ,MP) = S0 * COS_SSA * Ta

          !// Personal communication with Gunnar Myhre: 
          !// Agreed to let SSRD decrease linearly 
          !// to 1/4 of original value at 100 % cloud cover    
          SSRD(II,JJ,MP) = SSRD(II,JJ,MP) &
                           - SSRD(II,JJ,MP) * (cloudreducefactor*maxcloud)

          !// Calculate Solar radiation upwards (W/m2)
          !SSRU(II,JJ,MP) = albedo(I,J) * SSRD(II,JJ,MP)


       end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine SWradbdg_get
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine psadjust(AIRB,AREAXY,MP)
    !// --------------------------------------------------------------------
    !// Purpose: Adjust Ps, account for advection.
    !// Theory: None, really
    !// Author: Alf Grini
    !//
    !// To CTM3:
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP ! IJ-block number
    real(r8), intent(in)  :: AIRB(LPAR,IDBLK,JDBLK) ![kg] Air after advection
    real(r8), intent(in)  :: AREAXY(IPAR,JPAR)      ![m2] area of gridcells  

    !// Locals
    integer             :: i,j,II,JJ,l       !Counter variables
    !// --------------------------------------------------------------------

    !// Initialize ps
    SFP(:,:,MP) = 0._r8


    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        do L = 1, LPAR

          !// Add up air in kg
          SFP(II,JJ,MP) = SFP(II,JJ,MP) & !Old air mass in kg
                          + airB(L,II,JJ)  !Air mass in kg
        end do !// do L = 1, LPAR

        !//Change to Pa (=kgm/s2/m2)
        SFP(II,JJ,MP) = SFP(II,JJ,MP) / AREAXY(I,J) * G0

      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

  
    !// --------------------------------------------------------------------
  end subroutine psadjust
  !// ----------------------------------------------------------------------
        




  !// ----------------------------------------------------------------------
  subroutine dust_globalupdate(NDAY)
    !// --------------------------------------------------------------------
    !// Update global DUST parameters.
    !// Called from update_chemistry in main_oslo.f90.
    !//
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    use cmn_ctm, only: LYEAR
    use dsttvbds, only:dst_tvbds_ntp    !Adjust time varying dust datasets
    !// --------------------------------------------------------------------
    Implicit NONE
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)   :: NDAY
    !// --------------------------------------------------------------------

    !// Get calendar day (1-366) as real
    if (.not.LYEAR) then
       calday = mod(nday-1,365) + 1
    else
       calday = mod(nday-1,366) + 1
    end if

    !// Adjust time varying dust input
    call dst_tvbds_ntp(calday)


    !// --------------------------------------------------------------------
  end subroutine dust_globalupdate
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine oro_set(pland,sft,oro)
    !// --------------------------------------------------------------------
    !// Purpose: Set orography usable in MATCH with orography data from
    !// ECMWF model.
    !//
    !// Author: Alf Grini
    !//
    !// To CTM3:
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in)      ::sft   !Surface temp(K)
    real(r8), intent(in)      ::pland !land fraction (0-1)
    real(r8), intent(out)     ::oro   !oreography (0,1,2)
    !// oro=0 for ocean, 1 for land and 2 for sea-ice

    real(r8), parameter   :: templim_ice=277._r8  !Limit for ice (+4C) 
    real(r8), parameter   :: fraclim_sea=0.5_r8  !Limit for sea
    !// --------------------------------------------------------------------

    if(sft .lt. templim_ice .and. pland .lt. fraclim_sea) then
       !// ----------Check for sea ice-------------------------
       oro = 2._r8       !Sea ice
    else if(sft .gt. templim_ice .and. pland .lt. fraclim_sea) then
       !/// ----------Check for sea--------------------------
       oro = 0._r8
    else
       !// ----------LAND--------------------------------
       oro = 1._r8
    end if

    !// --------------------------------------------------------------------
  end subroutine oro_set
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine plevs0(ETAA, ETAB, pdel, pmid, pi, ps, sh, t, pt, tv, zi, zm)
    !// --------------------------------------------------------------------
    !// Get pressure, temperature and height properties for DUST.
    !// Author: Alf Grini
    !//
    !// To CTM3:
    !// Ole Amund Sovde, October 2009
    !// --------------------------------------------------------------------
    use cmn_parameters, only: AVOGNR, R_AIR, M_AIR, CP_AIR
    use pmgrid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    real(r8), intent(in)  :: sh(plev)    !Specific humidity (kg/kg)
    real(r8), intent(in)  :: t(plev)     !Temperature (K)
    real(r8), intent(in)  :: ps          !Surface pressure (Pa)
    real(r8), intent(in)  :: etaa(plevp) !Hybrid coordinate (hPa) 
    real(r8), intent(in)  :: etab(plevp) !hybrid coordinate (hPa)

    !// output
    real(r8), intent(out) :: pdel(plev)  !Pressure difference (Pa)
    real(r8), intent(out) :: pi(plevp)   !Interface pressure (Pa)
    real(r8), intent(out) :: pmid(plev)  !Midpoint pressure (Pa)
    real(r8), intent(out) :: pt(plev)    !Potential temperature (K)
    real(r8), intent(out) :: tv(plev)    !viritual temp (K)
    real(r8), intent(out) :: zm(plev)    !Midpoint height (m)
    real(r8), intent(inout) :: zi(plevp)   !Interface height (m)

    !// Local variables
    !real(r8), parameter :: CP_AIR=1004.64 !Heat capacity of air (UNITS ??)
    real(r8), parameter :: Pa2hpa=1.e-2_r8 !Number of hpa in one Pa
                                           !Used to convert ground pressure
    real(r8), parameter :: hpa2pa=1.e2_r8  !Number of Pa in one hPa
                                           !Used to convert back pressures
    real(r8) :: cappa                     !R/CP
    real(r8) :: rhlp1,rhlp2,rhlp3         !Help variables, see below
    real(r8) :: zg0                       !Inverse of g0
    real(r8) :: DZ,DZ_MID                 !Help variable for height of layer
    real(r8) :: p1,p2,pm                  !Help variables for pressures
    real(r8) :: p1a,p2a,p1b,p2b           !Help variables for hybrid coord
    real(r8) :: pmid_a,pmid_b             !Help variables for hybrid coord
    integer  :: lev_idx,lon_idx           !Counting variables
    !// --------------------------------------------------------------------

    !// Define help constants
    cappa = R_AIR/CP_AIR
    RHLP1 = 3._r8/5._r8               !For use in viritual temp
    RHLP2 = AVOGNR * 1.e-3_r8 / M_AIR !molecules in 1m3 air (?)
    RHLP3 = M_AIR/18._r8              !Molwt air/molwt water
    ZG0   = 1._r8/G0                  !Inverse of gravity

    !// initialize values at bottom of column

    !// Height of bottom of first layer above surface: 0 by definition
    zi(plevp) = 0._r8

    !// Pressure at surface equals surface pressure
    pi(plevp) = ps

    do lev_idx = plev,1,-1
       !Note, not plevp because we calculate for l+1 too

       !// Define hybrid layer coefficients
       !// Rember that it is turned around because Match and CTM use different
       !// definitions of what is L(top) and L(bottom)
       !// Here we use MATCH definitions and then L+1 is below L
       !// And remember that etaa and etab is in hPa mode
       P1A = ETAA(lev_idx+1)   !First value of hybrid coord (bottom of layer)
       P2A = ETAA(lev_idx)     !First value of hybrid coord (top of layer)
       P1B = ETAB(lev_idx+1)   !Second value of hybrid coord (bottom of layer)
       P2B = ETAB(lev_idx)     !Second value of hybrid coord (top of layer)

       !// Midpoint values
       PMID_A = 0.5_r8 * (P1A + P2A) !First value of hybrid coord (midpoint)
       PMID_B = 0.5_r8 * (P1B + P2B) !Second value of hybrid coord (midpoint)

       !// Pressure at bottom of layer (hPa)
       P1 = P1A + P1B*(Pa2hpa*PS)
       !// Pressure at bottom interface (hPa)
       pi(lev_idx+1) = P1
       !// Pressure at top of layer (hPa)
       P2 = P2A + P2B*(Pa2hpa*PS)
       !// Pressure at top interface (hPa)
       pi(lev_idx) = P2
       !// Pressure at layer midpoint (hPa)
       PM = PMID_A + PMID_B*(pa2hpa*PS)
       !// Pressure at midpoint (hPa)
       pmid(lev_idx) = PM

       !// Pressure difference;  Layer below-layer above gives positive
       !// difference in hPa
       pdel(lev_idx) = pi(lev_idx+1) - pi(lev_idx)

       !// Viritual temp
       tv(lev_idx) = T(lev_idx)*(1._r8 + RHLP1*sh(lev_idx))

       !// Potential temp (potential temp. refer to 1000hPa, not PS)
       pt(lev_idx) = T(lev_idx)*(1000._r8/PM)**cappa

       !// Height of every grid box
       !// Remember: P1>P2 because P1 is below P2 ==> DZ >0
       DZ = tv(lev_idx)*R_AIR*ZG0*Log(P1/P2)

       !// Height of interface (Remember that L+1 is below L)
       zi(lev_idx) = zi(lev_idx+1) + DZ
       !// Height up to midpoint of every grid box
       DZ_MID = tv(lev_idx)*R_AIR*ZG0*log(P1/PM)
       !// Height of midpoint
       zm(lev_idx) = zi(lev_idx+1) + DZ_MID

    end do
    !// Transfer all pressures to Pa
    pmid(:) = pmid(:)*hPa2Pa
    pi(:)   = pi(:)*hpa2Pa
    pdel(:) = pdel(:)*hpa2pa

    !// ------------------------------------------------------------------
  end subroutine plevs0
  !// ------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine plevs00(ETAA, ETAB, pdel, pmid, pi, ps, sh, t, pt, tv)
    !// --------------------------------------------------------------------
    !// Calculate pressure, virtual temperature and potential temperature.
    !//
    !// Ole Amund Sovde, January 2016
    !// --------------------------------------------------------------------
    use cmn_parameters, only: AVOGNR, R_AIR, M_AIR, CP_AIR
    use pmgrid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    real(r8), intent(in)  :: sh(plev)    !Specific humidity (kg/kg)
    real(r8), intent(in)  :: t(plev)     !Temperature (K)
    real(r8), intent(in)  :: ps          !Surface pressure (Pa)
    real(r8), intent(in)  :: etaa(plevp) !Hybrid coordinate (hPa) 
    real(r8), intent(in)  :: etab(plevp) !hybrid coordinate (hPa)

    !// output
    real(r8), intent(out) :: pdel(plev)  !Pressure difference (Pa)
    real(r8), intent(out) :: pi(plevp)   !Interface pressure (Pa)
    real(r8), intent(out) :: pmid(plev)  !Midpoint pressure (Pa)
    real(r8), intent(out) :: pt(plev)    !Potential temperature (K)
    real(r8), intent(out) :: tv(plev)    !viritual temp (K)

    !// Local variables
    real(r8), parameter :: Pa2hpa=1.e-2_r8 !Number of hpa in one Pa
                                           !Used to convert ground pressure
    real(r8), parameter :: hpa2pa=1.e2_r8  !Number of Pa in one hPa
                                           !Used to convert back pressures
    real(r8) :: cappa                     !R/CP
    real(r8) :: rhlp1                     !Help variables, see below
    real(r8) :: p1,p2,pm                  !Help variables for pressures
    real(r8) :: p1a,p2a,p1b,p2b           !Help variables for hybrid coord
    real(r8) :: pmid_a,pmid_b             !Help variables for hybrid coord
    integer  :: lev_idx                   !Counting variables
    !// --------------------------------------------------------------------

    !// Define help constants
    cappa = R_AIR/CP_AIR
    RHLP1 = 3._r8/5._r8               !For use in viritual temp


    !// Pressure at surface equals surface pressure
    pi(plevp) = ps

    do lev_idx = plev,1,-1
       !Note, not plevp because we calculate for l+1 too

       !// Define hybrid layer coefficients
       !// Rember that it is turned around because Match and CTM use different
       !// definitions of what is L(top) and L(bottom)
       !// Here we use MATCH definitions and then L+1 is below L
       !// And remember that etaa and etab is in hPa mode
       P1A = ETAA(lev_idx+1)   !First value of hybrid coord (bottom of layer)
       P2A = ETAA(lev_idx)     !First value of hybrid coord (top of layer)
       P1B = ETAB(lev_idx+1)   !Second value of hybrid coord (bottom of layer)
       P2B = ETAB(lev_idx)     !Second value of hybrid coord (top of layer)

       !// Midpoint values
       PMID_A = 0.5_r8 * (P1A + P2A) !First value of hybrid coord (midpoint)
       PMID_B = 0.5_r8 * (P1B + P2B) !Second value of hybrid coord (midpoint)

       !// Pressure at bottom of layer (hPa)
       P1 = P1A + P1B*(Pa2hpa*PS)
       !// Pressure at bottom interface (hPa)
       pi(lev_idx+1) = P1
       !// Pressure at top of layer (hPa)
       P2 = P2A + P2B*(Pa2hpa*PS)
       !// Pressure at top interface (hPa)
       pi(lev_idx) = P2
       !// Pressure at layer midpoint (hPa)
       PM = PMID_A + PMID_B*(pa2hpa*PS)
       !// Pressure at midpoint (hPa)
       pmid(lev_idx) = PM

       !// Pressure difference;  Layer below-layer above gives positive
       !// difference in hPa
       pdel(lev_idx) = pi(lev_idx+1) - pi(lev_idx)

       !// Viritual temp
       tv(lev_idx) = T(lev_idx)*(1._r8 + RHLP1*sh(lev_idx))

       !// Potential temp (potential temp. refer to 1000hPa, not PS)
       pt(lev_idx) = T(lev_idx)*(1000._r8/PM)**cappa


    end do
    !// Transfer all pressures to Pa
    pmid(:) = pmid(:)*hPa2Pa
    pi(:)   = pi(:)*hpa2Pa
    pdel(:) = pdel(:)*hpa2pa

    !// ------------------------------------------------------------------
  end subroutine plevs00
  !// ------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine dust_master(BTT, AIRB, BTEM, DTCHM, NOPS, MP)
    !// --------------------------------------------------------------------
    !// Master interaction between Oslo CTM3 and DEAD.
    !// Called from master_oslo in main_oslo.f90
    !//
    !// While DEAD takes mass mixing ratio of dust into its calculations
    !// and converts to mass by calculating air mass and density, we
    !// send in the masses of dust and air, along with calculated volume.
    !// This is to make the flux diagnostics fully consistent with
    !// the actual change in tracers.
    !//
    !// Ole Amund Sovde, February 2015, October 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR
    use cmn_ctm, only: ETAA, ETAB, GMTAU, XGRD, YGRD, &
         MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, AREAXY, SOLDEC, PLAND
    use cmn_met, only: CLDFR, Q, SA, SFT, SD, UMS, VMS, ZOFLE
    use cmn_oslo, only: DV_IJ
    use dstsfc, only: flx_rad_sfc_set, tpt_gnd_soi_set, vwc_sfc_set
    use dsttvbds, only:dst_tvbds_ntp    !Adjust time varying dust datasets
    use pmgrid, only: plev,plevp       !Physical grid
    !// --------------------------------------------------------------------
    Implicit NONE
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in)   :: NOPS, MP
    real(r8), intent(in)  :: DTCHM
    real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: AIRB, BTEM
    !// Input / Output
    real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT

    !// Locals
    real(r8) :: DUSTT(LPAR,Ndustbins) ![kg] dust mass
    real(r8) :: AIRM(LPAR) ![kg] air mass
    real(r8) :: VOLU(LPAR) ![m3] air volume
    real(r8) :: AREA       ![m2] grid box area
    real(r8) :: &
         oro, &               ![0=ocean,1=land,2=sea ice]
         ts, &                !Surface temperture
         snowh, &             !Snow height
         fsds, &              !Solar downwelling radiation (surf) 
         ps, &                !Surface pressure
         fsus, &              !Solar upwelling radiation (surf)
         flds, &              !Surface termal (longwave) radiation downwards
                              ! (longwave downwelling surface)
         vwcsfc, &            !Volumetric water content
         latrad, lonrad, &    !lat/lon in radians
         bcol(Ndustbins)      !total colum burdens
    integer  :: I,J,II,JJ,L,N   !Counting variables

    !// [kg] Dry deposited
    real(r8), dimension(Ndustbins) :: drydep
    !// [kg] Produced
    real(r8), dimension(Ndustbins) :: production

    !// Arrays with changed dimensions
    real(r8), dimension(plevp) :: &
         ETAA_MATCH, & !hybrid coord (match)
         ETAB_MATCH, & !Hybrid coord (match)
         zi            ![m] height of interstices
    
    real(r8), dimension(plev)  :: &
         Q_match, &    ![kg/kg] specific humidity with dimensions switched
         T_match, &    ![K] Temperature with dimensions switched
         zm            ![m] Height of grid box center

    real(r8) :: &
         UMPS_sfc, & ![m/s] sfc wind in midpoint (zonal)
         VMPS_sfc    ![m/s] sfc wind in midpoint (meridional)
    !// --------------------------------------------------------------------

    !// Initialize meteorological variables for DUST
    if (NOPS .eq. 1) call dustinput_met(NOPS,MP)

    !// Get right surface pressure after advection
    !// Only every CHMCYCLE?
    call psadjust(AIRB,AREAXY,MP)


    !// Get radiation budget for short wave radiation, if it is not
    !// set from meteorology data.
    if (.not.LSSRD) &
         call swradbdg_get(gmtau,ygrd,xgrd,soldec,CLDFR,SA,MP)



    !// Loop through IJ-block
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1


        !// Reversing vertical order for all variables levels because 1 is
        !// first level in CTM and plevp is first level in MATCH.
        !// (plevp = plev+1, interfaces)
        do L = 1, plevp
          etaa_match(plevp-L+1) = etaa(L) !Turn around Etaa levels
          etab_match(plevp-L+1) = etab(L) !Turn around Etab levels
          zi(plevp-L+1) = ZOFLE(L,I,J) - ZOFLE(1,I,J)
        end do

        do L = 1, plev
          Q_match(L)    = Q(I,J,plev-L+1) !Turn around Spec humi
          T_match(L)    = BTEM(plev-L+1,II,JJ) !Turn around temp
          zm(L) = 0.5_r8 * (zi(L) + zi(L+1)) !Height of box midpoint
        end do

        UMPS_sfc = UMS(1,I,J) !sfc zonal wind
        VMPS_sfc = VMS(1,I,J) !sfc meridional wind

        !// We have to set oreography consistent with what is used in Match
        call oro_set(PLAND(i,j), SFT(i,j), ORO)


        !// DUST components from BTT (convert to mixing ratio)
        do n = 1, Ndustbins
           do L = 1, plev
              !// Retrieve dust masses [kg]
              DUSTT(plev-L+1,N) = BTT(L,dust_trsp_idx(N),II,JJ)
           end do !// do L = 1, plev
        end do !// do n = 1, Ndustbins

        do L = 1, plev
           !// Retrieve air mass and grid box volume
           AIRM(plev-L+1) = AIRB(L,II,JJ)
           VOLU(plev-L+1) = DV_IJ(L,II,JJ,MP)
        end do !// do L = 1, plev


        !// Pick out 2D/1D properties
        ps    = SFP(II,JJ,MP)      !// Setting surface pressure
        ts    = SFT(I,J)           !// Surface temperature
        fsds  = SSRD(II,JJ,MP)     !// Solar downwelling radiation (surf) 
        !// Solar upwelling from surface is calculated from albedo and downward
        !// radiation.
        !fsus  = SSRU(II,JJ,MP)     !// Solar upwelling radiation (surf)
        fsus  = SA(I,J) * SSRD(II,JJ,MP) !// Solar upwelling radiation (surf)
        flds  = STRD(II,JJ,MP)     !// Surface termal (longwave) rad. downwards
        vwcsfc= SWVL1(II,JJ,MP)!// Volumetric water content
        snowh = SNOW(II,JJ,MP)     !// Snow height
        !snowh = SD(I,J)     !// Snow height
        latrad= YGRD(J)            !// Latitude (rad)
        lonrad= XGRD(I)            !// Longitude (rad)
        area  = AREAXY(I,J)        !// Grid box area [m2]

        !// Should store DUST variables here, only each NOPS
        !// In principle, DUST should not need duplicate variables, but
        !// it was kept in CTM2 to make the code easily updateable.
        if (NOPS .eq. 1) then
           !// Set fields used by mobilization routines (dstsfc.F90)
           call flx_rad_sfc_set(J, I, &
                flds, &      ! I [W m-2] Longwave downwelling flux at surface
                fsds - fsus) ! I [W m-2] Solar flux absorbed by ground
           call tpt_gnd_soi_set(J, I, &
                ts, &        ! I [K] Ground temperature
                ts)          ! I [K] Soil temperature = ground temp.
           call vwc_sfc_set(J, I, &
                vwcsfc)      ! I [m3 m-3] Volumetric water content
        end if

        !// Do master column
        call dust_column(DUSTT, AIRM, VOLU, AREA, &
             I,J,II,JJ,MP, DTCHM, ETAA_MATCH, ETAB_MATCH, &
             Q_MATCH, T_MATCH, UMPS_sfc, VMPS_sfc, &
             zi, zm, &
             PS, ORO, ts, fsds, fsus, flds, vwcsfc, snowh, latrad, lonrad, &
             drydep,production)

        !// Drydep budget
        do n = 1, Ndustbins
           !// Save dry deposition [kg] in global diagnose
           call add2drydepdiag(II,JJ,MP,dust_trsp_idx(n),drydep(n))
           !// Drydep for all DUST tracers [kg]
           dust_drydep(N,I,J) = dust_drydep(N,I,J) + drydep(N)
        end do


        !// Production budget
        do n = 1, Ndustbins
           !// Production for all DUST tracers [kg]
           dust_production(N,I,J) = dust_production(N,I,J) + production(N)
        end do
 

        !// DUSTT back to BTT and re-reversing vertical order
        bcol(:) = 0._r8          !// Burden [kg] in this column
        do n = 1, Ndustbins
           do L = 1, plev  
              !// Put back new dust masses
              BTT(L,dust_trsp_idx(N),II,JJ) = DUSTT(plev-L+1,N)
              !// Find column burden
              bcol(N) = bcol(N) + DUSTT(plev-L+1,N)
           end do !// do L = 1, plev
        end do !// do n = 1, Ndustbins

        !// Save burden [kg*s]
        !// Burden is weighted with DTCHM to allow for a possible
        !// varying DTCHM (between NOPS, not varying between MP-blocks).
        do n = 1, Ndustbins
           dust_burden(N,I,J) = dust_burden(N,I,J) + bcol(N) * DTCHM
        end do

      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// Save total time (seconds) for average
    !// Only do this for one IJ-block, to count the time only once.
    if (MP .eq. 1) dustbudgetacctime = dustbudgetacctime + DTCHM

    !// --------------------------------------------------------------------
  end subroutine dust_master
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine dust_column(DUSTT,AIRM,VOLU,AREA, &
       I,J,II,JJ,MP, DTDST, ETAA_m, ETAB_m, &
       SH1, T1, UMS, VMS, &
       zi, zm, &
       PS, ORO, TS, FSDS, FSUS, FLDS, VWCSFC, SNOWH, latrad, lonrad, &
       drydep, production)
    !// --------------------------------------------------------------------
    !// Column routine to calculate DUST as a column. Based on physlic, which
    !// covers physics parameterizations and output fields.
    !// Previous purpose was to call one latitude slice of dust model (DEAD),
    !// but for CTM3 it calls one column only.
    !// Original file physlic.F was the physics driver for the MATCH model.
    !//
    !// Also, the routines now work on the dust masses instead of
    !// mass mixing ratio. This was corrected because the actual change
    !// in tracer masses were not fully consistent with the flux
    !// calculations. This correction is small.
    !//
    !// Ole Amund Sovde, February 2015, October 2009
    !// --------------------------------------------------------------------
    use pmgrid, only: plev,plevp      !Physical grid
    use dstdpsdry,only:dst_dps_dry ! [mdl] Dry deposition driver
    use dstsfc, only: flx_rad_sfc_set, tpt_gnd_soi_set, vwc_sfc_set
    use dstmbl, only: dst_mbl
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in) :: I,J,II,JJ,MP
    real(r8), intent(in), dimension(plevp) :: &
            ETAA_m, & !// hybrid coord (match)
            ETAB_m, & !// Hybrid coord (match)
            zi        !// [m] Height of intertices

    real(r8), intent(in), dimension(plev) :: &
            SH1, &  !// [kg/kg] specific humidity with dimensions switched
            T1, &   !// [K] Temperature with dimensions switched
            AIRM, & !// [kg] Air mass
            VOLU, & !// [m3] Air (grid box) volume
            ZM      !// [m] Height of box midpoint

    real(r8), intent(in) :: &
            UMS, &  !// [m/s] Surface wind in midpoint (zonal)
            VMS     !// [m/s] Surface wind in midpoint (meridional)

    real(r8), intent(in) :: &
            DTDST, &          !// Time step
            ORO, &            !// Surface "orography" (0.d0, 1.d0 or 2.d0)
            ts, &             !// Surface temperture
            snowh, &          !// Snow height
            ps, &             !// Surface pressure
            fsds, &           !// Solar downwelling radiation (surf) 
            fsus, &           !// Solar upwelling radiation (surf)
            flds, &           !// Surface termal (longwave) radiation downwards
            vwcsfc, &         !// Volumetric water content
            latrad, lonrad, & !// Lat/lon in radians
            AREA              !// [m2] Grid box area

    !// Input/output: DUST species
    real(r8), intent(inout) :: DUSTT(LPAR,Ndustbins)
    !// Output 
    real(r8), intent(out) :: &
         drydep(Ndustbins), &   !// [kg] Dry deposited
         production(Ndustbins)  !// [kg] Produced


    !// Locals
    integer :: L, N, k, &
         ioro, lchnk, ncol

    real(r8), dimension(plev) :: &
            shadv, &  !// internal specific humidity after advection on input.
            pdel, &   !// pressure difference across layers
            pmid, &   !// pressure at layer midpoints
            pt0, &    !// potential temperature
            tv        !// virtual temperature

    real(r8), dimension(plevp) :: &
            pint   !// Pressure at interfaces (Pa)

    real(r8), dimension(1) :: obuf
    real(r8) :: mno_lng, &      !// [m] Monin-Obukhov length
              wnd_frc, &      !// [m s-1] Surface friction velocity
              wnd_rfr, &      !// [m s-1] Wind speed at reference height
              fsns            !// [W m-2] Solar flux absorbed by ground
    !// --------------------------------------------------------------------

    !// Set absolute humidity after advection equal to humidity before
    !// advection
    shadv(:) = sh1(:)

    !// Initialize some variables; pressure and heights in interfaces
    !// and midpoints
    ioro = nint( oro )

    !// Calculate pressure and viritual & potential temp
    call plevs00(etaa_m,etab_m,pdel,pmid,pint,ps,sh1,t1,pt0,tv)

    !// Dry deposition processes for grid (I,J) (checked ok for CTM3)
    !// At surface, which is plev in DEAD.
    call dst_dps_dry(lchnk, ncol, obuf, &
          zm(plev), &            ! I [m] Midpoint height above surface
          J, &                  ! I [idx] Model latitude index
          I, &                  ! I [idx] Model longitude index
          mno_lng, &            ! O [m] Monin-Obukhov length
          oro, &                ! I [frc] Orography
          pdel, &               ! I [Pa] Pressure thickness
          pmid, &               ! I [Pa] Midlayer pressure
          shadv, &              ! I [kg kg-1] Water vapor mixing ratio
          DUSTT, &              ! I/O [kg] Dust mass
          AIRM, &               ! I [kg] Air mass
          VOLU, &               ! I [m3] Air volume
          AREA, &               ! I [m2] Grid box area
          snowh, &              ! I [m] Equivalent liquid water snow depth
          dtdst, &              ! I [s] Adjustment timestep
          t1, &                 ! I [K] Temperature
          pt0(plev), &          ! I [K] Potential temperature at surface
          ts, &                 ! I [K] Surface temperature
          wnd_frc, &            ! O [m s-1] Surface friction velocity
          vms, &                ! I [m s-1] Meridional wind component
          wnd_rfr, &            ! O [m s-1] Wind speed at reference height
          ums, &                ! I [m s-1] Zonal wind component 
          drydep)               ! O [kg] Dry deposited



!    !// Set fields used by mobilization routines
!    fsns = fsds - fsus ! I [W m-2] Solar flux absorbed by ground

!    !// Set surface properties (dstsfc.F90) (checked ok for CTM3)
!    call flx_rad_sfc_set( &
!          J, &                ! I [idx] Latitude index
!          I, &                ! I [idx] Longitude index
!          flds, &             ! I [W m-2] Longwave downwelling flux at surface
!          fsns)               ! I [W m-2] Solar flux absorbed by ground


!    !// (in dstsfc.F90) (checked ok for CTM3)
!    call tpt_gnd_soi_set( &
!         J, &                ! I [idx] Latitude index
!         I, &                ! I [idx] Longitude index
!         ts, &               ! I [K] Ground temperature
!         !//fxm: Assume soil temperature is same as surface temperature
!         ts)                 ! I [K] Soil temperature


!    !// (in dstsfc.F90) (checked ok for CTM3)
!    call vwc_sfc_set( &
!         J, &                ! I [idx] Latitude index
!         I, &                ! I [idx] Longitude index
!         vwcsfc)             ! I [m3 m-3] Volumetric water content


    !// Driver for dust mobilisation (i.e. production)
    call dst_mbl(lchnk,ncol,obuf, &
         calday, &             ! I [day] Day of year [1.0..366.0)
         zm(plev), &         ! I [m] Midpoint height above surface
         J, &                ! I [idx] Model latitude index
         I, &                ! I [idx] Model longitude index
         latrad, &           ! I [rdn] Latitude
         oro, &              ! I [frc] Orography
         pdel(plev), &       ! I [Pa] Pressure thickness surface
         pmid(plev), &       ! I [Pa] Pressure mid-layer surface
         pint(plevp), &      ! I [Pa] Surface pressure NB: plevp
         shadv(plev), &      ! I [kg kg-1] Water vapor mixing ratio
         DUSTT(plev,:), &    ! I/O [kg] Dust mass @ sfc
         AIRM(plev), &       ! I [kg] Air mass @ sfc
         VOLU(plev), &       ! I [m3] Air volume @ sfc
         AREA, &             ! I [m2] Grid box area
         snowh, &            ! I [m] Equivalent liquid water snow depth
         dtdst, &            ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
         t1(plev), &         ! I [K] Temperature
         pt0(plev), &        ! I [K] Potential temperature
         vms, &              ! I [m s-1] Meridional wind component
         ums, &              ! I [m s-1] Zonal wind component
         production)         ! O [kg] Dust produced


    !// --------------------------------------------------------------------
  end subroutine dust_column
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine dustInstBdg(NDAY,NDAY0,NMET,NOPS,DTOPS)
    !// --------------------------------------------------------------------
    !// Print semi-instant and average global diagnostics to screen.
    !// Semi-instant means since last time dustInstBdg was called,
    !// currently each NOPS.
    !//
    !// Also calculates average lifetimes during the period defined
    !// by the diagnostic calendar. This routine is mainly for
    !// debugging; its value is rather limited as standard output.
    !//
    !// To include printouts to screen:
    !//   - set LPrntBud for budgets and/or LPrntLif for lifetimes.
    !//
    !// Routine is called at the end of the NOPS-loop in pmain,
    !// but will only do calculations if LPrntBud or LPrntLif are
    !// set.
    !//
    !// Amund Sovde Haslerud, June 2016
    !// --------------------------------------------------------------------
    use cmn_ctm, only: AREAXY, STT, NROPSM, &
         MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_diag, only: STTTN0, NTND_LSSCAV, NTND_CNSCAV, NTND_CHEM
    use cmn_oslo, only: CONVWASHOUT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NDAY0,NMET,NOPS
    real(r8), intent(in) :: DTOPS
    !// Locals
    real(r8) :: burden,lscv,cscv,prod,ddep, chem, davgtnd, totL, totP, &
         totLS,totCN,totD,burden2
    integer :: I,J,K,L,MP,II,JJ,N
    !// Parameters
    logical, parameter :: LPrntBud = .false.
    logical, parameter :: LPrntLif = .false.
    !// --------------------------------------------------------------------

    !// Return if no printouts are asked for
    if (.not. (LPrntBud .or. LPrntLif)) return

    !// Time in seconds since diagnostic calendar started
    davgtnd = dtops * real(nops, r8) &
         + dtops * real(nmet-1, r8) * real(NROPSM, r8) &
         + 86400._r8 * real(NDAY - NDAY0, r8)

    !// Instant burden of DUST
    burden = 0._r8
    do K = 1, Ndustbins
       N = dust_trsp_idx(K)
       burden = burden + sum(STT(:,:,:,N))
    end do

    !// Accumulated wet deposition (large scale)
    lscv = 0._r8
    do K = 1, Ndustbins
       N = dust_trsp_idx(K)
       lscv = lscv + STTTN0(N,NTND_LSSCAV)
    end do

    !// Accumulated change in "chemistry" loop
    chem = 0._r8
    do K = 1, Ndustbins
       N = dust_trsp_idx(K)
       chem = chem + STTTN0(N,NTND_CHEM)
    end do

    !// Accumulated wet deposition (convective)
    cscv = 0._r8
    do K = 1, Ndustbins
       N = dust_trsp_idx(K)
       cscv = cscv + STTTN0(N,NTND_CNSCAV)
    end do

    !// Accumulated production [kg]
    prod = sum(dust_production)
    !// Accumulated dry deposition [kg]
    ddep = -sum(dust_drydep)

    !// Instant totals during this nops
    if (LPrntBud) &
         write(6,'(a,5es10.2)') 'DUST INST P/D/LS/CN/dB [kg]  ', &
            prod - lastprod,& 
            ddep - lastddep,&
            lscv - lastlscv,&
            cscv - lastcscv,&
            burden - lastburden


    !// I have tested that P-L is the same as chemical tendency.
    !// It was not when DEAD calculated mass mixing ratios, but after I corrected
    !// it to calculate dust mass directly, these are the same.



    !// Lifetime [days] based on instant burden prod/loss rates.
    totD = (ddep - lastddep) / dtops
    totLS = (lscv - lastlscv) / dtops
    totCN = (cscv - lastcscv) / dtops
    totL = totD + totLS + totCN
    totP = (prod - lastprod) / dtops

    if (LPrntLif) &
         write(6,'(a,5es10.2)') 'DUST LIFE P/D/LS/CN/TL [days]', &
            burden/totP / 86400._r8, &
            -burden/totD / 86400._r8, &
            -burden/totLS / 86400._r8, &
            -burden/totCN / 86400._r8, &
            -burden/totL / 86400._r8

    !// Save values for next time
    lastburden = burden
    lastlscv = lscv
    lastcscv = cscv
    lastprod = prod
    lastddep = ddep
    lastchem = chem

    !// --------------------------------------------------------------------
  end subroutine dustInstBdg
  !// ----------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine add2drydepdiag(ii,jj,mp,n,ddmass)
    !// ------------------------------------------------------------------
    !// Saves DUST tracer mass lost to drydep.
    !//
    !// Ole Amund Sovde, September 2015
    !// ------------------------------------------------------------------
    use cmn_ctm, only: all_mp_indices
    use cmn_oslo, only: SCAV_MAP_DRY, SCAV_DD
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    integer, intent(in) :: ii,jj,mp,n
    real(r8), intent(in) :: ddmass
    !// ------------------------------------------------------------------
    SCAV_MAP_DRY(N,II,JJ,MP) = SCAV_MAP_DRY(N,II,JJ,MP) + ddmass
    !// Save totals
    SCAV_DD(N,MP) = SCAV_DD(N,MP) + ddmass
    !// ------------------------------------------------------------------
  end subroutine add2drydepdiag
  !// ------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine dust_set_ssrd(r8data)
    !//---------------------------------------------------------------------
    !// Puts metdata of size (ipar,jpar) into MP-block structure SSRD.
    !//
    !// Ole Amund Sovde, January 2016
    !//---------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: r8data(ipar,jpar)
    !// Locals
    integer :: I, J, II, JJ, MP
    !//---------------------------------------------------------------------
    !// Flag that SSRD is set from metdata
    LSSRD = .true.
    do MP = 1, MPBLK
      !// Loop over latitude (J is global, JJ is block)
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        !// Loop over longitude (I is global, II is block)
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// Change structure
          SSRD(II,JJ,MP) = r8data(I,J)
        end do
      end do
    end do
    !//---------------------------------------------------------------------
  end subroutine dust_set_ssrd
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine dust_set_strd(r8data)
    !//---------------------------------------------------------------------
    !// Puts metdata of size (ipar,jpar) into MP-block structure STRD.
    !//
    !// Ole Amund Sovde, January 2016
    !//---------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: r8data(ipar,jpar)
    !// Locals
    integer :: I, J, II, JJ, MP
    !//---------------------------------------------------------------------
    !// Flag that STRD is set from metdata
    LSTRD = .true.
    do MP = 1, MPBLK
      !// Loop over latitude (J is global, JJ is block)
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        !// Loop over longitude (I is global, II is block)
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// Change structure
          STRD(II,JJ,MP) = r8data(I,J)
        end do
      end do
    end do
    !//---------------------------------------------------------------------
  end subroutine dust_set_strd
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine dust_set_SWVL1(r8data)
    !//---------------------------------------------------------------------
    !// Puts metdata of size (ipar,jpar) into MP-block structure SWVL1.
    !//
    !// Ole Amund Sovde, January 2016
    !//---------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: r8data(ipar,jpar)
    !// Locals
    integer :: I, J, II, JJ, MP
    !//---------------------------------------------------------------------
    !// Flag that SWVL1 is set from metdata
    LSWVL1 = .true.
    do MP = 1, MPBLK
      !// Loop over latitude (J is global, JJ is block)
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        !// Loop over longitude (I is global, II is block)
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// Change structure
          SWVL1(II,JJ,MP) = r8data(I,J)
        end do
      end do
    end do
    !//---------------------------------------------------------------------
  end subroutine dust_set_SWVL1
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine dust_set_mbl_name_ff(name,fudgefactor)
    !//---------------------------------------------------------------------
    !// Sets name of mobilisation map to read from dust input.
    !// Sets fudge factor matching mobilisation map.
    !//
    !// Ole Amund Sovde, March 2016
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: fudgefactor
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'dust_set_mbl_name'
    !// --------------------------------------------------------------------
    write(6,'(a)') f90file//':'//subr// &
         ': Setting dust mobilisation variable name: '//trim(name)
    mbl_name = name
    write(6,'(a,es20.6)') f90file//':'//subr// &
         ': Setting dust mobilisation fudge factor: ',fudgefactor
    mbl_fudgefactor = fudgefactor
    !//---------------------------------------------------------------------
  end subroutine dust_set_mbl_name_ff
  !//-----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine dustbdg2file(NDAY, NDAYI, NDAY0)
    !// --------------------------------------------------------------------
    !// Writes out DUST budget averages, more comprehensive than
    !// the CTM2 dustbdg2d: Can also write budgets for all dust tracers.
    !// Depends on core diagnostics, as it fetches the diagnostics of wet
    !// scavenging from the STTTND array.
    !// Also writes totals and lifetimes to screen.
    !//
    !// Ole Amund Sovde, March 2016
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR, NTDPAR, MODEL
    use cmn_ctm, only: XDGRD, YDGRD, IYEAR, AREAXY, &
         MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_diag, only: STTTND, NTND_LSSCAV
    use cmn_oslo, only: CONVWASHOUT
    use netcdf
    use ncutils, only:  handle_error, handle_err
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NDAYI, NDAY0
    
    !// Locals
    integer :: I, J, II, JJ, MP, L, N, K     !// Indices
    real(r8) :: davgdays, davgtnd  !// Time durations
    real(r8) :: burden,lscv,cscv,prod,ddep
    real(r8) :: ZAREA(IPAR,JPAR)
    integer :: strLen
    !// Total budgets (summed up over DUST tracers)
    real(r8), dimension(IPAR,JPAR) :: &
         dstTotBdgProd, &
         dstTotBdgDdep, &
         dstTotBdgLscv, &
         dstTotBdgCscv, & 
         dstTotBdgBrdn

    !// --------------------------------------------------------------------
    !// netCDF: Dimensions ID and counters
    integer :: &
         lat_dim_id, &             !// Dimension ID for latitude
         lon_dim_id, &             !// Dimension ID for longitude
         nbr_dim_id, &             !// Dimension ID for number of tracers
         time_dim_id, &            !// Dimension ID for time
         timespan_dim_id, &        !// Dimension ID for timespan
         str_len_id, &             !// Dimension ID for string length
         time_id, &                !// Variable ID for time
         timespan_id, &            !// Variable ID for average time span (day)
         lon_id, &                 !// Variable ID for longitude
         lat_id, &                 !// Variable ID for latitude
         prod_id, ddep_id, lscv_id, cscv_id, brdn_id, &
         tot_prod_id, tot_ddep_id, tot_lscv_id, tot_cscv_id, tot_brdn_id, &
         gridarea_id, &            !// Variable ID for gridarea
         dustnames_id              !// Variable ID for dustnames
    integer :: &
         srt_time(1), &            !// starting point for time array
         dim_lon_lat_time(3), &    !// Dimension ID for processes
         srt_lon_lat_time(3), &    !// Start array for lon/lat/time
         cnt_lon_lat_time(3), &       !// Counting array for lon/lat/time
         dim_nbr_lon_lat_time(4), &    !// Dimension ID for processes
         srt_nbr_lon_lat_time(4), &    !// Start array for nbr/lon/lat/time
         cnt_nbr_lon_lat_time(4)       !// Counting array for nbr/lon/lat/time

    !// Other locals
    integer :: &
         ncid, &                   !// File ID for nc file
         status, &                 !// Error status for nc file
         nlons, &                  !// Number of longitudes found in file
         nlats, &                  !// Number of latitudes found in file
         nbrs, &                   !// Number of tracers found in file
         nsteps                    !// Number of steps found in file 
    character(len=4)      :: cyear                  !// Year in character*4
    character(len=3)      :: cday                   !// Day in character*3
    character(len=80)     :: time_label             !// Label for variable "time"
    real(r8)              :: time                   !// Time in this timestep
    character(len=:), allocatable :: dustnames(:)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'dustbdg2file'
    !// --------------------------------------------------------------------

    !// Get budgets from STTTND (6=LSSCAV, 8=C_SCAV)
    !// These tendencies are summed up (kg/gridbox) from NDAY0 to NDAY
    davgtnd = 86400._r8 * real(NDAY+1 - NDAY0, r8)

    !// Check the time vs time for other diagnoses
    davgdays = dustbudgetacctime
    if (davgdays .ne. davgtnd) then
       write(6,'(a,2f12.3,2i8)') f90file//':'//subr//': '// &
            'wrong dt',davgdays,davgtnd, NDAY0, NDAY
       stop
    end if

    !// Large scale wash out [kg] (accumulated from NDAY0 to NDAY)
    dust_lscv(:,:,:) = 0._r8
    do K = 1, Ndustbins
       N = dust_trsp_idx(K)
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                !// Save as positive values: STTTND is BTT-BTTBCK.
                dust_lscv(K,I,J) = dust_lscv(K,I,J) &
                     - STTTND(I,J,L,N,NTND_LSSCAV)
             end do
          end do
       end do
    end do

    !// Convective wash out [kg] (accumulated from NDAY0 to NDAY)
    dust_cscv(:,:,:) = 0._r8
    do K = 1, Ndustbins
      N = dust_trsp_idx(K)
      do MP = 1, MPBLK
        do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = 1, LPAR
              !// Save as positive values: CONVWASHOUT is neg for removal.
               dust_cscv(K,I,J) = dust_cscv(K,I,J) &
                    - CONVWASHOUT(L,N,II,JJ,MP)
            end do
          end do
        end do
      end do
    end do

    !// Generate average; i.e. divide by time duration of accumulation,
    !// which is davgtnd: [kg] -> [kg/s]
    dust_production(:,:,:) = dust_production(:,:,:) / davgtnd
    dust_drydep(:,:,:) = dust_drydep(:,:,:) / davgtnd
    dust_lscv(:,:,:) = dust_lscv(:,:,:) / davgtnd
    dust_cscv(:,:,:) = dust_cscv(:,:,:) / davgtnd
    !// For burden, this means [kg*s] -> [kg]
    dust_burden(:,:,:) = dust_burden(:,:,:) / davgtnd

    !// Average prod/ddep/lscv/cscv/burden during period (already divided by dt)
    prod   = sum(dust_production)
    ddep   = -sum(dust_drydep)
    lscv   = -sum(dust_lscv)
    cscv   = -sum(dust_cscv)
    burden = sum(dust_burden)

    write(6,'(a,6es10.2)') 'DUST AVG P/D/LS/CN[kg/s]/B[kg]', &
         prod, ddep, lscv, cscv, burden
    write(6,'(a,2f11.3)')  'DUST LIFE B/P, B/L [days]', &
         burden/prod / 86400._r8,&
         -burden/(ddep+lscv+cscv) / 86400._r8

    !// Prepare for output
    !// --------------------------------------------------------------------

    !// For totals tendencies [kg/s] and brdn [kg]
    do J = 1, JPAR
       do I = 1, IPAR
          dstTotBdgProd(I,J) = sum(dust_production(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          dstTotBdgDdep(I,J) = sum(dust_drydep(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          dstTotBdgLscv(I,J) = sum(dust_lscv(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          dstTotBdgCscv(I,J) = sum(dust_cscv(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          dstTotBdgBrdn(I,J) = sum(dust_burden(:,I,J))
       end do
    end do

    !// Divide by area [kg] -> [kg/m2]
    ZAREA(:,:) = 1._r8 / AREAXY(:,:)
    do N = 1, ndustbins
       dust_production(N,:,:) = dust_production(N,:,:) * ZAREA(:,:)
       dust_drydep(N,:,:) = dust_drydep(N,:,:) * ZAREA(:,:)
       dust_lscv(N,:,:) = dust_lscv(N,:,:) * ZAREA(:,:)
       dust_cscv(N,:,:) = dust_cscv(N,:,:) * ZAREA(:,:)
       dust_burden(N,:,:) = dust_burden(N,:,:) * ZAREA(:,:)
    end do
    !// For totals and brdn
    dstTotBdgProd(:,:) = dstTotBdgProd(:,:) * ZAREA(:,:)
    dstTotBdgDdep(:,:) = dstTotBdgDdep(:,:) * ZAREA(:,:)
    dstTotBdgLscv(:,:) = dstTotBdgLscv(:,:) * ZAREA(:,:)
    dstTotBdgCscv(:,:) = dstTotBdgCscv(:,:) * ZAREA(:,:)
    dstTotBdgBrdn(:,:) = dstTotBdgBrdn(:,:) * ZAREA(:,:)

    !// The total time span for the average [days]
    davgdays = davgdays / 86400._r8

    !// Write next item to file (it is initialized as 0)
    dustaverages_written = dustaverages_written + 1


    !// --------------------------------------------------------------------
    !// Generate netCDF file
    status = nf90_noerr

    !// First time routine is called, we need to define variables
    if (dustaverages_written .eq. 1) then
       !// Create file
       write(6,'(a)') f90file//':'//subr//': creating '//trim(dustbdgfile)
       !// Clobber means you can overwrite existing data
       status=nf90_create(path=dustbdgfile,cmode=nf90_clobber,ncid=ncid)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error opening '//trim(dustbdgfile))

       !// File headers. The group should agree on these so that output from
       !// our model always has the same headers.
       status = nf90_put_att(ncid,nf90_global,'title', &
            trim(MODEL)//' dust averaged fields')
       if (status .ne. nf90_noerr) call handle_error(status,'attributing title')

       !// Define dimensions
       status=nf90_def_dim(ncid,'lon',IPAR,lon_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim lon')
       status=nf90_def_dim(ncid,'lat',JPAR,lat_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim lat')
       status=nf90_def_dim(ncid,'dst_nbr',ndustbins,nbr_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim dst_nbr')
       strLen = len(TNAME(dust_trsp_idx(1))) !// Length of TNAMEs
       status=nf90_def_dim(ncid,'str_len',strLen,str_len_id)
       if (status .ne. nf90_noerr) call &
            handle_error(status,'defining dim strLen')
       status=nf90_def_dim(ncid,'time',nf90_unlimited,time_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim time')

       !// Defining the combined id for a field (lon /lat /time)
       dim_lon_lat_time(1)=lon_dim_id
       dim_lon_lat_time(2)=lat_dim_id
       dim_lon_lat_time(3)=time_dim_id

       dim_nbr_lon_lat_time(1)=nbr_dim_id
       dim_nbr_lon_lat_time(2)=lon_dim_id
       dim_nbr_lon_lat_time(3)=lat_dim_id
       dim_nbr_lon_lat_time(4)=time_dim_id

       !// Defining the lon/lat/time-variable
       status = nf90_def_var(ncid,'lon',nf90_double,lon_dim_id,lon_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining lon')
       status = nf90_def_var(ncid,'lat',nf90_double,lat_dim_id,lat_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining lat')
       status = nf90_def_var(ncid,'time',nf90_double,time_dim_id,time_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining time')
       !// Defining average time span
       status=nf90_def_var(ncid,'timespan',nf90_double,time_dim_id,timespan_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining timespan')
       !// Defining dust names
       status = nf90_def_var(ncid,'dustname',nf90_char, &
                                  (/str_len_id,nbr_dim_id/),dustnames_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining dustname')
       !// Defining grid box area
       status=nf90_def_var(ncid,'gridarea',nf90_double, &
                                (/lon_dim_id,lat_dim_id/),gridarea_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining gridarea')

       !// Putting attributes to /lon/lat/ variables
       !// Don't change the units. It is important to keep these for GRADS
       status=nf90_put_att(ncid,lon_id,'units','degree_east')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit lon')
       status=nf90_put_att(ncid,lat_id,'units','degree_north')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit lat')

       write(cyear,'(I4)')IYEAR
       write(cday,'(I3)')NDAYI
       time_label='days since '//cyear//': day number '//trim(cday)

       !// Putting units attribute to time variable
       status=nf90_put_att(ncid,time_id,'units',time_label)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit time')
       !// Average time span
       status=nf90_put_att(ncid,timespan_id,'units','Average time span (days)')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit timespan')
       !// Gridbox area
       status=nf90_put_att(ncid,gridarea_id,'units','Grid box area [m2]')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit gridarea')
       !// Dust names
       status=nf90_put_att(ncid,dustnames_id,'units','Dust tracer names')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit dustnames')

       !// Define the tracer field variables for totals and their units
       status = nf90_def_var(ncid,'dustprod',nf90_float,dim_lon_lat_time,tot_prod_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining dustprod')
       status = nf90_put_att(ncid,tot_prod_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit dustprod')
       status = nf90_put_att(ncid,tot_prod_id,'longname','Total dust production')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname dustprod')

       status = nf90_def_var(ncid,'dustddep',nf90_float,dim_lon_lat_time,tot_ddep_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining dustddep')
       status = nf90_put_att(ncid,tot_ddep_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit dustddep')
       status = nf90_put_att(ncid,tot_ddep_id,'longname','Total dust drydep')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname dustddep')

       status = nf90_def_var(ncid,'dustwdep',nf90_float,dim_lon_lat_time,tot_lscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining dustwdep')
       status = nf90_put_att(ncid,tot_lscv_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit dustwdep')
       status = nf90_put_att(ncid,tot_lscv_id,'longname', &
            'Total dust largescale scavenged')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname dustwdep')

       status = nf90_def_var(ncid,'dustcnvw',nf90_float,dim_lon_lat_time,tot_cscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining dustcnvw')
       status = nf90_put_att(ncid,tot_cscv_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit dustcnvw')
       status = nf90_put_att(ncid,tot_cscv_id,'longname', &
            'Total dust convective scavenged')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname dustcnvw')

       status = nf90_def_var(ncid,'dustbrdn',nf90_float,dim_lon_lat_time,tot_brdn_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining dustbrdn')
       status = nf90_put_att(ncid,tot_brdn_id,'units','kg/m2')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit dustbrdn')
       status = nf90_put_att(ncid,tot_brdn_id,'longname','Total dust burden')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname dustbrdn')


       !// Define the tracer field variables for all DUST tracers
       if (LDUSTDIAG3D) then
          status = nf90_def_var(ncid,'dust_prod',nf90_float, &
               dim_nbr_lon_lat_time, prod_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining dust_prod')
          status = nf90_put_att(ncid,prod_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit dust_prod')
          status = nf90_put_att(ncid,prod_id,'longname', &
               'Dust production, all dust tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname dust_prod')

          status = nf90_def_var(ncid,'dust_ddep',nf90_float, &
               dim_nbr_lon_lat_time, ddep_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining dust_ddep')
          status = nf90_put_att(ncid,ddep_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit dust_ddep')
          status = nf90_put_att(ncid,ddep_id,'longname', &
               'Dust drydep, all dust tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname dust_ddep')

          status = nf90_def_var(ncid,'dust_lscv',nf90_float, &
               dim_nbr_lon_lat_time, lscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining dust_lscv')
          status = nf90_put_att(ncid,lscv_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit dust_lscv')
          status = nf90_put_att(ncid,lscv_id,'longname', &
               'Dust largescale scavenged, all dust tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname dust_lscv')

          status = nf90_def_var(ncid,'dust_cscv',nf90_float, &
               dim_nbr_lon_lat_time, cscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining dust_cscv')
          status = nf90_put_att(ncid,cscv_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit dust_cscv')
          status = nf90_put_att(ncid,cscv_id,'longname', &
               'Dust convective scavenged, all dust tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname dust_cscv')

          status = nf90_def_var(ncid,'dust_brdn',nf90_float, &
               dim_nbr_lon_lat_time, brdn_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining dust_brdn')
          status = nf90_put_att(ncid,brdn_id,'units','kg/m2')
          if (status .ne. nf90_noerr)&
               call handle_error(status,'error unit dust_brdn')
          status = nf90_put_att(ncid,brdn_id,'longname', &
               'Dust burden, all dust tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname dust_brdn')
       end if !// if (LDUSTDIAG3D) then

       !// End definition mode
       status = nf90_enddef(ncid)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error nf90_enddef')
 
       !// Putting the lon/lat/areaxy variables
       status = nf90_put_var(ncid,lon_id,XDGRD)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting lon')
       status = nf90_put_var(ncid,lat_id,YDGRD)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting lat')
       status = nf90_put_var(ncid,gridarea_id,AREAXY)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting gridarea')
       !// dust names
       allocate( character(strLen) :: dustnames(ndustbins) )
       do N = 1, ndustbins
          dustnames(N) = TNAME(dust_trsp_idx(N))
       end do
       status = nf90_put_var(ncid,dustnames_id,dustnames)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dustnames')
       deallocate( dustnames )
       write(6,'(a)') f90file//':'//subr//': created '//trim(dustbdgfile)

    else !// THE FILE HAS BEEN USED BEFORE
       !// Open the existing file
       write(6,'(a)') f90file//':'//subr//': opening '//trim(dustbdgfile)
       status = nf90_open(dustbdgfile, nf90_write, ncid)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error in open')
       !// Inquire dimension ids
       status = nf90_inq_dimid(ncid,'lat',lat_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting lat_id')
       status = nf90_inq_dimid(ncid,'lon',lon_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting lon_id')
       status = nf90_inq_dimid(ncid,'dst_nbr',nbr_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting nbr_id')
       status = nf90_inq_dimid(ncid,'time',time_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting time_id')
     
       !// Inquire dimensions
       status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim lat')
       if (nlats .ne. JPAR)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': JPAR mismatch on existing file: ',nlats,JPAR
          stop
       endif
       status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim lon')
       if (nlons .ne. IPAR)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': IPAR mismatch on existing file: ',nlons,IPAR
          stop
       endif
       status = nf90_Inquire_Dimension(ncid,nbr_dim_id,len=nbrs)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim nbr')
       if (nbrs .ne. ndustbins)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': nbrs mismatch on existing file: ',nbrs,ndustbins
          stop
       endif
       status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim time')
       if (nsteps+1 .ne. dustaverages_written)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': nsteps mismatch on existing file: ', &
               nsteps, dustaverages_written
          stop
       endif

       !// Get ids for totals
       status = nf90_inq_varid(ncid,'dustprod',tot_prod_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_prod_id')
       status = nf90_inq_varid(ncid,'dustddep',tot_ddep_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_ddep_id')
       status = nf90_inq_varid(ncid,'dustwdep',tot_lscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_lscv_id')
       status = nf90_inq_varid(ncid,'dustcnvw',tot_cscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_cscv_id')
       status = nf90_inq_varid(ncid,'dustbrdn',tot_brdn_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_brdn_id')

       !// Get ids for all dust bins
       if (LDUSTDIAG3D) then
          status = nf90_inq_varid(ncid,'dust_prod',prod_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting prod_id')
          status = nf90_inq_varid(ncid,'dust_ddep',ddep_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting ddep_id')
          status = nf90_inq_varid(ncid,'dust_lscv',lscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting lscv_id')
          status = nf90_inq_varid(ncid,'dust_cscv',cscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting cscv_id')
          status = nf90_inq_varid(ncid,'dust_brdn',brdn_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting brdn_id')
       end if !// if (LDUSTDIAG3D) then

       !// Get variable id for time
       status=nf90_inq_varid(ncid,'time',time_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting time_id')
       status = nf90_inq_varid(ncid,'timespan',timespan_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting timespan_id')

    end if !// Has the file already been in use or not ?

    !// --------------------------------------------------------------------
    !// Defining how far to count for each time a data set is added
    !// Adding 2D fields
    cnt_lon_lat_time = (/IPAR , JPAR,  1/)
    !// Defining where to start adding the new time step
    srt_lon_lat_time = (/1, 1, dustaverages_written/)

    !// Adding 3D fields
    cnt_nbr_lon_lat_time = (/ndustbins, IPAR , JPAR,  1/)
    !// Defining where to start adding the new time step
    srt_nbr_lon_lat_time = (/1, 1, 1, dustaverages_written/)


    srt_time(1) = dustaverages_written    !// Start value for new time step
    time = real(NDAY - NDAYI + 1, r8)     !// Time in r8 format (days)
  
    status = nf90_put_var(ncid,time_id,time,start=srt_time)
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting time')

    !// Also put out average time span
    status = nf90_put_var(ncid,timespan_id,davgdays,start=srt_time)
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting avg time span')


    !// Put tracer fields for totals
    status = nf90_put_var(ncid, tot_prod_id, dstTotBdgProd(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dustprod')

    status = nf90_put_var(ncid, tot_ddep_id, dstTotBdgDdep(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dustddep')

    status = nf90_put_var(ncid, tot_lscv_id, dstTotBdgLscv(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dustwdep')

    status = nf90_put_var(ncid, tot_cscv_id, dstTotBdgCscv(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dustcnvw')

    status = nf90_put_var(ncid, tot_brdn_id, dstTotBdgBrdn(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dustbrdn')


    !// Put fields for all tracers
    if (LDUSTDIAG3D) then
       !// Put out dust_prod
       status = nf90_put_var(ncid, prod_id, &  !// File id, variable id
            dust_production(:,:,:), &       !// Dust produced
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dust_prod')

       !// Put out dust_ddep
       status = nf90_put_var(ncid, ddep_id, &  !// File id, variable id
            dust_drydep(:,:,:), &           !// Dust drydepped
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dust_ddep')

       !// Put out dust_lscv
       status = nf90_put_var(ncid, lscv_id, &  !// File id, variable id
            dust_lscv(:,:,:), &             !// Dust ls scav
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dust_lscv')

       !// Put out dust_cscv
       status = nf90_put_var(ncid, cscv_id, &  !// File id, variable id
            dust_cscv(:,:,:), &             !// Dust cnv scav
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dust_cscv')

       !// Put out dust_brdn
       status = nf90_put_var(ncid, brdn_id, &  !// File id, variable id
            dust_burden(:,:,:), &           !// Dust burden
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting dust_brdn')
    end if !// if (LDUSTDIAG3D) then
  
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)
    write(6,'(a)') f90file//':'//subr//': updated '//trim(dustbdgfile)


    !// Re-initialize
    dustbudgetacctime = 0._r8

    dust_production(:,:,:) = 0._r8
    dust_drydep(:,:,:)     = 0._r8
    dust_lscv(:,:,:)       = 0._r8
    dust_cscv(:,:,:)       = 0._r8
    dust_burden(:,:,:)     = 0._r8

    !// Re-initialize last-values when dustbudget is zeroed
    lastburden = 0._r8
    lastlscv = 0._r8
    lastcscv = 0._r8
    lastprod = 0._r8
    lastddep = 0._r8
    lastchem = 0._r8

    !// --------------------------------------------------------------------
  end subroutine dustbdg2file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
end module dust_oslo
!//=========================================================================
