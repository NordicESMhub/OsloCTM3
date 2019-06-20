!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Amund Sovde Haslerud, May 2016, April 2015
!//=========================================================================
!// CTM diagnostics for scavenging (dry & wet)
!//=========================================================================
module diagnostics_scavenging
  !//-----------------------------------------------------------------------
  !// MODULE: diagnostics_scavenging
  !// DESCRIPTION: Routines for scavenging diagnostics.
  !// ----------------------------------------------------------------------
  !// Diagnostics for the Oslo CTM3.
  !// Contains:
  !//   subroutine scav_diag_init
  !//   subroutine scav_diag_put_ddep
  !//   subroutine scav_diag_brd
  !//   subroutine scav_diag_ls
  !//   subroutine scav_diag_cn
  !//   subroutine scav_diag_collect_daily
  !//   subroutine scav_diag_2fileA
  !//   subroutine scav_diag_2fileB
  !//   subroutine scav_diag_put_gsto
  !//   subroutine scav_diag_put_fsto
  !//   subroutine scav_diag_nmet_output_nc
  !//
  !// Stefanie Falk, July 2018
  !//   Added output and computation of daily stomatal conductance and flux.
  !//   Added output of stomatal conductance and flux each timestep (NMETxNOPS).
  !// Ole Amund Sovde, April 2015
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, MPBLK, IDBLK, JDBLK
  use cmn_oslo, only: SCAV_LS, SCAV_CN, SCAV_DD, SCAV_BRD, SCAV_DIAG, &
       SCAV_MAP_WLS, SCAV_MAP_WCN, SCAV_MAP_DRY, &
       GSTO3_AVG, FSTO3_AVG, VRAO3_AVG, VRBO3_AVG, VRCO3_AVG
  use cmn_sfc, only: VGSTO3, NLCAT
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'diagnostics_scavenging.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public scav_diag_init, scav_diag_put_ddep, scav_diag_brd, scav_diag_ls, &
       scav_diag_cn, scav_diag_collect_daily, &
       scav_diag_2fileA, scav_diag_2fileB, &
       scav_diag_put_gsto, scav_diag_put_drydepvelo, &
       scav_diag_put_fsto, scav_diag_nmet_output_nc
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine scav_diag_init(FLAG)
    !// --------------------------------------------------------------------
    !// Initialize diagnose arrays.
    !// Called from init_oslo.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: FLAG
    !// --------------------------------------------------------------------

    !// IJ-arrays
    SCAV_LS(:,:) = 0._r8
    SCAV_CN(:,:) = 0._r8
    SCAV_DD(:,:) = 0._r8
    SCAV_BRD(:,:) = 0._r8
    !// Global daily (should only be initialized at model start)
    if (FLAG) SCAV_DIAG(:,:,:) = 0._r8
    !// Maps
    SCAV_MAP_WLS(:,:,:,:) = 0._r8
    SCAV_MAP_WCN(:,:,:,:) = 0._r8
    SCAV_MAP_DRY(:,:,:,:) = 0._r8
    !// Dry deposition velocities
    GSTO3_AVG(:,:,:) = 0._r8
    VRAO3_AVG(:,:,:) = 0._r8
    VRBO3_AVG(:,:,:) = 0._r8
    VRCO3_AVG(:,:,:) = 0._r8
    !// Leaf levelStomatal flux
    FSTO3_AVG(:,:,:) = 0._r8

    !// --------------------------------------------------------------------
  end subroutine scav_diag_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine scav_diag_put_ddep(DDDIAG,DV,II,JJ,MP)
    !// --------------------------------------------------------------------
    !// Save drydep amount [kg] from chemistry. Values are
    !// accumulated over the chemistry time step, and the units
    !// of DDDIAG are [molec/cm3].
    !//
    !// Ole Amund Sovde, March 2014
    !// --------------------------------------------------------------------
    use cmn_size, only: TRACER_ID_MAX
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    use cmn_oslo, only: chem_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: DDDIAG(TRACER_ID_MAX), DV
    integer,intent(in) :: II,JJ,MP
    !// Locals
    real(r8)  :: rfac
    integer :: N
    !// --------------------------------------------------------------------
    do N = 1, NPAR
       rfac = DV * 1.e6_r8 / AVOGNR * TMASS(N) * 1.e-3_r8
       !// DDDIAG can be zero. Only collect for transported species.
       SCAV_MAP_DRY(N,II,JJ,MP) = SCAV_MAP_DRY(N,II,JJ,MP) &
            + DDDIAG(chem_idx(N)) * rfac !// molecules/cm3 > kg/gridbox
       !// Save totals
       SCAV_DD(N,MP) = SCAV_DD(N,MP) + DDDIAG(chem_idx(N)) * rfac 
    end do
    !// --------------------------------------------------------------------
  end subroutine scav_diag_put_ddep
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine scav_diag_brd(NSUB, NOPS, NMET, BTT,MP)
    !// --------------------------------------------------------------------
    !// Save tracer masses each NOPS; an average will be calculated
    !// at the end of each day. These burdens can be used to calculate
    !// lifetimes.
    !//
    !// Ole Amund Sovde, June 2013
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT
    integer, intent(in) :: NSUB, NOPS, NMET, MP
    !// Locals
    integer :: N 
    !// --------------------------------------------------------------------
    if (NSUB.ne.1) return
    !// Add up burden to make average
    do N = 1, NPAR
       SCAV_BRD(N,MP) = SCAV_BRD(N,MP) + sum(BTT(:,N,:,:))
    end do
    !// --------------------------------------------------------------------
  end subroutine scav_diag_brd
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine scav_diag_ls(NMET,BTT,BTTBCK,MP)
    !// --------------------------------------------------------------------
    !// Save tracer masses lost due to large scale wet scavenging.
    !// Accumulate all time steps though the day.
    !//
    !// Ole Amund Sovde, June 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: LBCOC, LSOA
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TCHENA
    use bcoc_oslo, only: bcsnow_diagwetrm
    use soa_oslo, only: soa_diag_lsscav
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT, BTTBCK
    integer, intent(in) :: NMET, MP
    !// Locals
    integer :: I,J,II,JJ,N, NN, NWET, NTNR(NPAR)
    real(r8) :: tdep
    !// --------------------------------------------------------------------
    !// Find which components are washed out
    NTNR(:) = 0
    NWET = 0
    do N = 1, NPAR
       if (TCHENA(N) .gt. 0._r8) then
          NWET = NWET + 1
          NTNR(NWET) = N
       end if
    end do
    if (NWET .eq. 0) return

    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1
        do NN = 1, NWET
          N = NTNR(NN)
          tdep = sum(BTTBCK(:,N,II,JJ) - BTT(:,N,II,JJ))
          !// Add to map
          SCAV_MAP_WLS(N,II,JJ,MP) = SCAV_MAP_WLS(N,II,JJ,MP) + tdep
          !// Add to totals
          SCAV_LS(N,MP) = SCAV_LS(N,MP) + tdep
        end do
      end do
    end do

    !// BCsnow - diagnose scavenged BC
    if (LBCOC) call bcsnow_diagwetrm(BTT,BTTBCK,MP)

    !// SOA diagnose large scale scav in double precision
    if (LSOA) call soa_diag_lsscav(BTT,BTTBCK,MP)

    !// --------------------------------------------------------------------
  end subroutine scav_diag_ls
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine scav_diag_cn(NMET,BTT,BTTBCK,MP)
    !// --------------------------------------------------------------------
    !// Save tracer masses lost due to convective wet scavenging.
    !// Accumulate all time steps though the day.
    !//
    !// Ole Amund Sovde, June 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: LBCOC
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TCHENA
    use bcoc_oslo, only: bcsnow_diagwetrm
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT, BTTBCK
    integer, intent(in) :: NMET, MP
    !// Locals
    integer :: I,J,II,JJ,N, NN, NWET, NTNR(NPAR)
    real(r8) :: tdep
    !// --------------------------------------------------------------------
    !// Find which components are washed out
    NTNR(:) = 0
    NWET = 0
    do N = 1, NPAR
       if (TCHENA(N) .gt. 0._r8) then
          NWET = NWET + 1
          NTNR(NWET) = N
       end if
    end do
    if (NWET .eq. 0) return

    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1
        do NN = 1, NWET
          N = NTNR(NN)
          tdep = sum(BTTBCK(:,N,II,JJ) - BTT(:,N,II,JJ))
          !// Add to map
          SCAV_MAP_WCN(N,II,JJ,MP) = SCAV_MAP_WCN(N,II,JJ,MP) + tdep
          !// Add to totals
          SCAV_CN(N,MP) = SCAV_CN(N,MP) + tdep
        end do
      end do
    end do

    !// BCsnow - diagnose scavenged BC
    if (LBCOC) call bcsnow_diagwetrm(BTT,BTTBCK,MP)

    !// --------------------------------------------------------------------
  end subroutine scav_diag_cn
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine scav_diag_collect_daily(KDAY)
    !// --------------------------------------------------------------------
    !// Collect daily scavenged tracer masses of large scale and convective
    !// scavenging, and also daily average tracer burden.
    !//
    !// SCAV_LS, SCAV_CN and SCAV_BRD are accumulated data, and zeroed
    !// only when BUDGETS calendar initiates writing to file
    !// (typically at the beginning of a month).
    !// The SCAV_DIAG is thus accumulated, but corrected for it when
    !// written to file.
    !//
    !// Ole Amund Sovde, October 2014, June 2013
    !// --------------------------------------------------------------------
    use cmn_ctm, only: NROPSM, NRMETD, STT
    use cmn_chem, only: TNAME
    use cmn_diag, only: JDO_T
    use utilities, only: get_free_fileid
    use cmn_oslo, only: DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: KDAY
    !// Locals
    real(r8) :: scav(NPAR), burden_all(NPAR), t1, t2, t3, t4
    integer :: MP, N, TD, PDAY
    !// --------------------------------------------------------------------

    !// Wet scavenging for LS
    scav(:) = 0._r8
    do MP = 1, MPBLK
       do N = 1, NPAR
          scav(N) = scav(N) + SCAV_LS(N,MP)
       end do
    end do
    SCAV_DIAG(:,1,KDAY) = scav(:)

    !// Wet scavenging for CNV
    scav(:) = 0._r8
    do MP = 1, MPBLK
       do N = 1, NPAR
          scav(N) = scav(N) + SCAV_CN(N,MP)
       end do
    end do
    SCAV_DIAG(:,2,KDAY) = scav(:)

    !// Dry deposition (scavenging)
    scav(:) = 0._r8
    do MP = 1, MPBLK
       do N = 1, NPAR
          scav(N) = scav(N) + SCAV_DD(N,MP)
       end do
    end do
    SCAV_DIAG(:,3,KDAY) = scav(:)

    !// Accumulated burden
    burden_all(:) = 0._r8
    do MP = 1, MPBLK
       do N = 1, NPAR
          !// Divide by (NROPSM*NRMETD) to get accumulated daily mean
          burden_all(N) = burden_all(N) &
               + SCAV_BRD(N,MP) / real(NROPSM*NRMETD, r8)
       end do
    end do
    SCAV_DIAG(:,4,KDAY) = burden_all(:)

    if (KDAY .gt. 59 .and. DINM(2) .eq. 28._r8) then
       !// Not leap year, so we have to jump past 60 in JDO_T
       TD = KDAY+1
    else
       TD = KDAY
    end if
    !// Daily average of dry deposition velocites -> devide by (NROPSM*NRMETD)
    VRAO3_AVG = VRAO3_AVG / (NROPSM*NRMETD)
    VRBO3_AVG = VRBO3_AVG / (NROPSM*NRMETD)
    VRCO3_AVG = VRCO3_AVG / (NROPSM*NRMETD)
    !// Daily average of stomatal conductance -> devide by (NROPSM*NRMETD)
    GSTO3_AVG = GSTO3_AVG / (NROPSM*NRMETD)
    !// Daily average of stomatal flux -> devide by (NROPSM*NRMETD)
    FSTO3_AVG = FSTO3_AVG / (NROPSM*NRMETD)

    !// Previous day
    if (KDAY .eq. 1) then
       !// End of last year
       PDAY = int(sum(DINM))
    else
       PDAY = KDAY - 1
    end if
    !// If JDO_T is not set at 1 Jan, and the run continues from
    !// last year, we should subtract what is accumulated at the
    !// end of last year. When starting at 1Jan, this value is zero,
    !// so we can only test on JDO_T.

    do N = 1, NPAR
!       if (KDAY .eq. 1) then
!          t1 = scav_diag(N,1,KDAY)
!          t2 = scav_diag(N,2,KDAY)
!          t3 = scav_diag(N,3,KDAY)
!          t4 = scav_diag(N,4,KDAY)
!       else if (JDO_T(TD) .eq. 0) then
       if (JDO_T(TD) .eq. 0) then
          !// KDAY is current day, before updating calendar.
          !// Which means we have to subtract accumulated data
          !// from previous day to get daily accumulated data.
          t1 = scav_diag(N,1,KDAY) - scav_diag(N,1,PDAY)
          t2 = scav_diag(N,2,KDAY) - scav_diag(N,2,PDAY)
          t3 = scav_diag(N,3,KDAY) - scav_diag(N,3,PDAY)
          t4 = scav_diag(N,4,KDAY) - scav_diag(N,4,PDAY)
       else
          !// JDO_T is set; accumulated arrays were zeroed at the
          !// end of the previous day, so we must not subtract
          !// scav_diag for previous day.
          t1 = scav_diag(N,1,KDAY)
          t2 = scav_diag(N,2,KDAY)
          t3 = scav_diag(N,3,KDAY)
          t4 = scav_diag(N,4,KDAY)
       end if

       if ((t1+t2+t3) .gt. 0._r8) &
            write(6,'(a,4es10.2,es12.4)') &
            'Daily lifetimes L/C/L+C/D [days] / B [kg] '// &
            tname(n)//' ', t4/t1, t4/t2, t4/(t1+t2), t4/t3, t4

    end do

    !// --------------------------------------------------------------------
  end subroutine scav_diag_collect_daily
  !// ----------------------------------------------------------------------
 !// ----------------------------------------------------------------------
  subroutine scav_diag_put_drydepvelo(NMET,VRAO3,VRBO3,VRCO3)
    !// --------------------------------------------------------------------
    !// Compute the daily average dry deposition velocitites.
    !//
    !// VRaO3, VRbO3, and VRcO3 arein units of ms-1.
    !//
    !// Stefanie Falk, July 2018
    !// --------------------------------------------------------------------
    
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), dimension(IPAR,JPAR,NLCAT), intent(in) :: VRAO3,VRBO3,VRCO3
    integer,intent(in) :: NMET
    
    !// --------------------------------------------------------------------
    !// Collect totals
    !// average will be calculated in scav_diag_collect_daily
    VRAO3_AVG = VRAO3_AVG + VRAO3
    VRBO3_AVG = VRBO3_AVG + VRBO3
    VRCO3_AVG = VRCO3_AVG + VRCO3
    !// --------------------------------------------------------------------
  end subroutine scav_diag_put_drydepvelo
  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------
  subroutine scav_diag_put_gsto(NMET,VGSTO3)
    !// --------------------------------------------------------------------
    !// Compute the daily average stomatal uptake from VGstO3.
    !//
    !// VGstO3 is in units of mmol m-2s-1.
    !//
    !// Stefanie Falk, July 2018
    !// --------------------------------------------------------------------
    
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), dimension(IPAR,JPAR,NLCAT), intent(in) :: VGSTO3
    integer,intent(in) :: NMET
    
    !// --------------------------------------------------------------------
    !// Collect totals
    !// average will be calculated in scav_diag_collect_daily
    GSTO3_AVG = GSTO3_AVG + VGSTO3
    !// --------------------------------------------------------------------
  end subroutine scav_diag_put_gsto
  !// ----------------------------------------------------------------------

!// ----------------------------------------------------------------------
  subroutine scav_diag_put_fsto(BTT,NMET,VGSTO3,MP)
    !// --------------------------------------------------------------------
    !// Compute the daily average stomatal uptake from FGstO3.
    !//
    !// FGstO3 is in units of mmol m-2s-1, but it will be changed 
    !// to nmol m-2s-1 when written to output file.
    !//
    !// Stefanie Falk, July 2018
    !// --------------------------------------------------------------------
    !use cmn_size, only: IDBLK, JDBLK, NPAR, LPAR
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_oslo, only: AIRMOLEC_IJ, DV_IJ
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), dimension(IPAR,JPAR,NLCAT), intent(in) :: VGSTO3
    integer,  intent(in) :: NMET, MP
    real(r8), intent(in)  :: BTT(LPAR,NPAR,IDBLK,JDBLK)
    !// For looping
    integer :: I, J, II, JJ
    real(r8) :: RDUM, VMR_O3
    !// Local
    real(r8), dimension(IPAR,JPAR,NLCAT) :: fsto3_tmp
    !// --------------------------------------------------------------------
    
    !// Initialize fsto3_tmp
    fsto3_tmp(:,:,:) = 0._r8 
    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// DV_IJ and AIRMOLEC_IJ are zero during the first call(?)
          !// Prevent NAN-spread
          if ( DV_IJ(1,II,JJ,MP) .eq. 0._r8 .or. AIRMOLEC_IJ(1,II,JJ,MP).eq. 0._r8 ) then
             fsto3_tmp(I,J,:) = 0._r8
          else
             !// Convert O3 [kg/gridbox] to [molec/cm3]
             RDUM =  1.e-3_r8 * AVOGNR / TMASS(1) / DV_IJ(1,II,JJ,MP)
             !// Compute stomatal flux by multiplying conductance with 
             !// VMR_O3 at ground level
             VMR_O3 = (BTT(1,1,II,JJ)*RDUM)/AIRMOLEC_IJ(1,II,JJ,MP)
             !// fsto3_tmp(I,J) units [mmol m-2 s-1]
             fsto3_tmp(I,J,:) = VGSTO3(I,J,:)*VMR_O3
          end if
       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    !// --------------------------------------------------------------------
    !// Collect totals
    !// average will be calculated in scav_diag_collect_daily
    FSTO3_AVG = FSTO3_AVG + fsto3_tmp
   
    !// --------------------------------------------------------------------
  end subroutine scav_diag_put_fsto
  !// ----------------------------------------------------------------------
  
  !// ----------------------------------------------------------------------
  subroutine scav_diag_nmet_output_nc(JYEAR,JMONTH,JDATE,NDAY,NMET,NOPS)
    !// --------------------------------------------------------------------
    !// Called from pmain before the calendar update.
    !//
    !// Write values of FSTO3 and GSTO3 at each NMET x NOPS
    !// to netcdf file.
    !//
    !// FSTO3_AVG and GSTO3_AVG are accumulated. 
    !// Both will be to deaccumulate before they are written to file.
    !// FSTO3_AVG and GSTO3_AVG are in units of mmol m-2 s-1.
    !//
    !// Stefanie Falk, July 2018
    !// --------------------------------------------------------------------
    use netcdf
    use cmn_diag, only: nc4deflate_global, nc4shuffle_global
    use cmn_ctm, only: XDGRD, YDGRD, ZGRD, NRMETD
    use cmn_met, only: MYEAR
    use cmn_sfc, only: LDDEPmOSaic
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JYEAR, JMONTH, JDATE, NDAY, NMET, NOPS
    !// Locals
    character(len=80) :: filename
    character(len=3)  :: cday           !Name of start day in character
    character(len=2)  :: cmonth         !Name of start month in character
    character(len=4)  :: cyear          !Name of start year in character
    !// Dimensions ID and counters
    integer :: lat_dim_id               !Dimension ID for latitude
    integer :: lon_dim_id               !Dimension ID for longitude
    integer :: nlcat_dim_id             !Dimension ID for NLCAT
    integer :: time_dim_id              !Dimension ID for time
    integer :: lat_id                   !Variable ID for latitude
    integer :: lon_id                   !Variable ID for longitude
    integer :: nlcat_id                 !Variable ID for NLCAT
    integer :: time_id                  !Variable ID for time
    integer :: dim_lon_lat_nlcat_time(4)!Dimension ID for processes
    integer :: srt_lon_lat_nlcat_time(4)!Start array for lon/lat/nlcat/time
    integer :: cnt_lon_lat_nlcat_time(4)!Counting array for lon/lat/nlcat/time
    integer :: srt_time(1)              !starting point for time array
    integer :: gsto3_inst_id            !Variable ID for stomatal conductance
    integer :: fsto3_inst_id            !Variable ID for stomatal flux
  
    !// Other locals
    integer :: ncid                     !File ID for nc file
    integer :: status                   !Error status for nc file
    integer :: nlons                    !Number of longitudes found in file
    integer :: nlats                    !Number of latitudes found in file
    integer :: nsteps                   !Number of steps found in file 
    character(len=80) :: time_label     !Label for variable "time"
    real(r8):: time                     !Time in this timestep
    integer :: nbr_steps                !Step number
    integer :: ncats                    !Numer of land use categories found in file
    integer :: K
    integer, dimension(NLCAT) :: LCAT =  &   !Landuse category numbering
         (/(K, K=1,NLCAT, 1)/)         

    real(r8), save, dimension(IPAR,JPAR,NLCAT) :: GSTO3_priv=0._r8, FSTO3_priv=0._r8
    real(r8), dimension(IPAR,JPAR,NLCAT,1) :: GSTO3_inst, FSTO3_inst
    
    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'scav_diag_nmet_output_nc'
    !//---------------------------------------------------------------------
    !// No need to initialize if new dry deposition scheme is not used
    if (.not.LDDEPmOSaic) return
    
    !// Time step
    nbr_steps = NOPS + (NMET-1)* 24 / NRMETD !3
    
    !// Initialize and reset GSTO3_priv and FSTO_priv every first timestep
    if ((nbr_steps.eq.1)) then
       GSTO3_priv=0._r8
       FSTO3_priv=0._r8

    end if

    !//---------------------------------------------------------------------
    !// Deaccumulate GSTO3_AVG and FSTO_AVG
    write(6,'(a)') f90file//':'//subr// ': Deaccumulating stomatal fields'
  
    GSTO3_inst(:,:,:,1) = GSTO3_AVG-GSTO3_priv
    FSTO3_inst(:,:,:,1) = FSTO3_AVG-FSTO3_priv
    !// Save previous step
    GSTO3_priv = GSTO3_AVG
    FSTO3_priv = FSTO3_AVG

    !// Output
    !// [mmol m-2 s-1] -> [nmol m-2 s-1]
    FSTO3_inst = FSTO3_inst*1.e6_r8
   
    !//---------------------------------------------------------------------

    !// Filename
    write(cday(1:3),'(i3.3)') NDAY
    write(cyear(1:4),'(i4.4)') MYEAR
    filename = 'scavenging_daily_stomata_'//cyear//'_'//cday//'.nc'
    
    if (nbr_steps .eq. 1) then
       !First time this is done. Need to define variables
       write(6,'(a71)') '--------------------------------------------'// &
            '---------------------------'
       write(6,'(a)') f90file//':'//subr// &
            ' Will write data to netCDF every 1 hour'
       !// Create file
       write(6,*)'creating file: ',filename,status
       !// Use nf90_netcdf4 to create netcdf4:
       status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
       if (status/=nf90_noerr) call handle_error(status,'in creating file')

       !//File headers
       status=nf90_put_att(ncid,nf90_global,'title','1h output of stomata conductance/flux from Oslo CTM3')
       if (status/=nf90_noerr) call handle_error(status,f90file//':'//subr//': file header')
       status=nf90_put_att(ncid,nf90_global,'Modelinfo1','Oslo CTM3 is a 3D Chemical Transport Model')
       if (status/=nf90_noerr) call handle_error(status,f90file//':'//subr//': modelifo1')
       status=nf90_put_att(ncid,nf90_global,'Modelinfo2','Oslo CTM3 is driven by ECMWF meteorological data')
       if (status/=nf90_noerr) call handle_error(status,f90file//':'//subr//': modelinfo2')
       status=nf90_put_att(ncid,nf90_global,'contactinfo','For errors, contact CICERO')
       if (status/=nf90_noerr) call handle_error(status,f90file//':'//subr//': contactinfo')

       !// Define sizes
       !// Define dimensions (JM, IM, NLCAT, time)
       status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lat dim')
       status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lon dim')
       !// Defining nlcat dimension
       status = nf90_def_dim(ncid,"NLCAT",NLCAT,nlcat_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define NLCAT dim')
       !// Defining dimension time of length unlimited
       status = nf90_def_dim(ncid,"time",nf90_unlimited,time_dim_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define time dim')
      

       !// Defining the combined id for a field (lon /lat/ nlcat /time)
       dim_lon_lat_nlcat_time(1)=lon_dim_id
       dim_lon_lat_nlcat_time(2)=lat_dim_id
       dim_lon_lat_nlcat_time(3)=nlcat_dim_id
       dim_lon_lat_nlcat_time(4)=time_dim_id

       !// Defining the lon/lat/time-variable
       status = nf90_def_var(ncid,"lon",nf90_float,lon_dim_id,lon_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lon variable')
       status = nf90_def_var(ncid,"lat",nf90_float,lat_dim_id,lat_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lat variable')
       status = nf90_def_var(ncid,"time",nf90_float,time_dim_id,time_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define time variable')
       !// Defining the nlcat-variable
       status = nf90_def_var(ncid,"NLCAT",nf90_int,nlcat_dim_id,nlcat_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define nlcat variable')
       

       !// Putting attributes to /lon/lat variables
       status = nf90_put_att(ncid,lon_id,'units','degree_east')
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units lon')
       status = nf90_put_att(ncid,lat_id,'units','degree_north')
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units lat')

       !// Putting attributes to nlcat variable
       status = nf90_put_att(ncid,nlcat_id,'description', &
            'Landuse categories used in Oslo CTM3.' // char(10) // &
            '01 - Needleleaftree temp./bor.(CF)' // char(10) // &
            '02 - Deciduoustree temp./bor. (DF)' // char(10) // &
            '03 - Needleleaftree med. (NF)' // char(10) // &
            '04 - Broadleaftree (BF)' // char(10) // &
            '05 - Crops (TC)' // char(10) // &
            '06 - Moorland (SNL)' // char(10) // &
            '07 - Grassland (GR)' // char(10) // &
            '08 - Scrubs med. (MS)' // char(10) // &
            '09 - Wetlands (WE)' // char(10) // &
            '10 - Tundra (TU)' // char(10) // &
            '11 - Desert (DE)' // char(10) // &
            '12 - Water (W)' // char(10) // &
            '13 - Urban (U)' // char(10) // &
            '14 - Ice/Snow (ICE)') 
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description nlcat')  
      
       !// Putting attributes to time variable
       time_label='hours since yyyy-mm-dd 00:00:00'
       write(time_label(13:23),'(i4.4,a1,i2.2,a1,i2.2)') &
            MYEAR,'-',JMONTH,'-',JDATE
       status = nf90_put_att(ncid,time_id,'units',time_label)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units time')

       !// Daily average of stomatal conductance
       status = nf90_def_var(ncid,'GstO3_inst', &
            nf90_double, dim_lon_lat_nlcat_time, gsto3_inst_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable GstO3_inst ')
       status = nf90_def_var_deflate(ncid, gsto3_inst_id, &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable GstO3_inst ')
       status = nf90_put_att(ncid, gsto3_inst_id,'unit','mmol m-2 s-1')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute unit GstO3_inst ')

       !// Daily average of stomatal flux
       status = nf90_def_var(ncid,'FstO3_inst', &
            nf90_double, dim_lon_lat_nlcat_time, fsto3_inst_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable FstO3_inst ')
       status = nf90_def_var_deflate(ncid, fsto3_inst_id, &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable FstO3_inst ')
       status = nf90_put_att(ncid, fsto3_inst_id,'unit','nmol m-2 s-1')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute unit FstO3_inst ')

       !//---------------------------------------------------------------------
       !// End definition mode
       status = nf90_enddef(ncid)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': end defmode')
       !//---------------------------------------------------------------------
       !// Putting the lon/lat/nlcat variables
       status = nf90_put_var(ncid,lon_id,XDGRD)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting lon')
       status = nf90_put_var(ncid,lat_id,YDGRD)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting lat')
       !// Putting the NLCAT variables
       status = nf90_put_var(ncid,nlcat_id,LCAT)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting nlcat')
       
    else                                      !// THE FILE HAS BEEN USED BEFORE
       !// Open the existing file
       status = nf90_open(filename, nf90_write, ncid)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': open existing file')

       !// Inquire dimension ids
       status = nf90_inq_dimid(ncid,"lat",lat_dim_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': getting lat')
       status = nf90_inq_dimid(ncid,"lon",lon_dim_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': getting lon')
       status = nf90_inq_dimid(ncid,"NLCAT",nlcat_dim_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': getting NLCAT')
       status = nf90_inq_dimid(ncid,"time",time_dim_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': getting time')

       !// Inquire dimensions
       status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq lat dim')
       if (nlats/=JPAR) then
          write(6,*) f90file//':'//subr//': reports JM = ',nlats, JPAR
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq lon dim')
       if (nlons/=IPAR) then
          write(6,*) f90file//':'//subr//': reports IM = ',nlons,IPAR
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,nlcat_dim_id,len=ncats)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq NLCAT dim')
       if (ncats/=NLCAT) then
          write(6,*) f90file//':'//subr//': reports NLCAT = ',ncats, NLCAT
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq time dim')
       if (nsteps.lt.0) then
          write(6,*) f90file//':'//subr//': reports already added nsteps = ',nsteps
          stop
       end if

       !// Check variable id
       status = nf90_inq_varid(ncid,'GstO3_inst',gsto3_inst_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq tracer varid')
       status = nf90_inq_varid(ncid,'FstO3_inst',fsto3_inst_id)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq tracer varid')

       !// Get variable id for time
       status = nf90_inq_varid(ncid,"time",time_id) 
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq time varid')
    end if
    
    !// PUT THE VARIABLES TO FILE
    !// --------------------------------------------------------------------

    !// For the tracer fields:
    !// Defining how far to count for each time a data set is added
    cnt_lon_lat_nlcat_time = (/IPAR , JPAR , NLCAT, 1/)
    !// Defining where to start adding the new time step
    srt_lon_lat_nlcat_time = (/1, 1, 1, nbr_steps/)

    !// Start value for new time step
    srt_time(1) = nbr_steps
    time = nbr_steps-1

    status = nf90_put_var(ncid,time_id,time,start=srt_time)
    if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting time')

    !// Put variable tracer field
    status = nf90_put_var(ncid, &               !File id
         gsto3_inst_id, &                       !field_id for netCDF file (should match id set in def_var) 
         gsto3_inst, &                          !Tracer field (REAL4)
         start=srt_lon_lat_nlcat_time, &        !starting point for writing
         count=cnt_lon_lat_nlcat_time )         !Counts how many bytes written
    if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting GstO3_inst data')
    status = nf90_put_var(ncid, &               !File id
         fsto3_inst_id, &                       !field_id for netCDF file (should match id set in def_var) 
         fsto3_inst, &                          !Tracer field (REAL4)
         start=srt_lon_lat_nlcat_time, &        !starting point for writing
         count=cnt_lon_lat_nlcat_time )         !Counts how many bytes written
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting FstO3_inst data')

    !// close netcdf file
    status = nf90_close(ncid)
    if (status/=nf90_noerr) call handle_error(status,'close file')

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    

    !// --------------------------------------------------------------------
  end subroutine scav_diag_nmet_output_nc
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  subroutine scav_diag_2fileA(YEAR,MONTH,DAY,NDAY)
    !// --------------------------------------------------------------------
    !// Called from pmain AFTER calendar update.
    !//
    !// Write daily scavenged tracer masses of large scale and convective
    !// scavenging, and also daily average tracer burden.
    !//
    !// The SCAV_DIAG is accumulated data for the period defined
    !// by the BUDGETS calendar, but in the output it is
    !// converted to daily totals (this was not in the first version
    !// of these files).
    !//
    !// Amund Sovde Haslerud, June 2017
    !//   Output as netCDF4.
    !// Ole Amund Sovde, October 2014, June 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: IPARW, JPARW, LPARW, TNAMELEN
    use cmn_chem, only: TNAME, TMASS
    use cmn_met, only: METTYPE, metCYCLE, metREVNR, MET_ROOT, MPATH1
    use cmn_diag, only: JDO_T, RUNTITLE, metTypeInfo, resolutionInfo, &
         nc4deflate_global, nc4shuffle_global
    use cmn_ctm, only: LDLYSCAV
    use utilities, only: get_free_fileid
    use cmn_oslo, only: chem_idx, DINM
    use netcdf
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: YEAR, MONTH, DAY, NDAY
    !// Locals
    integer :: I, D, PD, TD, DYR, YDAY
    character(len=80) :: filename
    character(len=4) :: CYEAR
    real(r8), dimension(NPAR,4,366) :: SCAV_OUT
    integer, dimension(366) :: dayOfYear
    !//---------------------------------------------------------------------
    integer, parameter :: version = 3
    !//---------------------------------------------------------------------
    !// netCDF variables
    integer :: &
         status, &                 !Status for netcdf file 0=OK
         ncid, &                   !file id for output netcdf file
         ncomps_dim_id, &          !Dimension id for NPAR
         nlcat_dim_id, &           !Dimension id for NLCAT
         days_dim_id, &            !Dimension id for days
         days_id, &                !Dimension id for days of year (1:366)
         tracer_name_len_dim_id, & !Dimension id for tname charater length
         tracer_idx_id, &          !ID for all tracer number
         nlcat_id, &               !ID for land use categories
         dim_ncomps_days_id(2), &
         version_id, &             !ID for file version number
         tracer_molw_id, &         !ID for all tracer molecular weights
         tracer_name_id, &         !ID for all tracer names
         native_lon_id, &          !Variable id for native longitude size
         native_lat_id, &          !Variable id for native latitude size
         native_lev_id, &          !Variable id for native level size
         lon_id, &                 !Variable id for longitude
         lat_id, &                 !Variable id for latitude
         lev_id, &                 !Variable id for longitude interstices
         scav_ls_id, &
         scav_cnv_id, &
         scav_dry_id, & 
         scav_brd_id
         
    character(len=TNAMELEN), dimension(NPAR) :: tracer_name
    real(r8), dimension(NPAR) :: tracer_molw
    integer, dimension(NPAR) :: tracer_idx
    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'scav_diag_2fileA'
    !//---------------------------------------------------------------------

    !// Will only put calculate from 1Jan (or model start date if it differ),
    !// until the day before DAY. If DAY==1, a new year is coming up, and
    !// we use the number of days in the previous year
    DYR = sum(DINM)

    !// YEAR may have been updated by CALENDR before calling
    !// so we check if the new day is 1 Jan
    if (MONTH.eq.1 .and. DAY.eq.1) then
       write(CYEAR(1:4),'(i4.4)') YEAR-1
       YDAY = DYR
    else
       write(CYEAR(1:4),'(i4.4)') YEAR
       YDAY = DAY - 1
    end if


    !// Change from accumulated data to daily totals
    !// Loop through calendar to find output dates
    !// When JDO_T>0, the accumulated data are zeroed AFTER
    !// written to file.
    SCAV_OUT(:,:,:) = 0._r8        !// Will be calculated
    do D = 1, YDAY
       !// Previous day
       if (D .eq. 1) then
          PD = DYR
       else
          PD = D - 1
       end if
       !// Day number for JDO_T
       if (D .gt. 59 .and. DINM(2) .eq. 28._r8) then
          !// Not leap year, so we have to jump past 60 in JDO_T
          TD = D+1
       else
          TD = D
       end if
       !// Tendencies calendar tells us when to subtract value from
       !// previous day.
       if (JDO_T(TD) .eq. 0) then
          !// Subtract values from day before
          SCAV_OUT(:,:,D) = SCAV_DIAG(:,:,D) - SCAV_DIAG(:,:,PD)
       else
          !// Data is only accumulated during the current day
          SCAV_OUT(:,:,D) = SCAV_DIAG(:,:,D)
       end if
    end do

    !// Old write to file
    !filename = 'scavenging_daily_totals_'//CYEAR//'.dta'
    !open(fnr,file=filename,form='unformatted')
    !write(fnr) NPAR, YEAR, version
    !write(fnr) TNAME
    !write(fnr) chem_idx
    !write(fnr) TMASS
    !write(fnr) SCAV_OUT
    !close(fnr)

    !//---------------------------------------------------------------------
    !// Write netCDF4 file
    filename = 'scavenging_daily_totals_'//CYEAR//'.nc'
    !// Open netCDF4 file for writing
    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')
    !//---------------------------------------------------------------------

    !//---------------------------------------------------------------------
    !//File headers
    status=nf90_put_att(ncid,nf90_global, 'title', &
         'Oslo CTM3 total daily scavenging')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': file header')
    status=nf90_put_att(ncid,nf90_global,'model_info', &
         'Oslo CTM3 is a 3D Chemical Transport Model')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': model_info')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology',metTypeInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology_path',MPATH1)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology_path')
    status=nf90_put_att(ncid,nf90_global,'resolution_info',resolutionInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': resolution_info')
    status=nf90_put_att(ncid,nf90_global,'runtitle',trim(runtitle))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': runtitle')
    status=nf90_put_att(ncid,nf90_global,'contact_info','CICERO')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': contact_info')

    !// Define sizes

    !// Define NCOMPS = NPAR
    status = nf90_def_dim(ncid,"NCOMPS",NPAR,ncomps_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NCOMPS dim')
    !// All the component IDs
    status = nf90_def_var(ncid,"tracer_idx",nf90_int,ncomps_dim_id,tracer_idx_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_idx variable')
    status = nf90_put_att(ncid,tracer_idx_id,'description', &
         'ID numbers for species used in Oslo CTM3.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_idx')

    !// Define length of tracer name string (use TNAME for this)
    status = nf90_def_dim(ncid,"tracer_name_len",TNAMELEN,tracer_name_len_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name_len dim')

        !// Defining the DAYS variable
    status = nf90_def_dim(ncid,"DAYS",366,days_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define DAYS dim')
    status = nf90_def_var(ncid,"DAYS",nf90_int,days_dim_id,days_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define DAYS variable')
    !// Attributes
    status = nf90_put_att(ncid,days_id,'units','Day of year (1:366)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units DAYS')


    !// Info about native size IPARW
    status = nf90_def_var(ncid,"IPARW",nf90_int,native_lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define IPARW variable')
    status = nf90_put_att(ncid,native_lon_id,'description', &
         'Meteorological data native longitudinal resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description IPARW')
    !// Info about native size JPARW
    status = nf90_def_var(ncid,"JPARW",nf90_int,native_lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define JPARW variable')
    status = nf90_put_att(ncid,native_lat_id,'description', &
         'Meteorological data native latitudinal resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description JPARW')
    !// Info about native size LPARW
    status = nf90_def_var(ncid,"LPARW",nf90_int,native_lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define LPARW variable')
    status = nf90_put_att(ncid,native_lev_id,'description', &
         'Meteorological data vertical resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description LPARW')

    !// Info about size IPAR
    status = nf90_def_var(ncid,"IPAR",nf90_int,lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define IPAR variable')
    status = nf90_put_att(ncid,lon_id,'description', &
         'Longitudinal resolution of simulation.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description IPAR')
    !// Info about native size JPAR
    status = nf90_def_var(ncid,"JPAR",nf90_int,lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define JPAR variable')
    status = nf90_put_att(ncid,lat_id,'description', &
         'Latitudinal resolution of simulation.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description JPARW')
    !// Info about native size LPAR
    status = nf90_def_var(ncid,"LPAR",nf90_int,lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define LPAR variable')
    status = nf90_put_att(ncid,lev_id,'description', &
         'Vertical resolution of simulation.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description LPAR')

    !// File version number
    status = nf90_def_var(ncid,"VERSION",nf90_int,version_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define VERSION variable')
    status = nf90_put_att(ncid,version_id,'description', &
         'Output file version number')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description VERSION')

    !// Molecular masses of transported components (r8)
    status = nf90_def_var(ncid,"tracer_molweight",nf90_double,ncomps_dim_id,tracer_molw_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_molweight variable')
    status = nf90_put_att(ncid,tracer_molw_id,'units','g/mol')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units tracer_molweight')
    status = nf90_put_att(ncid,tracer_molw_id,'description', &
         'Molecular weights.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_molweight')

    !// Tracer names
    status = nf90_def_var(ncid,"tracer_name",nf90_char, &
         (/tracer_name_len_dim_id,ncomps_dim_id/),tracer_name_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define trace_name variable')
    status = nf90_put_att(ncid,tracer_name_id,'description', &
         'Tracer/species names.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_name')

    !// Defining the combined id for a 2d field (lon,lat)
    dim_ncomps_days_id(1) = ncomps_dim_id
    dim_ncomps_days_id(2) = days_dim_id

    !// LS, CNV, DRY and BRD
    !// Totals wet scavenged, large scale (r8), deflate netcdf4
    if (LDLYSCAV(4)) then
       status = nf90_def_var(ncid,"SCAV_LS",nf90_double,dim_ncomps_days_id,scav_ls_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define SCAV_LS variable')
       status = nf90_def_var_deflate(ncid,scav_ls_id,nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate SCAV_LS variable')
       status = nf90_put_att(ncid,scav_ls_id,'units','kg(tracer)')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units SCAV_LS')
       status = nf90_put_att(ncid,scav_ls_id,'description', &
            'Daily total mass scavenged by large scale scavenging.')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description SCAV_LS')
    end if
    !// Totals wet scavenged, convective (r8), deflate netcdf4
    if (LDLYSCAV(5)) then
       status = nf90_def_var(ncid,"SCAV_CNV",nf90_double,dim_ncomps_days_id,scav_cnv_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define SCAV_CNV variable')
       status = nf90_def_var_deflate(ncid,scav_cnv_id,nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate SCAV_CNV variable')
       status = nf90_put_att(ncid,scav_cnv_id,'units','kg(tracer)')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units SCAV_CNV')
       status = nf90_put_att(ncid,scav_cnv_id,'description', &
            'Daily total mass scavenged by convectivescavenging.')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description SCAV_CNV')
    end if
    if (LDLYSCAV(6)) then
       !// Totals wet scavenged, convective (r8), deflate netcdf4
       status = nf90_def_var(ncid,"SCAV_DRY",nf90_double,dim_ncomps_days_id,scav_dry_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define SCAV_DRY variable')
       status = nf90_def_var_deflate(ncid,scav_dry_id,nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate SCAV_DRY variable')
       status = nf90_put_att(ncid,scav_dry_id,'units','kg(tracer)')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units SCAV_DRY')
       status = nf90_put_att(ncid,scav_dry_id,'description', &
            'Daily total mass scavenged by dry deposition.')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description SCAV_DRY')
    end if
    !// Average burdens (r8), deflate netcdf4
    if (LDLYSCAV(3)) then
       status = nf90_def_var(ncid,"SCAV_BRD",nf90_double,dim_ncomps_days_id,scav_brd_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define SCAV_BRD variable')
       status = nf90_def_var_deflate(ncid,scav_brd_id,nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate SCAV_BRD variable')
       status = nf90_put_att(ncid,scav_brd_id,'units','kg(tracer)')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units SCAV_BRD')
       status = nf90_put_att(ncid,scav_brd_id,'description', &
            'Daily average total mass (burden) of tracer.')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description SCAV_BRD')
    end if
    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------

    !// Putting the day of year
    do I = 1, 366
       dayOfYear(I) = I
    end do
    status = nf90_put_var(ncid,days_id,dayOfYear)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting days')

    !// Info about native size
    status = nf90_put_var(ncid,native_lon_id,IPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting IPARW')
    status = nf90_put_var(ncid,native_lat_id,JPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting JPARW')
    status = nf90_put_var(ncid,native_lev_id,LPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting LPARW')

    !// Info about model run size
    status = nf90_put_var(ncid,lon_id,IPAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting IPAR')
    status = nf90_put_var(ncid,lat_id,JPAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting JPAR')
    status = nf90_put_var(ncid,lev_id,LPAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting LPAR')

    !// File version number
    status = nf90_put_var(ncid,version_id,version)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting VERSION')

    !// All the transported components
    status = nf90_put_var(ncid,tracer_idx_id,chem_idx)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//':putting tracer_idx')

    !// Molecular masses of transported components (r8)
    status = nf90_put_var(ncid,tracer_molw_id,TMASS)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_molw')

    !// Tracer names
    status = nf90_put_var(ncid,tracer_name_id,TNAME)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_name')

    !// Scavenged LS
    if (LDLYSCAV(4)) then
       status = nf90_put_var(ncid,scav_ls_id,SCAV_OUT(:,1,:))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SCAV_LS')
    end if
    !// Scavenged CNV
    if (LDLYSCAV(5)) then
       status = nf90_put_var(ncid,scav_cnv_id,SCAV_OUT(:,2,:))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SCAV_CNV')
    end if
    !// Scavenged DRY
    if (LDLYSCAV(6)) then
       status = nf90_put_var(ncid,scav_dry_id,SCAV_OUT(:,3,:))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SCAV_DRY')
    end if
    !// Burden average BRD
    if (LDLYSCAV(3)) then
       status = nf90_put_var(ncid,scav_brd_id,SCAV_OUT(:,4,:))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SCAV_BRD')
    end if
    
    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) &
         call handle_error(status,'close file: '//trim(filename))
    !//---------------------------------------------------------------------


    !// Reset arrays: This is only done when JDO_T is true for the
    !// upcoming day.
    if (D .gt. 59 .and. DINM(2) .eq. 28._r8) then
       !// Not leap year, so we have to jump past 60 in JDO_T
       TD = DAY+1
    else
       TD = DAY
    end if
    !// Just in case (when called from correct place in pmain,
    !// this is not necessary)
    if (JDO_T(TD) .eq. 0) return

    SCAV_LS(:,:) = 0._r8
    SCAV_CN(:,:) = 0._r8
    SCAV_DD(:,:) = 0._r8
    SCAV_BRD(:,:) = 0._r8
    

    !// I think it is correct that SCAV_DIAG is only initialised
    !// at model start, and not at the beginning of each year.

    !// --------------------------------------------------------------------
  end subroutine scav_diag_2fileA
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine scav_diag_2fileB(YEAR,MONTH,DAY,NDAY)
    !// --------------------------------------------------------------------
    !// Write 2D daily scavenged tracer masses of large scale and
    !// convective scavenging, plus dry deposition.
    !// Written as double precision.
    !//
    !// Amund Sovde Haslerud, June 2017
    !//   Output as netCDF4.
    !// Ole Amund Sovde, September 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: IPARW, JPARW, LPARW, TNAMELEN
    use cmn_ctm, only: NROPSM, NRMETD, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         XDGRD, YDGRD, XDEDG, YDEDG, AREAXY, LDLYSCAV
    use cmn_chem, only: TNAME, TMASS
    use cmn_met, only: METTYPE, metCYCLE, metREVNR, MET_ROOT, MPATH1
    use cmn_diag, only: RUNTITLE, metTypeInfo, resolutionInfo, &
         nc4deflate_global, nc4shuffle_global
    use utilities, only: get_free_fileid
    use cmn_oslo, only: chem_idx
    use cmn_sfc, only: LDDEPmOSaic
    use netcdf
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: YEAR, MONTH, DAY, NDAY
    !// Locals
    character(len=80) :: filename
    character(len=4) :: CYEAR, CDATE
    real(r8) :: RDUM(IPAR,JPAR,NPAR)
    integer :: I,J,II,JJ,N,MP,K
    !// --------------------------------------------------------------------
    integer, parameter :: version = 3
    !//---------------------------------------------------------------------
    !// netCDF variables
    !// Dimensions ID and counters
    integer :: lon_dim_id             !Dimension ID for longitude
    integer :: lat_dim_id             !Dimension ID for latitude
    integer :: ilon_dim_id            !Dimension ID for longitude interstices
    integer :: ilat_dim_id            !Dimension ID for latitude interstices
    integer :: ncomps_dim_id          !Dimension ID for NPAR
    integer :: nlcat_dim_id           !Dimension ID for NLCAT
    integer :: tracer_name_len_dim_id !Dimension ID for tname charater length
    integer :: lon_id                 !Variable ID for longitude
    integer :: lat_id                 !Variable ID for latitude
    integer :: ilon_id                !Variable ID for longitude interstices
    integer :: ilat_id                !Variable ID for latitude interstices
    integer :: nlcat_id               !Variable ID for nlcat
    integer :: tracer_idx_id          !ID for all tracer number
    integer :: dim_lon_lat_id(2)      !IF for lon/lat
    integer :: dim_lon_lat_nlcat_id(3)!ID for lon/lat/nlcat
    integer :: version_id             !ID for file version number
    integer :: tracer_molw_id         !ID for all tracer molecular weights
    integer :: tracer_name_id         !ID for all tracer names
    integer :: native_lon_id          !Variable ID for native longitude size
    integer :: native_lat_id          !Variable ID for native latitude size
    integer :: native_lev_id          !Variable ID for native level size
    integer :: current_lon_id         !Variable ID for simulation longitude size
    integer :: current_lat_id         !Variable ID for simulation latitude size
    integer :: current_lev_id         !Variable ID for simulation level size
    integer :: areaxy_id              !Variable ID for grid area
    integer :: date_year_id           !Variable ID for year
    integer :: date_month_id          !Variable ID for month
    integer :: date_day_id            !Variable ID for day
    integer :: scav_ls_id(NPAR), &
         scav_cnv_id(NPAR), &
         scav_dry_id(NPAR)
    integer :: &
         gsto3_avg_id, &
         vrao3_avg_id, &
         vrbo3_avg_id, &
         vrco3_avg_id, &
         fsto3_avg_id

    !// Other locals
    integer :: status                 !Status for netcdf file 0=OK
    integer :: ncid                   !file id for output netcdf file
    integer, dimension(NLCAT) :: LCAT =  &   !Landuse category numbering
         (/(K, K=1,NLCAT, 1)/) 

    character(len=TNAMELEN), dimension(NPAR) :: tracer_name
    real(r8), dimension(NPAR) :: tracer_molw
    integer, dimension(NPAR) :: tracer_idx
    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'scav_diag_2fileB'
    !//---------------------------------------------------------------------

    !// YEAR/MONTH/DAY have not yet been updated by calendar, hence
    !// applies for current day.
    write(CYEAR(1:4),'(i4.4)') YEAR
    write(CDATE(1:4),'(i2.2,i2.2)') MONTH,DAY

    !//---------------------------------------------------------------------
    !// Write netCDF4 file
    filename = 'scavenging_daily_2d_'//CYEAR//CDATE//'.nc'
    !// Open netCDF4 file for writing
    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')
    !//---------------------------------------------------------------------

    !//---------------------------------------------------------------------
    !//File headers
    status=nf90_put_att(ncid,nf90_global,'title','Oslo CTM3 daily scavenging maps')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': file header')
    status=nf90_put_att(ncid,nf90_global,'model_info', &
         'Oslo CTM3 is a 3D Chemical Transport Model')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': model_info')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology',metTypeInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology_path',MPATH1)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology_path')
    status=nf90_put_att(ncid,nf90_global,'resolution_info',resolutionInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': resolution_info')
    status=nf90_put_att(ncid,nf90_global,'runtitle',trim(runtitle))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': runtitle')
    status=nf90_put_att(ncid,nf90_global,'contact_info','CICERO')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': contact_info')

    !// Define sizes
    !// Define lat/lon
    status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat dim')

    status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon dim')

    !// Defining nlcat
    status = nf90_def_dim(ncid,"NLCAT",NLCAT,nlcat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NLCAT dim')
    
    !// Define NCOMPS = NPAR
    status = nf90_def_dim(ncid,"NCOMPS",NPAR,ncomps_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NCOMPS dim')

    !// Define ilat/ilon
    status = nf90_def_dim(ncid,"ilat",JPAR+1,ilat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat dim')
    status = nf90_def_dim(ncid,"ilon",IPAR+1,ilon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//':define ilon dim')

    !// Define length of tracer name string (use TNAME for this)
    status = nf90_def_dim(ncid,"tracer_name_len",TNAMELEN,tracer_name_len_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name_len dim')

    !// All the component IDs
    status = nf90_def_var(ncid,"tracer_idx",nf90_int,ncomps_dim_id,tracer_idx_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_idx variable')
    status = nf90_put_att(ncid,tracer_idx_id,'description', &
         'Component ID in Oslo CTM3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_idx')

    !// Info about YEAR
    status = nf90_def_var(ncid,"YEAR",nf90_int,date_year_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define YEAR variable')
    status = nf90_put_att(ncid,date_year_id,'description', &
         'Year of accumulated scavenged data')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description Year')
    !// Info about MONTH
    status = nf90_def_var(ncid,"MONTH",nf90_int,date_month_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define MONTH variable')
    status = nf90_put_att(ncid,date_month_id,'description', &
         'Month of accumulated scavenged data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description MONTH')
    !// Info about DAY
    status = nf90_def_var(ncid,"DAY",nf90_int,date_day_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define DAY variable')
    status = nf90_put_att(ncid,date_day_id,'description', &
         'Day in month of the accumulated scavenged data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description DAY')



    !// Info about native size IPARW
    status = nf90_def_var(ncid,"IPARW",nf90_int,native_lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define IPARW variable')
    status = nf90_put_att(ncid,native_lon_id,'description', &
         'Meteorological data native longitudinal resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description IPARW')
    !// Info about native size JPARW
    status = nf90_def_var(ncid,"JPARW",nf90_int,native_lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define JPARW variable')
    status = nf90_put_att(ncid,native_lat_id,'description', &
         'Meteorological data native latitudinal resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description JPARW')
    !// Info about native size LPARW
    status = nf90_def_var(ncid,"LPARW",nf90_int,native_lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define LPARW variable')
    status = nf90_put_att(ncid,native_lev_id,'description', &
         'Meteorological data native vertical resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description LPARW')

    !// Info about size IPAR
    status = nf90_def_var(ncid,"IPAR",nf90_int,current_lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define IPAR variable')
    status = nf90_put_att(ncid,current_lon_id,'description', &
         'Longitudinal resolution of simulation.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description IPAR')
    !// Info about native size JPAR
    status = nf90_def_var(ncid,"JPAR",nf90_int,current_lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define JPAR variable')
    status = nf90_put_att(ncid,current_lat_id,'description', &
         'Latitudinal resolution of simulation.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description JPAR')
    !// Info about native size LPAR
    status = nf90_def_var(ncid,"LPAR",nf90_int,current_lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define LPAR variable')
    status = nf90_put_att(ncid,current_lev_id,'description', &
         'Vertical resolution of simulation.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description LPAR')

    !// File version number
    status = nf90_def_var(ncid,"VERSION",nf90_int,version_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define VERSION variable')
    status = nf90_put_att(ncid,version_id,'description', &
         'Output file version number')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description VERSION')

    !// Molecular masses of transported components (r8)
    status = nf90_def_var(ncid,"tracer_molweight",nf90_double,ncomps_dim_id,tracer_molw_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_molweight variable')
    status = nf90_put_att(ncid,tracer_molw_id,'units','g/mol')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units tracer_molweight')
    status = nf90_put_att(ncid,tracer_molw_id,'description', &
         'Molecular weights.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_molweight')

    !// Tracer names
    status = nf90_def_var(ncid,"tracer_name",nf90_char, &
         (/tracer_name_len_dim_id,ncomps_dim_id/),tracer_name_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define trace_name variable')
    status = nf90_put_att(ncid,tracer_name_id,'description', &
         'Tracer/species names.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_name')

    !// Defining the combined id for a 2d field (lon,lat)
    dim_lon_lat_id(1) = lon_dim_id
    dim_lon_lat_id(2) = lat_dim_id
    !// Defining the combined id for a 3d field (lon,lat,nlcat)
    dim_lon_lat_nlcat_id(1) = lon_dim_id
    dim_lon_lat_nlcat_id(2) = lat_dim_id
    dim_lon_lat_nlcat_id(3) = nlcat_dim_id
   

    !// Defining the lon variable
    status = nf90_def_var(ncid,"lon",nf90_double,lon_dim_id,lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon variable')
    !// Attributes
    status = nf90_put_att(ncid,lon_id,'units','degrees east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lon')
    status = nf90_put_att(ncid,lon_id,'description','Value at grid box center.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description lon')

    !// Defining the lat variable
    status = nf90_def_var(ncid,"lat",nf90_double,lat_dim_id,lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat variable')
    !// Attributes
    status = nf90_put_att(ncid,lat_id,'units','degrees north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lat')
    status = nf90_put_att(ncid,lat_id,'description','Value at grid box center.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description lat')

    !// Defining ilon variable (lon on interstices)
    status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    !// Attributes
    status = nf90_put_att(ncid,ilon_id,'units','degrees east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilon')
    status = nf90_put_att(ncid,ilon_id,'description','Value at eastern edge of grid box.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ilon')

    !// Defining ilat variable (lat on interstices)
    status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')
    !// Attributes
    status = nf90_put_att(ncid,ilat_id,'units','degrees north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilat')
    status = nf90_put_att(ncid,ilat_id,'description','Value at southern edge of grid box.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ilat')

    !// Defining the nlcat-variable
    status = nf90_def_var(ncid,"NLCAT",nf90_int,nlcat_dim_id,nlcat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define nlcat variable')
    status = nf90_put_att(ncid,nlcat_id,'description', &
         'Landuse categories used in Oslo CTM3.' // char(10) // &
            '01 - Needleleaftree temp./bor.(CF)' // char(10) // &
            '02 - Deciduoustree temp./bor. (DF)' // char(10) // &
            '03 - Needleleaftree med. (NF)' // char(10) // &
            '04 - Broadleaftree (BF)' // char(10) // &
            '05 - Crops (TC)' // char(10) // &
            '06 - Moorland (SNL)' // char(10) // &
            '07 - Grassland (GR)' // char(10) // &
            '08 - Scrubs med. (MS)' // char(10) // &
            '09 - Wetlands (WE)' // char(10) // &
            '10 - Tundra (TU)' // char(10) // &
            '11 - Desert (DE)' // char(10) // &
            '12 - Water (W)' // char(10) // &
            '13 - Urban (U)' // char(10) // &
            '14 - Ice/Snow (ICE)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description nlcat')       

    !// Grid area (r8), deflate netcdf4
    status = nf90_def_var(ncid,"gridarea",nf90_double,dim_lon_lat_id,areaxy_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define gridarea variable')
    status = nf90_def_var_deflate(ncid,areaxy_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate gridarea variable')
    status = nf90_put_att(ncid,areaxy_id,'unit','m2')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit gridarea')



    !// Each component
    !// Name convention? ls_O3/cnv_O3/dry_O3?
    !// SCAV_LS
    if (LDLYSCAV(4)) then
       do N = 1, NPAR
          status = nf90_def_var(ncid,'ls_'//trim(TNAME(N)), &
               nf90_double, dim_lon_lat_id, scav_ls_id(N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define SCAV_LS variable '//trim(TNAME(N)))
          status = nf90_def_var_deflate(ncid, scav_ls_id(N), &
               nc4shuffle, 1, nc4deflate)
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define deflate SCAV_LS variable '//trim(TNAME(N)))
          status = nf90_put_att(ncid,scav_ls_id(N),'unit','kg(tracer)')
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': attribute unit SCAV_LS '//trim(TNAME(N)))
       end do
    end if
    
    !// SCAV_CNV
    if (LDLYSCAV(5)) then
       do N = 1, NPAR
          status = nf90_def_var(ncid,'cnv_'//trim(tname(N)), &
               nf90_double, dim_lon_lat_id, scav_cnv_id(N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define SCAV_CNV variable '//trim(TNAME(N)))
          status = nf90_def_var_deflate(ncid, scav_cnv_id(N), &
               nc4shuffle, 1, nc4deflate)
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define deflate SCAV_CNV variable '//trim(TNAME(N)))
          status = nf90_put_att(ncid, scav_cnv_id(N), 'unit', 'kg(tracer)')
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': attribute unit SCAV_CNV '//trim(TNAME(N)))
       end do
    end if

    !// SCAV_DRY
    if (LDLYSCAV(6)) then
       do N = 1, NPAR
          status = nf90_def_var(ncid,'dry_'//trim(tname(N)), &
               nf90_double, dim_lon_lat_id, scav_dry_id(N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define SCAV_DRY variable '//trim(TNAME(N)))
          status = nf90_def_var_deflate(ncid, scav_dry_id(N), &
               nc4shuffle, 1, nc4deflate)
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define deflate SCAV_DRY variable '//trim(TNAME(N)))
          status = nf90_put_att(ncid,scav_dry_id(N),'unit','kg(tracer)')
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': attribute unit SCAV_DRY '//trim(TNAME(N)))
       end do
    end if

    if(LDDEPmOSaic .and. LDLYSCAV(7)) then
       !// Daily average of stomatal conductance
       status = nf90_def_var(ncid,'GstO3_avg', &
            nf90_double, dim_lon_lat_nlcat_id, gsto3_avg_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable GstO3_avg ')
       status = nf90_def_var_deflate(ncid, gsto3_avg_id, &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable GstO3_avg ')
       status = nf90_put_att(ncid, gsto3_avg_id,'unit','mmol m-2 s-1')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute unit GstO3_avg ')

       !// Daily average of aerodynamical dry deposition velocity
       status = nf90_def_var(ncid,'VRaO3_avg', &
            nf90_double, dim_lon_lat_nlcat_id, vrao3_avg_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable VRaO3_avg ')
       status = nf90_def_var_deflate(ncid, vrao3_avg_id, &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable VRaO3_avg ')
       status = nf90_put_att(ncid, vrao3_avg_id,'unit','m s-1')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute unit VRaO3_avg ')

       !// Daily average of quasi-laminar dry deposition velocity
       status = nf90_def_var(ncid,'VRbO3_avg', &
            nf90_double, dim_lon_lat_nlcat_id, vrbo3_avg_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable VRbO3_avg ')
       status = nf90_def_var_deflate(ncid, vrbo3_avg_id, &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable VRbO3_avg ')
       status = nf90_put_att(ncid, vrbo3_avg_id,'unit','m s-1')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute unit VRbO3_avg ')

       !// Daily average of canopy dry deposition velocity
       status = nf90_def_var(ncid,'VRcO3_avg', &
            nf90_double, dim_lon_lat_nlcat_id, vrco3_avg_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable VRcO3_avg ')
       status = nf90_def_var_deflate(ncid, vrco3_avg_id, &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable VRcO3_avg ')
       status = nf90_put_att(ncid, vrco3_avg_id,'unit','m s-1')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute unit VRcO3_avg ')

       !// Daily average of stomatal flux
       status = nf90_def_var(ncid,'FstO3_avg', &
            nf90_double, dim_lon_lat_nlcat_id, fsto3_avg_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable FstO3_avg ')
       status = nf90_def_var_deflate(ncid, fsto3_avg_id, &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable FstO3_avg ')
       status = nf90_put_att(ncid, fsto3_avg_id,'unit','nmol m-2 s-1')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute unit FstO3_avg ')
    end if
    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------

    !// Putting the lon/lat/lev variables
    status = nf90_put_var(ncid,lon_id,XDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lon')
    status = nf90_put_var(ncid,lat_id,YDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lat')

    status = nf90_put_var(ncid,ilon_id,XDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilon')
    status = nf90_put_var(ncid,ilat_id,YDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilat')

    !// Putting the NLCAT variables
    status = nf90_put_var(ncid,nlcat_id,LCAT)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting nlcat')

    !// Info about date
    status = nf90_put_var(ncid,date_year_id,YEAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting YEAR')
    status = nf90_put_var(ncid,date_month_id,MONTH)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting MONTH')
    status = nf90_put_var(ncid,date_day_id,DAY)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting DAY')

    !// Info about native size
    status = nf90_put_var(ncid,native_lon_id,IPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting IPARW')
    status = nf90_put_var(ncid,native_lat_id,JPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting JPARW')
    status = nf90_put_var(ncid,native_lev_id,LPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting LPARW')

    !// Info about model run size
    status = nf90_put_var(ncid,current_lon_id,IPAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting IPAR')
    status = nf90_put_var(ncid,current_lat_id,JPAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting JPAR')
    status = nf90_put_var(ncid,current_lev_id,LPAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting LPAR')

    !// File version number
    status = nf90_put_var(ncid,version_id,version)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting VERSION')

    !// All the transported components
    status = nf90_put_var(ncid,tracer_idx_id,chem_idx)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//':putting tracer_idx')

    !// The land use categories
    status = nf90_put_var(ncid,nlcat_id,LCAT)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//':putting nlcat_id')

    !// Molecular masses of transported components (r8)
    status = nf90_put_var(ncid,tracer_molw_id,TMASS)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_molw')

    !// Tracer names
    status = nf90_put_var(ncid,tracer_name_id,TNAME)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_name')

    !// Grid area (r8)
    status = nf90_put_var(ncid,areaxy_id,AREAXY)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting AREAXY')

    !// Scavenged LS: SCAV_MAP_WLS
    if (LDLYSCAV(4)) then
       do MP = 1, MPBLK
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do N = 1, NPAR
                   RDUM(I,J,N) = SCAV_MAP_WLS(N,II,JJ,MP)
                end do
             end do
          end do
       end do
       
       do N = 1, NPAR
          status = nf90_put_var(ncid,scav_ls_id(N),RDUM(:,:,N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': putting SCAV_MAP_LS')
       end do
    end if

    !// Scavenged CNV: SCAV_MAP_WCN
    if (LDLYSCAV(5)) then
       do MP = 1, MPBLK
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do N = 1, NPAR
                   RDUM(I,J,N) = SCAV_MAP_WCN(N,II,JJ,MP)
                end do
             end do
          end do
       end do
       
       do N = 1, NPAR
          status = nf90_put_var(ncid,scav_cnv_id(N),RDUM(:,:,N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': putting SCAV_MAP_CNV')
       end do
    end if
    
    !// Scavenged DRY: SCAV_MAP_DRY
    if (LDLYSCAV(6)) then
       do MP = 1, MPBLK
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do N = 1, NPAR
                   RDUM(I,J,N) = SCAV_MAP_DRY(N,II,JJ,MP)
                end do
             end do
          end do
       end do
       
       do N = 1, NPAR
          status = nf90_put_var(ncid,scav_dry_id(N),RDUM(:,:,N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': putting SCAV_MAP_DRY')
       end do
    end if

    !// Dry deposition velocities
    if(LDDEPmOSaic .and. LDLYSCAV(7)) then
       status = nf90_put_var(ncid,gsto3_avg_id,GSTO3_AVG)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting GstO3_avg')
       status = nf90_put_var(ncid,vrao3_avg_id,VRAO3_AVG)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting VRaO3_avg')
       status = nf90_put_var(ncid,vrbo3_avg_id,VRBO3_AVG)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting VRbO3_avg')
       status = nf90_put_var(ncid,vrco3_avg_id,VRCO3_AVG)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting VRcO3_avg')
       !// [mmol m-2 s-1] -> [nmol m-2 s-1]
       status = nf90_put_var(ncid,fsto3_avg_id,FSTO3_AVG*1.e6_r8)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting FstO3_avg')
    end if
    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) &
         call handle_error(status,'close file: '//trim(filename))
    !//---------------------------------------------------------------------


    !// Reset
    SCAV_MAP_WLS(:,:,:,:) = 0._r8
    SCAV_MAP_WCN(:,:,:,:) = 0._r8
    SCAV_MAP_DRY(:,:,:,:) = 0._r8
    GSTO3_AVG(:,:,:) = 0._r8
    VRAO3_AVG(:,:,:) = 0._r8
    VRBO3_AVG(:,:,:) = 0._r8
    VRCO3_AVG(:,:,:) = 0._r8
    FSTO3_AVG(:,:,:) = 0._r8
    !// --------------------------------------------------------------------
  end subroutine scav_diag_2fileB
  !// ----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
end module diagnostics_scavenging
!//=========================================================================
