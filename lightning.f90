!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Lightning emissions.
!//=========================================================================
module lightning
  !// ----------------------------------------------------------------------
  !// MODULE: lightning
  !// DESCRIPTION: Module conatining routines to calculate lightning NOx
  !//              emissions. Several options are available; from the first
  !//              CTM3 version (GMD2012) as described by Søvde et al. (2012,
  !//              GMD, doi:10.5194/gmd-7-175-2014) to the UCI2015 version
  !//              (not documented yet).
  !//
  !// I was not happy with the UCI2015 version, so I have combined
  !// parts of it with the GMD2012 version; it produces better
  !// annual cycle in total flash rates, and also acceptable
  !// horizontal distribution.
  !//
  !// The lightning subroutine calculates the 3-D lightning source
  !// for use in the SOURCE subroutine or as source term in Oslo
  !// chemistry. It must be called for each met field after metdata
  !// is updated.
  !//
  !// Contains:
  !//   - subroutine LIGHTNING_OAS2015
  !//   - subroutine LIGHTNING_GMD2012 (LIGHTNING_CDH)
  !//   - subroutine LIGHTNING_UCI2015
  !//   - subroutine LIGHTNING_ALLEN2002
  !//   - subroutine LIGHTDIST
  !//
  !// Ole Amund Sovde, March 2015, November 2014, August 2011
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8, r4
  use cmn_parameters, only: secYear
  use cmn_chem, only: INFILE_LIGHTNING !// Ott etal data
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

   integer, parameter :: NNLIGHT=3200
  integer, parameter :: NLTYPE=4
  real(r8), dimension(NNLIGHT,NLTYPE) :: LITPROFILE

  !// Indices to indicate latitudes to regard as polar
  integer :: JPS, JPN

  !// Diagnose amount of lightning emissions
  real(r8) :: totlitn, totlitfrq

  !// Scaling factors.
  !// The total accumulated modelled flash rates for ocean
  !// and land, and the accumulated time, is needed to calculate
  !// the scaling factors. As will be described in the subroutines,
  !// the scaling factors convert from modeled flash rate to
  !// climatologically normalised fraction of observed annual
  !// flash rates.
  real(r8) :: accFocean, accFland, accDT

  !// Logical to flag if entrainment is included or not.
  logical :: LentuExist

  !// Seconds in a year
  real(r8), parameter :: dtyear = secYear
  !// Flag to turn on flashrate.dta production
  logical, parameter :: getFlashFiles = .false.
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'lightning.f90'
  !// ----------------------------------------------------------------------
  !// Save variables defined above
  save
  !// All is private
  private
  !// except the routines
  public LIGHTNING_OAS2015, LIGHTNING_UCI2015, LIGHTNING_GMD2012,  &
         LIGHTNING_ALLEN2002
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine getScaleFactors(scaleOcean, scaleLand)
    !// --------------------------------------------------------------------
    !// Scaling parameters depending on meteorological data.
    !//
    !// Price etal equations have units 1/minutes, while observations
    !// are given in 1/s. This conversion is taken care of by the factors
    !// scaleLand and scaleOcean. Thus, after multiplying with the
    !// scaling factors, the units of flash rate is fraction of
    !// all annual lightnings per second. Multiplying by time step
    !// and annual amount of emission (e.g. 5Tg(N)), gives the LNOx
    !// emissions for the time step.
    !// --------------------------------------------------------------------
    use cmn_size, only: IPARW, IDGRD
    use cmn_met, only: metTYPE, metCYCLE, metREVNR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(out) :: scaleOcean, scaleLand
    !// --------------------------------------------------------------------
    !// ECMWF IFS cy36r1
    !// T42L60: 1997 - 2010 cy36r1 average scaling factor (Ocean, Land)
    real(r8), parameter :: Oecmwf_ifs_c36r1_T42 = 2.304304e-13_r8
    real(r8), parameter :: Lecmwf_ifs_c36r1_T42 = 4.329157e-16_r8
    !// T159L60: (T159 2007 / T42 2007 * clim.T42)
    real(r8), parameter :: Oecmwf_ifs_c36r1_T159 = 3.594025e-14_r8
    real(r8), parameter :: Lecmwf_ifs_c36r1_T159 = 7.061138e-17_r8

    !// ECMWF oIFS cy38r1
    !// T159L60: 1995 - 2005
    real(r8), parameter :: Oecmwf_oifs_c38r1_T159 = 4.012695e-14_r8
    real(r8), parameter :: Lecmwf_oifs_c38r1_T159 = 9.542578e-17_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'getScaleFactors'
    !// --------------------------------------------------------------------

    !// Initialise
    scaleOcean = 0._r8
    scaleLand  = 0._r8

    if (trim(metTYPE) .eq. 'ECMWF_IFS') then
       !// ECMWF IFS
       !// -----------------------------------------------------------------
       if (metCYCLE .eq. 36 .and. metREVNR .eq. 1) then
          !// Cycle 36r1
          if (IPARW .eq. 128) then
             if (IDGRD .eq. 1) then
                scaleOcean = Oecmwf_ifs_c36r1_T42
                scaleLand  = Lecmwf_ifs_c36r1_T42
             else if (IDGRD .eq. 2) then
                !// Scaling factors based on cy36 2007 meteorology
                scaleOcean = Oecmwf_ifs_c36r1_T42 * 1.642020_r8
                scaleLand  = Lecmwf_ifs_c36r1_T42 * 1.484240_r8
             else if (IDGRD .eq. 4) then
                !// Scaling factors based on cy36 2007 meteorology
                scaleOcean = Oecmwf_ifs_c36r1_T42 * 2.309359_r8
                scaleLand  = Lecmwf_ifs_c36r1_T42 * 2.010008_r8
             else
                write(6,'(a,2i7)') f90file//':'//subr// &
                     ': IPARW/IDGRD is unknown: ',iparw,idgrd
                write(6,'(a)')    '  metTYPE: '//trim(metTYPE)
                write(6,'(a,i5)') '  metCYCLE:',metCYCLE
                write(6,'(a,i5)') '  metREVNR:',metREVNR
                stop 'STOP in '//subr
             end if
          else if (IPARW .eq. 320) then
             if (IDGRD .eq. 1) then
                scaleOcean = Oecmwf_ifs_c36r1_T159
                scaleLand  = Lecmwf_ifs_c36r1_T159
             else if (IDGRD .eq. 2) then
                !// Scaling factors based on cy36 2007 meteorology
                scaleOcean = Oecmwf_ifs_c36r1_T159 * 2.210675_r8
                scaleLand  = Lecmwf_ifs_c36r1_T159 * 2.011487_r8
             else if (IDGRD .eq. 4) then
                !// Scaling factors based on cy36 2007 meteorology
                scaleOcean = Oecmwf_ifs_c36r1_T159 * 3.350172_r8
                scaleLand  = Lecmwf_ifs_c36r1_T159 * 2.907207_r8
             else
                write(6,'(a,2i7)') f90file//':'//subr// &
                     ': IPARW/IDGRD is unknown: ',iparw,idgrd
                write(6,'(a)')    '  metTYPE: '//trim(metTYPE)
                write(6,'(a,i5)') '  metCYCLE:',metCYCLE
                write(6,'(a,i5)') '  metREVNR:',metREVNR
                stop 'STOP in '//subr
             end if
          else
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': scaleOcean and scaleLand not '// &
                  'defined for current horizontal resolution'
             stop 'STOP in '//subr
          end if
       else  !// if (metCYCLE .eq. 36 .and. metREVNR .eq. 1) then
          write(6,'(a,2i7)') f90file//':'//subr// &
               ': metCYCLE and metREVNR not not defined'
          write(6,'(a)')    '  metTYPE: '//trim(metTYPE)
          write(6,'(a,i5)') '  metCYCLE:',metCYCLE
          write(6,'(a,i5)') '  metREVNR:',metREVNR
          stop 'STOP in '//subr
       end if

    else if (trim(metTYPE) .eq. 'ECMWF_oIFS' .or. &
             trim(metTYPE) .eq. 'ECMWF_oIFSnc4') then

       !// ECMWF OpenIFS
       if (metCYCLE .eq. 38 .and. metREVNR .eq. 1) then
          !// Cycle 38r1
          if (IPARW .eq. 320) then
             if (IDGRD .eq. 1) then
                scaleOcean = Oecmwf_oifs_c38r1_T159
                scaleLand  = Lecmwf_oifs_c38r1_T159
             else if (IDGRD .eq. 2) then
                !// Scaling factors based on cy38r1 2005 meteorology
                scaleOcean = Oecmwf_oifs_c38r1_T159 * 2.319733_r8
                scaleLand  = Lecmwf_oifs_c38r1_T159 * 2.041230_r8
             else if (IDGRD .eq. 4) then
                !// Scaling factors based on cy38r1 2005 meteorology
                scaleOcean = Oecmwf_oifs_c38r1_T159 * 5.221729_r8
                scaleLand  = Lecmwf_oifs_c38r1_T159 * 4.143434_r8
             else
                write(6,'(a,2i7)') f90file//':'//subr// &
                     ': IPARW/IDGRD is unknown: ',iparw,idgrd
                write(6,'(a)')    '  metTYPE: '//trim(metTYPE)
                write(6,'(a,i5)') '  metCYCLE:',metCYCLE
                write(6,'(a,i5)') '  metREVNR:',metREVNR
                stop 'STOP in '//subr
             end if
          else
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': scaleOcean and scaleLand not '// &
                  'defined for current horizontal resolution'
             stop 'STOP in '//subr
          end if
       else if (metCYCLE .eq. 40 .and. metREVNR .eq. 1) then
          !// Cycle 40r1 - assume same as cy38 for now
          if (IPARW .eq. 320) then
             if (IDGRD .eq. 1) then
                scaleOcean = Oecmwf_oifs_c38r1_T159
                scaleLand  = Lecmwf_oifs_c38r1_T159
             else if (IDGRD .eq. 2) then
                !// Scaling factors based on cy38r1 2005 meteorology
                scaleOcean = Oecmwf_oifs_c38r1_T159 * 2.319733_r8
                scaleLand  = Lecmwf_oifs_c38r1_T159 * 2.041230_r8
             else if (IDGRD .eq. 4) then
                !// Scaling factors based on cy38r1 2005 meteorology
                scaleOcean = Oecmwf_oifs_c38r1_T159 * 5.221729_r8
                scaleLand  = Lecmwf_oifs_c38r1_T159 * 4.143434_r8
             else
                write(6,'(a,2i7)') f90file//':'//subr// &
                     ': IPARW/IDGRD is unknown: ',iparw,idgrd
                write(6,'(a)')    '  metTYPE: '//trim(metTYPE)
                write(6,'(a,i5)') '  metCYCLE:',metCYCLE
                write(6,'(a,i5)') '  metREVNR:',metREVNR
                stop 'STOP in '//subr
             end if
          else
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': scaleOcean and scaleLand not '// &
                  'defined for current horizontal resolution'
             stop 'STOP in '//subr
          end if
       else  !// if (metCYCLE .eq. 38 .and. metREVNR .eq. 1) then
          write(6,'(a,2i7)') f90file//':'//subr// &
               ': metCYCLE and metREVNR not not defined'
          write(6,'(a)')    '  metTYPE: '//trim(metTYPE)
          write(6,'(a,i5)') '  metCYCLE:',metCYCLE
          write(6,'(a,i5)') '  metREVNR:',metREVNR
          stop 'STOP in '//subr
       end if

    else
       write(6,'(a,2i7)') f90file//':'//subr//': metTYPE not defined!'
       write(6,'(a)')    '  metTYPE: '//trim(metTYPE)
       write(6,'(a,i5)') '  metCYCLE:',metCYCLE
       write(6,'(a,i5)') '  metREVNR:',metREVNR
       stop 'STOP in '//subr
    end if !// if (trim(metType)...)

    !// --------------------------------------------------------------------
  end subroutine getScaleFactors
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine LIGHTNING_OAS2015(NDAY,NDAYI,dt_met,LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Lightning NOx
    !//
    !// Algorithm:
    !//   This routine combines parts of GMD2012 method with the
    !//   new UCI2015 treatment. It has better annual cycle in
    !//   flash rates.
    !//
    !// Ole Amund Sovde, March 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: IPARW,IPAR,JPAR,LPAR,LPARW, LWEPAR, IDGRD
    use cmn_ctm, only: ETAA, ETAB, XDGRD, YDGRD, XDEDG, YDEDG, &
         AREAXY, PLAND
    use cmn_chem, only: NEMLIT, LITSRC, LITFAC
    use cmn_met, only: CWETE, CENTU, ZOFLE, T, PRECCNV, P
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input / Output
    integer, intent(in) :: NDAY, NDAYI
    logical, intent(in) :: LNEW_MONTH
    real(r8),  intent(in) :: dt_met

    !// Locals
    integer :: I, J, L
    real(r8) :: CTH      !// Cloud top height (top of convective mass flux)
    real(r8) :: fPR      !// Instant Price&Rind flash rate

    !// Model flash rates over land and ocean. Their units do not
    !// matter, since they are scaled to match observations, but
    !// if you are curious the Price et al. equations are [fl/min].
    real(r8) :: flashL(IPAR,JPAR), flashO(IPAR,JPAR)

    !// We find the flash rates normalised to a climatological year.
    !// The sum of normFrate*dt is 1 for a year when it represents the
    !// climatological year. This is to distribute lightning emissions
    !// (e.g. 5Tg(N)) using normFrate.
    real(r8) :: normFrate(IPAR, JPAR)

    !// normFrate is the sum of flashL and flashO scaled by their
    !// respective factors, which also take care of the normalisation.
    real(r8) :: scaleOcean, scaleLand

    !// Vertical distribution of normFrate(i,j)
    real(r8) :: vertprof(LWEPAR), finc, zfmass(LWEPAR)

    !// Additional factors to calculate flash rates
    real(r8) :: FACT1, FACT2

    !// Meteorological variables
    real(r8) :: LNCWET(LWEPAR), LNCENTU(LWEPAR), TEMP(LWEPAR)
    real(r8) :: ZL(LWEPAR+1)

    real(r8) :: PBOT !// Pressure at bottom of grid box
    real(r8) :: MFLX !// Mass flux
    real(r8) :: MENT !// Mass entrainment [kg/s]
    real(r8) :: minT, maxT !// min/max temperatures in cloud

    !// Diagnostics to get total emissions
    real(r8) :: CNV_SUM, ntgyr
    real(r8) :: NO_SUM(IPAR, JPAR)
    real(r8) :: fLand, fOcean !// Total unscaled flash rates

    !// Counters & indices
    integer :: III, JJJ, IOS, ifnr, NLCL, NLCO !// Counters & indices
    integer :: LFREE, CBL, CTL, LFLX           !// Level indices

    !// Pressures for defining shallow or deep convection.
    real(r8), parameter :: pconvection = 500._r8

    !// Observed lightning flash rated (1/s).
    !// These are taken from  LIS/OTD data LISOTD_HRAC_V2.3.2014.hdf,
    !// as produced by Cecil et al. (2014, J. Atmos. Res.,
    !// vol 135-136, 404-414, doi:10.1016/j.atmosres.2012.06.028).
    !// With an annual global flash rate of 46fl/s, using land
    !// fractions separates this into 9.1fl/s over ocean and
    !// 36.9 over land. If we were to define all grid boxes having
    !// land within 300km reach (roughly similar to Christian et al.
    !// (2003, JGR, 108(D1),doi:10.1029/2002JD002347), these numbers
    !// would be 4.25fl/s and 41.75fl/s, respectively.
    !// Christian et al. (2003) suggested about 5fl/s over ocean.
    !// However, we will use a stricter definition of land, so
    !// we stick to the former.
    real(r8), parameter :: obsFocean=9.1_r8, obsFland=36.9_r8
    real(r8), parameter :: obsFall = obsFocean + obsFland


    !// Scaling parameters depending on meteorological data.
    !//
    !// Price etal equations have units 1/minutes, while observations
    !// are given in 1/s. This conversion is taken care of by the factors
    !// scaleLand and scaleOcean. Thus, after multiplying with the
    !// scaling factors, the units of flash rate is fraction of
    !// all annual lightnings per second. Multiplying by time step
    !// and annual amount of emission (e.g. 5Tg(N)), gives the LNOx
    !// emissions for the time step.
!    !//
!    !// 1997 - 2010 cy36r1 average scaling factor
!    real(r8), parameter :: scaleOceanT42 = 2.304304e-13_r8
!    real(r8), parameter :: scaleLandT42  = 4.329157e-16_r8
!
!    !//
!    !// T159L60 cy36r1 (T159 2007 / T42 2007 * clim.T42)
!    real(r8), parameter :: scaleOceanT159 = 3.594025e-14_r8
!    real(r8), parameter :: scaleLandT159  = 7.061138e-17_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'LIGHTNING_OAS2015'
    !// --------------------------------------------------------------------

    !// If no lightning, return
    if (NEMLIT .eq. 0)  return

    !// Initialize
    CNV_SUM = 0._r8
    NO_SUM(:,:)  = 0._r8

    call getScaleFactors(scaleOcean, scaleLand)

    !// Get correct scaling factors
!    if (IPARW .eq. 128) then
!       if (IDGRD .eq. 1) then
!          scaleOcean = scaleOceanT42
!          scaleLand  = scaleLandT42
!       else if (IDGRD .eq. 2) then
!          !// Scaling factors based on cy36 2007 meteorology
!          scaleOcean = scaleOceanT42 * 1.642020_r8
!          scaleLand  = scaleLandT42  * 1.484240_r8
!       else if (IDGRD .eq. 4) then
!          !// Scaling factors based on cy36 2007 meteorology
!          scaleOcean = scaleOceanT42 * 2.309359_r8
!          scaleLand  = scaleLandT42  * 2.010008_r8
!       else
!          print*,'lightning.f90: IDGRD is unknown:',idgrd
!          stop 'in LIGHTNING_OAS2015'
!       end if
!    else if (IPARW .eq. 320) then
!       if (IDGRD .eq. 1) then
!          scaleOcean = scaleOceanT159
!          scaleLand  = scaleLandT159
!       else if (IDGRD .eq. 2) then
!          !// Scaling factors based on cy36 2007 meteorology
!          scaleOcean = scaleOceanT159 * 2.210675_r8
!          scaleLand  = scaleLandT159  * 2.011487_r8
!       else if (IDGRD .eq. 4) then
!          !// Scaling factors based on cy36 2007 meteorology
!          scaleOcean = scaleOceanT159 * 3.350172_r8
!          scaleLand  = scaleLandT159  * 2.907207_r8
!       else
!          print*,'lightning.f90: IDGRD is unknown:',idgrd
!          stop 'in LIGHTNING_OAS2015'
!       end if
!    else
!       print*,'* lightning.f90: scaleOcean and scaleLand not '// &
!              'defined for current horizontal resolution'
!       stop 'in LIGHTNING_OAS2015'
!    end if

    !// Flash rates
    normFrate(:,:) = 0._r8
    FlashL(:,:) = 0._r8
    FlashO(:,:) = 0._r8

    !// Lightning source
    LITSRC(:,:,:) = 0._r8

    !=================================================================
    ! Read lightning CDF from Ott et al [JGR, 2010]. (ltm, 1/25/11)
    !=================================================================
    if ( LNEW_MONTH .and. NDAY.eq.NDAYI ) then

       !// Initialize lightning profiles
       LITPROFILE(:,:) = 0._r8
       totlitfrq = 0._r8
       totlitn = 0._r8
       accFocean = 0._r8
       accFland  = 0._r8
       accDT     = 0._r8

       !// Echo info
       write(*,*) '   - INIT_LIGHTNING: Reading '//TRIM( INFILE_LIGHTNING )

       !// Get unused file ID
       IFNR = get_free_fileid()

       !// Open file containing lightning PDF data
       open( IFNR, file=trim( INFILE_LIGHTNING ), status='OLD', iostat=IOS)
       if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:1'
         
       !// Read 12 header lines
       do III = 1, 12
          read( IFNR, '(a)', IOSTAT=IOS ) 
          if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:2'
       end do
         
       !// Read NNLIGHT types of lightning profiles
       do III = 1, NNLIGHT
          read( IFNR,*,IOSTAT=IOS) (LITPROFILE(III,JJJ),JJJ=1,NLTYPE)
       end do
         
       !// Close file
       close( IFNR )

       !// Generate flash files for calculating scaling factors
       if (getFlashFiles) then
          open(IFNR,file='flashrate.dta',form='unformatted')
          write(IFNR) IPAR,JPAR
          write(IFNR) XDGRD,YDGRD, AREAXY, PLAND,XDEDG,YDEDG
          close(IFNR)
       end if

       !// Possible extended definition of land grid boxes when
       !// considering lightning.
       !call distland(landlit)

       if (maxval(CENTU) .gt. 0._r8) then
          LentuExist = .true.
       else
          LentuExist = .false.
       end if

    end if !// if ( LNEW_MONTH .and. NDAY.eq.NDAYI ) then

    !// Locate layer of approximately 1500m above surface, which is
    !// about layer 12 in native vertical L60. If collapsing this
    !// will be layer 9 (as in UCI standard code).
    LFREE = 12 - (LPARW-LPAR)


    !=================================================================
    ! Flashrate: Horizontal and temporal distribution
    !=================================================================
!$omp parallel private (I,J,L) &
!$omp          private (CBL,CTL,LFLX) &
!$omp          private (CTH,PBOT) &
!$omp          private (fPR,finc,FACT1,FACT2,minT,maxT,MFLX,MENT) &
!$omp          private (LNCWET,LNCENTU, ZL,TEMP) &
!$omp          private (vertprof) &
!$omp          shared (CWETE,CENTU,T,ZOFLE,PRECCNV) &
!$omp          shared (P,ETAA,ETAB) &
!$omp          shared (LITFAC,LITSRC,normFrate,NO_SUM) &
!$omp          shared (LentuExist, LFREE, JPS, JPN, PLAND) &
!$omp          shared (scaleLand, scaleOcean) &
!$omp          shared (flashO, flashL) &
!$omp          default(NONE)
!$omp do
    do J = 1,JPAR
       do I = 1,IPAR


          !// Store updraft convection.
          !// I do not know the reason behind this parameterisation yet.
          !// Find convection at 850hPa. For all layers above this, we
          !// do not allow convective mass flux to be larger.
          LNCWET(1:(LFREE-1)) = 0._r8
          LNCENTU(1:(LFREE-1)) = 0._r8
          do L = LFREE, LWEPAR
             LNCWET(L) = CWETE(I,J,L)
             LNCENTU(L) = CENTU(I,J,L)
          end do

          !// Use rain as screening. There may be some convection
          !// within large scale rain, without producing convective
          !// rain falling to the ground. It is not clear to what extent
          !// that can happen for vigorous convection. Therefore I
          !// check only convective rain at the surface; it seems to
          !// be a filter that works well to get temporal distribution
          !// fairly correct.
          if (PRECCNV(I,J,1) .eq. 0._r8) cycle

          !// There may be some instances where convective mass flux
          !// occur in single layers or separated into several intervals.
          !// This is due to the nature of the IFS model physics.
          !// If e.g. level 10 and 20 are the only layers with
          !// convective mass flux, there should be no lightning.
          call filterLNC(I, J, LFREE, LNCWET, LNCENTU, CBL, CTL)


          !// Check if there are any entraining updrafts;
          !// if none, then cycle to next box.
          if (CTL .eq. 0) cycle


          !// There are several scaling factors or additional filters
          !// that should be included, to sort out the flash rates.
          !//
          !// 1. Temperature in cloud: Cloud should consist of both
          !//    liquid and ice. This can be checked with temperature.
          !// 2. Cloud fraction scaling: The convective cell does not
          !//    cover the whole grid box.
          !// 3. Updraft velocity: Lightning needs large vertical
          !//    velocities in the convective cell.
          !// 4. Area scaling: A convective cloud with a certain cloud
          !//    top height may cover a larger area at low latitudes,
          !//    giving higher flash rate there.
          !// 5. Flux scaling: This may take horizontal extent of the
          !//    convective cloud into account.
          !// Final: scaleLand and scaleOcean: to scale to observed flash
          !//        rates (i.e. normalised to climatological year).

          !// Retrieve temperature as 1D array
          do L = 1, LWEPAR
             TEMP(L) = T(I,J,L)
          end do



          !// 1. Temperature
          !// ------------------------------------------------------------
          !// Make sure the cloud is mixed phase, containing both
          !// liquid and ice. We assume this occurs when minimum
          !// temperature in the column is colder than -40C and
          !// the maximum temperature is warmer than 0C.
          !// To allow for some uncertainty, UCI scales the first
          !// factor linearly from 0 at -30C to 1 at -40C, and the
          !// latter from 1 at 0C to 0 at -20C.
          !//
          !// I think the first scaling is ok, but the latter scaling
          !// in the UCI method, uses surface temperature. In any case
          !// I think allowing lightning down to -20C sounds very weird.
          !//
          !// A better choice is to require the cloud to have max
          !// temperatures above -2.5C before producing lightning;
          !// taking into account some uncertainty in temperature: The
          !// metdata temperature is gridbox average temperatures, not
          !// necessarily the correct updraft temperatures.
          !// Strict temperature control may not very physical, but could
          !// be an alternative.
          !//
          !// Surface temperature scalings?
          !// Using the surface temperature to scale lightning, will
          !// probably give wrong diurnal cycle; there should be more
          !// flashes at local night time. Instead of such a scaling,
          !// I use a linear scale from 0 at -2.5C to 1 at 7.5C.
          !// In this way we should increase the more vigorous lightning
          !// and reduce the lightning from relatively cold clouds.

          minT = minval(TEMP(CBL:CTL))
          maxT = maxval(TEMP(CBL:CTL))

          !// Linear scaling: T>=7.5: FACT1=1, T<-2.5: FACT1=0
          FACT1 = max(0._r8, min(0.1_r8*(maxT - 270.5_r8), 1._r8) )
          !// Check if criteria is met
          if ( FACT1 .eq. 0._r8 ) cycle

          !// Linear scaling for cloud minimum temperature:
          !// T=<-40C: FACT2=1, T>-30C: FACT2=0
          FACT2 = max(0._r8, min(-0.1_r8*(minT - 243._r8), 1._r8) )
          !// Check if criteria is met
          if ( FACT2 .eq. 0._r8 ) cycle


          !// 2. Cloud fraction scaling?
          !// ------------------------------------------------------------
          !// The flux itself, in values of kg/s may be less meaningful;
          !// when the same flux covers a smaller area, convection should
          !// be more rigorous.
          !// But convection covering a larger area is likely producing
          !// more lightning. One way to reduce lightning frequency over
          !// small areas is to multiply with cloud fraction.
          !//
          !// But, unfortunately, cloud fraction is an instant field,
          !// whereas convective mass flux is based on accumulated values,
          !// so there may be times where we have convective mass flux
          !// but no clouds. Testing so far indicates that including
          !// this factor makes the horizontal distribution worse.



          !// 3. Updraft velocity
          !// ------------------------------------------------------------
          !// For a thunderstorm convective cell, updraft velocities
          !// typically reach > 5m/s. However, since we do not know the
          !// size of the convective cell (cloud fraction cannot be
          !// used since it is instant and flux is accumulated), we
          !// could in principle use the grid box mass. That would
          !// not give a  "real" updraft velocity, but how fast the
          !// the mass in the grid box would move to match the flux.
          !// -> This would be resolution dependent, so the velocity
          !//    would have different meaning in different resolutions.
          !//
          !// A possible work-around is to calculate the fraction of
          !// moved mass to some fixed mass value that does not change
          !// across resolutions. As a good proxy, we can use the total
          !// mass per level. The total mass per level should
          !// change minimally between resolutions.
          !// But unfortunately, finding a certain limit that work well
          !// in different resolutions is tricky. I tried it, but could
          !// not get it working.


          !// Restrict lightning to be produced only when convection
          !// reaches a height where pressure is less than some
          !// limit (pconvection).
          !// Also require that entrainment occurs above this level (it
          !// usually does).
          LFLX = 0
          do L = CBL, CTL
             PBOT = ETAA(L) + P(I,J)*ETAB(L)
             if (LentuExist) then
                MENT = sum(LNCENTU(L:CTL))
                if (PBOT .le. pconvection .and. MENT .gt. 0._r8) then
                   LFLX = L
                   exit           !// Exit L-loop
                end if
             else
                if (PBOT .le. pconvection) then
                   LFLX = L
                   exit           !// Exit L-loop
                end if
             end if
          end do

          !// Skip lightning if convection does not reach the limit.
          if (LFLX .eq. 0) cycle


          !// 4. Area scaling?
          !// ------------------------------------------------------------
          !// We find frequency based on cloud top height, not considering
          !// the horizontal extent of the convection. A cloud with a
          !// certain size should have the same frequency at Equator as at
          !// high latitudes,
          !// However, at Equator the convective cloud may generally be
          !// larger, or there may even be several convective clouds in
          !// the same grid box. So it can be argued that the modelled
          !// frequency should be weighted by the fractional area relative
          !// to Equator.
          !// Again, cloud fraction would solve this issue.
          !// To take into account the horizontal extent of the convective
          !// plume, we might scale with the convective mass flux at some
          !// level. They should not be used together, because the flux
          !// do carry some information about size.


          !// 5. Flux scaling
          !// ------------------------------------------------------------
          !// If the flux is large, the convective cloud may be larger and
          !// should produce more lightning. We can scale the flash rates
          !// with the flux at a certain level, e.g. at 500hPa.
          !// This produces more lightning at high latitudes than we want.
          !//
          !// I have also found that it could be smart to screen away
          !// instances where the entrainment below this level
          !// (i.e. at LFLX-1) is small, e.g. 10% of the flux at LFLX.
          !// However, it seems that in cy36r1 had some changes to the
          !// entrainment from 2009-2012. I remember the ECMWF corrected
          !// a few things, but do not know if that should affect us.
          !// So I will skip this screening for now, only require
          !// entrainment > 0.
          !// Flux into LFLX level
          !MFLX = LNCWET(LFLX)
          !// Entrainment into LFLX-1 level
          !MENT = LNCENTU(LFLX-1)


          !// Final: scaleLand and scaleOcean
          !// ------------------------------------------------------------
          !// Price etal equations have units 1/minutes. The final scaling
          !// convert to 1/s, and also convert to fraction of climatological
          !// mean, so that the sum for a year on average will be 1.
          !// This is taken care of by the factors scaleLand and
          !// scaleOcean.


          !// Find layer bottom height above ground [km]
          !// ------------------------------------------------------------
          do L = 1, LWEPAR+1
             ZL(L) = 1.e-3_r8 * (ZOFLE(L,I,J) - ZOFLE(1,I,J))
          end do

          !// Now use Price et al. equations
          !// ------------------------------------------------------------
          !// Cloud top height above ground is top of level CTL [km]
          CTH = ZL(CTL+1)

          if (PLAND(I,J) .gt. 0.25_r8) then
             fPR = CTH**4.92_r8 !// 3.44d-5 included in scaleLand
             finc = FACT1 * FACT2
             normFrate(I,J) = finc * fPR * scaleLand
             flashL(I,J) = finc * fPR !// Save unscaled land
          else
             fPR = CTH**1.73_r8 !// 6.4d-4 included in scaleOcean
             finc = FACT1 * FACT2
             normFrate(I,J) = finc * fPR * scaleOcean
             flashO(I,J) = finc * fPR !// Save unscaled oceanic
          end if



          !---------------------------------------------------------
          ! NOx emission, vertical distribution
          !---------------------------------------------------------

          ! If there's lightning in the column
          if ( normFrate(I,J) .gt. 0._r8 ) then

             !// Partition NOx emissions vertically using the PDFs from
             !// Ott et al. (2010) and regional definitions from
             !// Allen et al. (2010)
             call lightdist(I, J, CTL, ZL, normFrate(I,J), vertprof)

             !// Assign emissions into array for later use. 
             !// (units: fraction of annual emission per second)
             do L=1,CTL
                LITSRC(L,I,J) = vertprof(L)
                !// Sum up NO to test 5 Tg(N)/year
                NO_SUM(I,J) = NO_SUM(I,J) + (vertprof(L)*LITFAC(1))
             end do

          end if

       end do
    end do !// do J = 1,JPAR
!$omp end do
!$omp end parallel

    !=================================================================
    ! Diagnostics for lightning production
    !=================================================================

    !// Total unscaled land and ocean flashes
    fLand  = sum(flashL)
    fOcean = sum(flashO)

    !// Flash files
    if (getFlashFiles) then
       IFNR = get_free_fileid()
       open(IFNR,file='flashrate.dta',form='unformatted', position='append')
       write(IFNR) real(flashL, r4), real(flashO, r4)
       close(IFNR)
    end if

    !// How many grid boxes experienced lightning?
    NLCL = 0
    NLCO = 0
    do J = 1, JPAR
       do I = 1, IPAR
          if (flashL(I,J) .gt. 0._r8) NLCL = NLCL + 1
          if (flashO(I,J) .gt. 0._r8) NLCO = NLCO + 1
       end do
    end do
    write(*,'(A,3i7)') 'Lightning: # grid boxes w/flash (L/O/L+O): ', &
         NLCL, NLCO, NLCL+NLCO

    !// Total global fraction of annual total flashrate. Must
    !// multiply by sum up to make annual total 1
    cnv_sum = sum(normFrate) * dt_met

    !// Accumulate this fraction
    totlitfrq = totlitfrq + cnv_sum

    !// Sum up L-NOx until next update of LITSRC (duration dt_met)
    totlitn = totlitn + sum(no_sum) * dt_met * 1.e-9_r8 * 14._r8/30._r8
    !// Emissions this step, scaled to Tg(N)/yr
    ntgyr = sum(NO_SUM) * 1.e-9_r8 * 14._r8/30._r8 * dtyear

    !// Accumulate flashes during time step and total time.
    accFland = accFland + fLand * dt_met
    accFocean = accFocean + fOcean * dt_met
    accDT = accDT + dt_met

    !// Write total flash rate for ocean and land
    !  write(*,'(A,f12.3,1x,f12.3)')
    ! &     'Lightning: Flashrates (1/s) ocean/land:',
    ! &     fOcean*scaleOcean, fLand*scaleLand

    !// Fraction of annual total flashrates
    write(*,'(A,es23.17)') 'Lightning: Fraction of annual flashes: ', cnv_sum

    !// Emissions during this time step (as Tg(N)/yr) and
    !// also the accumulated amount
    write(*,'(A,2f12.5)') 'Lightning: Tg(N)/yr & accumulated Tg(N):', &
         ntgyr, totlitn

    !// Accumualted flash frequency fraction of annual total, both
    !// new and old values. Their difference is cnv_sum.
    write(*,'(A,es23.17,1x,es23.17)') 'Lightning: Accu.flashfreq n/o: ', &
         totlitfrq,totlitfrq-cnv_sum

    !// Scaling factors. These are the values you need to modify
    !// scaleOean and scaleLand if you use other meteorological data.
    !// IMPORTANT:
    !// Remember that you need to use at least an annual mean, but
    !// preferably a mean over several years of data.
    !//
    !// Note also that for a single year of meteorological data,
    !// this printout may differ slightly from the applied scaling
    !// factors. This does not mean that the scalings are wrong,
    !// it may be that the applied scalings were calculated using
    !// a different year or several years.
    !//
    !// Printout is done after calculating lightning for the 
    !// last meteorological time step of each day, giving the
    !// average scalings for days from NDAYI through NDAY.
    !//
    !// Will find scaling so that
    !//   annual_modelflashes * scale = 1.
    !// separated into land and ocean:
    !//   scaleO = obsFocean/obsFall / annual_model_Oflashes
    !//   scaleL = obsFland/obsFall  / annual_model_Lflashes
    if (mod(accDT,86400._r8) .eq. 0._r8) &
         write(*,'(A,i4,es16.7,1x,es16.7)') &
         'Lightning: Accu.scal. NDAY/ocean/land: ',NDAY, &
         obsFocean/obsFall / (accFocean/accDT * dtyear), &
         obsFland/obsFall / (accFland/accDT * dtyear)

    !// Sensiblitity check
    if (ntgyr .lt. 1._r8 .or. ntgyr .gt. 12._r8) then
       print*,'*** WARNING: lightning.f90: Low/High LNOx: '// &
              'Check scaleLand and scaleOcean!'
    end if

    !// --------------------------------------------------------------------
  end subroutine LIGHTNING_OAS2015
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine LIGHTNING_UCI2015(NDAY,NDAYI,dt_met,LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Comment by Ole Amund Sovde:
    !// This method is NOT to be used in CTM3. While the method
    !// produces fairly ok annual horizontal distribution of flashes,
    !// the temporal distribution is almost gone. In addition, the
    !// number of lightning days are very high in many places, giving
    !// small amounts of lightning emissions very often.
    !// Better when filtering convective events using the filterLNC.
    !//
    !// About this routine and why to dismiss it:
    !// This is from the qcode version 7.1. It uses the Price and Rind
    !// (1992) equations slightly differently than before. Their method
    !// was updated so that that the flash frequency in a 2x2 window
    !// should be equal to the flash frequency summed up over its four
    !// 1x1 boxes.
    !// Clearly the max level of convection will be the same.
    !// Say we have fluxes F1, F2, F3, F4 and cloud tops
    !// H1, H2, H3, H4. Ideally, we could have:
    !//    FL = (F1*H1^m + F2*H2^m + F3*H3^m + F4*H4^m)/(F1+F2+F3+F4)
    !// where m is the Price etal exponents.
    !//
    !// But how should we determine the individual F-values? When
    !// running 2x2 we only have F=F1+F2+F3+F4 available.
    !// The UCI method tries to look sum up several "clouds", where
    !// it is assumed that each time the flux is reduced, a new
    !// cloud top is reached. Thus, if the flux decreases in 10
    !// consecutive layers, all 10 will count. Probably the largest
    !// will give the most contribution. UCI also considers the
    !// temperature of the layer above each cloud top, which needs
    !// to be cold.
    !// To me this seems to cover only a part of the picture. What they
    !// do is a quasi-calculation of "detrainment", and seems not
    !// completely correct.
    !//
    !// I have tried to calculate the true detrainment from LFLX
    !// upwards, and then scale the normalized flashrate with the
    !// sum of detrainment, but to no help. Dropping it seems better.
    !//
    !// For UCI method, using daily flashes f(i,j,d), scaled to match
    !// 1.45d9 fl/yr, i.e. sum(f) = 1. We multiply f by 5d9*30/14 to
    !// get total NO, and then plot e.g. Oslo (i=5,j=54 in T42).
    !// This gives about 160 lightning days, where ~60 have 1d3kg(NO)
    !// and ~5 have 1d4kg(NO).
    !// 100kg(NO) per day occurs very frequently. Is that a problem?
    !// Probably not very, but if the background is clean, it may very
    !// well have an impact.
    !//
    !// 
    !// The global lightning flash rate is observed to be 46 flashes/s,
    !// of which 76% (35 flashes/s) occurs over land. This is derived
    !// from the space-based OTD and LIS sensors, which operated over
    !// 1998-2005. We scale the flash rate derived from cloud top height
    !// and other criteria to match this mean rate over land and ocean.
    !// The scale factors need to be recalculated for each change of
    !// meteorology version and resolution, or whenever other parts of
    !// the algorithm change. CDH calculated factors for EC cycle 36r1,
    !// T42L60, over 1999-2009.
    !// --------------------------------------------------------------------
    !// Horizontal distribution is based on:
    !//   Price & Rind (1992), doi:10.1029/92JD00719
    !// Vertical distribution is based on:
    !//   Ott et al (2010), doi:10.1029/2009JD011880
    !//
    !// Ole Amund Sovde, March 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: IPARW,IPAR,JPAR,LPAR,LPARW, LWEPAR, IDGRD
    use cmn_ctm, only: ETAA, ETAB, XDGRD, YDGRD, XDEDG, YDEDG, &
         AREAXY, PLAND
    use cmn_chem, only: NEMLIT, LITSRC, LITFAC
    use cmn_met, only: CWETE, CENTU, ZOFLE, T, PRECCNV, P
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input / Output
    integer, intent(in) :: NDAY, NDAYI
    logical, intent(in) :: LNEW_MONTH
    real(r8),  intent(in) :: dt_met

    !// Locals
    integer :: I, J, L

    !// Model flash rates over land and ocean. Their units do not
    !// matter, since they are scaled to match observations, but
    !// if you are curious the Price et al. equations are [fl/min].
    real(r8) :: flashL(IPAR,JPAR), flashO(IPAR,JPAR)

    !// We find the flash rates normalised to a climatological year.
    !// The sum of normFrate*dt is 1 for a year when it represents the
    !// climatological year. This is to distribute lightning emissions
    !// (e.g. 5Tg(N)) using normFrate.
    real(r8) :: normFrate(IPAR, JPAR)

    !// normFrate is the sum of flashL and flashO scaled by their
    !// respective factors, which also take care of the normalisation.
    real(r8) :: scaleOcean, scaleLand

    !// Vertical distribution of normFrate(i,j)
    real(r8) :: vertprof(LWEPAR)

    !// Additional factors to calculate flash rates
    real(r8) :: FACT1, FACT2, finc, fPR, MENT, PBOT

    !// Meteorological variables
    real(r8) :: LNCWET(LWEPAR), TEMP(LWEPAR)
    real(r8) :: ZL(LWEPAR+1)
    real(r8) :: LNCENTU(LWEPAR) !// Used in filterLNC

    !// Diagnostics to get total emissions
    real(r8) :: CNV_SUM, ntgyr
    real(r8) :: NO_SUM(IPAR, JPAR)
    real(r8) :: fLand, fOcean !// Total unscaled flash rates

    !// Counters & indices
    integer :: III, JJJ, IOS, ifnr, NLCL, NLCO !// Counters & indices
    integer :: LFREE, CBL, CTL, LFLX      !// Level indices

    !// Logical switch to define polar grid boxes
    logical :: LPOLES

    !// Pressures for defining shallow or deep convection.
    real(r8), parameter :: pconvection = 500._r8

    !// Observed lightning flash rate (1/s), explained in subroutine
    !// LIGHTNING_GMD2012.
    real(r8), parameter :: obsFocean=11._r8, obsFland=35._r8
    real(r8), parameter :: obsFall = obsFocean + obsFland

    !// Scaling parameters depending on meteorological data.
    !// Price and Rind equations have units 1/minutes, while observations
    !// are given in 1/s. This conversion is taken care of by the factors
    !// scaleLand and scaleOcean. Thus, after multiplying with the
    !// scaling factors, the units of flash rate is 1/s.
    !//
    !// 1997 - 2010 cy36r1 average scaling factor
    !// These differ slightly from UCI numbers
    real(r8), parameter :: scaleOceanT42 = 2.406103e-22_r8
    real(r8), parameter :: scaleLandT42  = 5.919869e-25_r8

    !//
    !// T159L60 cycle 36r1 (not computed)
    real(r8), parameter :: scaleOceanT159 = 2.e-22_r8
    real(r8), parameter :: scaleLandT159  = 6.e-25_r8

    !// Use filterLNC?
    logical, parameter :: LuseFilterLNC=.false.
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'LIGHTNING_UCI2015'
    !// --------------------------------------------------------------------

    !// If no lightning, return
    if (NEMLIT .eq. 0)  return

    !// Initialize
    CNV_SUM = 0._r8
    NO_SUM(:,:)  = 0._r8

    !// Get correct scaling factors
    if (IPARW .eq. 128) then
       if (IDGRD .eq. 1) then
          scaleOcean = scaleOceanT42
          scaleLand  = scaleLandT42
       else
          write(6,'(a,i5)') f90file//':'//subr// &
               ': IDGRD is unknown:',idgrd
          stop 'STOP in '//subr
       end if
    else if (IPARW .eq. 320) then
       write(6,'(a,i5)') f90file//':'//subr// &
            ': IDGRD is unknown:',idgrd
       stop 'STOP in '//subr
    else
       write(6,'(a,i5)') f90file//':'//subr// &
            ': scaleOcean and scaleLand not '// &
            'defined for current horizontal resolution'
       stop 'STOP in '//subr
    end if

    !// Flash rates
    normFrate(:,:) = 0._r8
    FlashL(:,:) = 0._r8
    FlashO(:,:) = 0._r8

    !// Lightning source
    LITSRC(:,:,:) = 0._r8

    !=================================================================
    ! Read lightning CDF from Ott et al [JGR, 2010]. (ltm, 1/25/11)
    !=================================================================
    if ( LNEW_MONTH .and. NDAY.eq.NDAYI ) then

       !// Initialize lightning profiles
       LITPROFILE(:,:) = 0._r8
       totlitfrq = 0._r8
       totlitn = 0._r8
       accFocean = 0._r8
       accFland  = 0._r8
       accDT     = 0._r8
       !// Poles: south of 60S and north of 60N are counted as ocean
       J  = 1
       do while (YDEDG(J) .lt. -60._r8)
          J = J + 1
       end do
       JPS = J
       JPN = JPAR+1 - JPS

       !// Echo info
       write(*,*) '   - INIT_LIGHTNING: Reading '//TRIM( INFILE_LIGHTNING )

       !// Get unused file ID
       IFNR = get_free_fileid()

       !// Open file containing lightning PDF data
       open( IFNR, FILE=trim( INFILE_LIGHTNING ), STATUS='OLD', IOSTAT=IOS)
       if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:1'
         
       !// Read 12 header lines
       do III = 1, 12
          read( IFNR, '(a)', IOSTAT=IOS ) 
          if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:2'
       end do
         
       !// Read NNLIGHT types of lightning profiles
       do III = 1, NNLIGHT
          read( IFNR,*,IOSTAT=IOS) (LITPROFILE(III,JJJ),JJJ=1,NLTYPE)
       end do
         
       !// Close file
       close( IFNR )

       !// Flash files
       if (getFlashFiles) then
          open(IFNR,file='flashrate.dta',form='unformatted')
          write(IFNR) IPAR,JPAR
          write(IFNR) XDGRD,YDGRD, AREAXY, PLAND, XDEDG,YDEDG
          close(IFNR)
       end if

       if (maxval(CENTU) .gt. 0._r8) then
          LentuExist = .true.
       else
          LentuExist = .false.
       end if


    end if !// if ( LNEW_MONTH .and. NDAY.eq.NDAYI ) then


    !// Locate layer of approximately 1500m or 850hPa, which is
    !// about layer 12 in native vertical L60. If collapsing this
    !// will be layer 9 (as in UCI standard code).
    LFREE = 12 - (LPARW-LPAR)

!$omp parallel private (I,J,L) &
!$omp         private (LPOLES,CTL,LNCWET,ZL,TEMP) &
!$omp         private (fPR,FACT1,FACT2) &
!$omp         private (vertprof) &
!$omp         private (LFLX,CBL,LNCENTU,MENT,PBOT) &
!$omp         shared (CENTU,PRECCNV) &
!$omp         shared (CWETE,T,ZOFLE) &
!$omp         shared (P,ETAA,ETAB) &
!$omp         shared (LITFAC,LITSRC,normFrate,NO_SUM) &
!$omp         shared (LentuExist, LFREE, JPS, JPN, PLAND) &
!$omp         shared (scaleLand, scaleOcean) &
!$omp         shared (flashO, flashL) &
!$omp         default(NONE)
!$omp do
    do J = 1,JPAR
       LPOLES = J.lt.JPS .or. J.gt.JPN
       do I = 1,IPAR
          !---------------------------------------------------------
          ! Flashrate: Horizontal and temporal distribution
          !---------------------------------------------------------


          !// Store updraft convection.
          LNCWET(1:(LFREE-1)) = 0._r8
          do L = LFREE, LWEPAR
             LNCWET(L) = CWETE(I,J,L)
          end do

          !// Extra stuff if filterLNC is to be used
          if (LuseFilterLNC) then
             LNCENTU(1:(LFREE-1)) = 0._r8
             do L = LFREE, LWEPAR
                LNCENTU(L) = CENTU(I,J,L)
             end do
             if (PRECCNV(I,J,1) .eq. 0._r8) cycle
             call filterLNC(I, J, LFREE, LNCWET, LNCENTU, CBL, CTL)
             if (CTL .eq. 0) cycle
          end if



          !// Find convection at ~1500m. For all layers above this,
          !// do not allow convective mass flux to be larger.
          !// OAS: I think this is some obscure way to locate levels
          !// where there is detrainment, which can be treated as
          !// a potential cloud top, and that all those levels make up
          !// several cloud tops which may produce lightning.
          if (LuseFilterLNC) then
             LFLX = 0
             do L = CBL, CTL
                PBOT = ETAA(L) + P(I,J)*ETAB(L)
                if (LentuExist) then
                   MENT = sum(LNCENTU(L:CTL))
                   if (PBOT .le. pconvection .and. MENT .gt. 0._r8) then
                      LFLX = L
                      exit           !// Exit L-loop
                   end if
                else
                   if (PBOT .le. pconvection) then
                      LFLX = L
                      exit           !// Exit L-loop
                   end if
                end if
             end do
             !// Skip lightning if convection does not reach the limit.
             if (LFLX .eq. 0) cycle

             !// If filterLNC is to be used, adjust LNCWET ala UCI.
             do L = CBL+1, CTL
                LNCWET(L) = min(LNCWET(L), LNCWET(L-1))
             end do
          else
             !// Standard UCI2015 treatment.
             !// OAS: Actually, this may be slightly flawed, since if
             !// convection starts above LFREE, LNCWET will become zero.
             do L = LFREE+1, LWEPAR
                LNCWET(L) = min(LNCWET(L), LNCWET(L-1))
             end do
             !// locate the top layer with convection
             CTL = 0
             do L=1,LWEPAR-1
                if ( LNCWET(L) .gt. 0._r8) then
                   CTL = L+1 !// this is actually the layer above top conv.
                end if
             end do
             !// Override
             if (CTL .gt. LWEPAR) CTL = LWEPAR
             !// Check if there are any entraining updrafts
             !// If none, then cycle to next box
             if (CTL .eq. 0) cycle
          end if



          !// Retrieve temperature as 1D array
          do L = 1,LWEPAR
             TEMP(L) = T(I,J,L)
          end do
            
          !// Make sure the cloud top is colder than -40C
          !// and the surface is warmer than 0C. This guarantees
          !// that the cloud contains mixed liquid and ice phases.
          FACT1 = max(0._r8, min(0.05_r8*(TEMP(1) - 253._r8), 1._r8) )
          if ( FACT1 .eq. 0._r8) cycle
          if ( TEMP(CTL) .ge. 243._r8 ) cycle

          !// Find layer bottom height above ground [km]
          do L = 1,LWEPAR+1
             ZL(L) = 1.e-3_r8 * (ZOFLE(L,I,J) - ZOFLE(1,I,J))
          end do


          if (.not. LuseFilterLNC) then
             CBL = LFREE
          end if

          !// Use Price and Rind formulas for land and ocean
          !// flash frequency.
          !// Use additional scale factors for land and ocean to match
          !// the OTD-LIS climatology for flashes over land (35/s) 
          !// and ocean (11/s)
          !//
          !// Price and Rind equations have units 1/minutes, which we
          !// convert to 1/s. This is taken care of by the factors
          !// scaleLand and scaleOcean.
          if (PLAND(I,J) .gt. 0.25_r8 .and. .not.LPOLES) then
             !// Flashrate is initialised before IJ-loop
             fPR = 0._r8
             do L = CBL, CTL-1
                if (LNCWET(L+1) .lt. LNCWET(L)) then
                   FACT2 = max(0._r8, min(-0.1_r8*(TEMP(L+1)-243._r8), 1._r8) )
                   fPR = fPR &
                       + FACT2 * (LNCWET(L) - LNCWET(L+1)) * ZL(L)**4.9_r8
                else
                   !// Actually, this is a bit silly, because if LNCWET
                   !/// at L+1 /= L, they are the same, as set at the
                   !// beginning of the routine.
                   fPR = fPR &
                       + (LNCWET(L) - LNCWET(L+1)) * ZL(L)**4.9_r8
                end if
             end do
             normFrate(I,J) = fPR * FACT1 * scaleLand
             flashL(I,J) = fPR * FACT1 !// Save unscaled oceanic
          else
             fPR = 0._r8
             do L = CBL, CTL-1
                if (LNCWET(L+1) .lt. LNCWET(L)) then
                   FACT2 = max(0._r8, min(-0.1_r8*(TEMP(L+1)-243._r8), 1._r8) )
                   fPR = fPR &
                       + FACT2 * (LNCWET(L) - LNCWET(L+1)) * ZL(L)**1.73_r8
                else
                   fPR = fPR &
                       + (LNCWET(L) - LNCWET(L+1)) * ZL(L)**1.73_r8
                end if
             end do
             normFrate(I,J) = fPR * FACT1 * scaleOcean
             flashO(I,J) = fPR * FACT1 !// Save unscaled oceanic
          end if


          !---------------------------------------------------------
          ! NOx emission, vertical distribution
          !---------------------------------------------------------

          ! If there's lightning in the column
          if ( normFrate(I,J) .gt. 0._r8 ) then

             !// Partition NOx emissions vertically using the PDFs from
             !// Ott et al. (2010) and regional definitions from
             !// Allen et al. (2010)
             call lightdist(I, J, CTL, ZL, normFrate(I,J), vertprof)

             !// Assign emissions into array for later use. 
             !// (units: fraction of annual emission per second)
             do L=1,CTL
                LITSRC(L,I,J) = vertprof(L)
                !// Sum up NO to test 5 Tg(N)/year
                NO_SUM(I,J) = NO_SUM(I,J) + (vertprof(L)*LITFAC(1))
             end do

          end if !// if ( normFrate(I,J) .gt. 0._r8 ) then

       end do !// do I = 1,IPAR
    end do !// do J = 1,JPAR
!$omp end do
!$omp end parallel

    !// Total unscaled land and ocean flashes
    fLand  = sum(flashL)
    fOcean = sum(flashO)

    !// Flash files
    if (getFlashFiles) then
       IFNR = get_free_fileid()
       open(IFNR,file='flashrate.dta',form='unformatted', position='append')
       write(IFNR) real(flashL, r4), real(flashO, r4)
       close(IFNR)
    end if

    !// How many grid boxes experienced lightning?
    NLCL = 0
    NLCO = 0
    do J = 1, JPAR
       do I = 1, IPAR
          if (flashL(I,J) .gt. 0._r8) NLCL = NLCL + 1
          if (flashO(I,J) .gt. 0._r8) NLCO = NLCO + 1
       end do
    end do
    write(*,'(A,3i7)') 'Lightning: # grid boxes w/flash (L/O/L+O): ', &
         NLCL, NLCO, NLCL+NLCO

    !// Total global fraction of annual total flashrate. Must
    !// multiply by sum up to make annual total 1
    cnv_sum = sum(normFrate) * dt_met
    !// Accumulate this fraction
    totlitfrq = totlitfrq + cnv_sum

    !// Sum up L-NOx until next update of LITSRC (duration dt_met)
    totlitn = totlitn + sum(no_sum) * dt_met * 1.e-9_r8 * 14._r8/30._r8
    !// Emissions this step, scaled to Tg(N)/yr
    ntgyr = sum(NO_SUM) * 1.e-9_r8 * 14._r8/30._r8 * dtyear

    !// Accumulate flashes during time step and total time.
    accFland = accFland + fLand * dt_met
    accFocean = accFocean + fOcean * dt_met
    accDT = accDT + dt_met

    !// Write total flash rate for ocean and land
    !  write(*,'(A,f12.3,1x,f12.3)')
    ! &     'Lightning: Flashrates (1/s) ocean/land:',
    ! &     fOcean*scaleOcean, fLand*scaleLand

    !// Fraction of annual total flashrates
    write(*,'(A,es23.17)') 'Lightning: Fraction of annual flashes: ', cnv_sum

    !// Emissions during this time step (as Tg(N)/yr) and
    !// also the accumulated amount
    write(*,'(A,2f12.5)') 'Lightning: Tg(N)/yr & accumulated Tg(N):', &
          ntgyr, totlitn

    !// Accumualted flash frequency fraction of annual total, both
    !// new and old values. Their difference is cnv_sum.
    write(*,'(A,es23.17,1x,es23.17)') 'Lightning: Accu.flashfreq n/o: ', &
         totlitfrq,totlitfrq-cnv_sum

    !// Scaling factors. These are the values you need to modify
    !// scaleOean and scaleLand if you use other meteorological data.
    !// IMPORTANT:
    !// Remember that you need to use at least an annual mean, but
    !// preferably a mean over several years of data.
    !//
    !// Note also that for a single year of meteorological data,
    !// this printout may differ slightly from the applied scaling
    !// factors. This does not mean that the scalings are wrong,
    !// it may be that the applied scalings were calculated using
    !// a different year or several years.
    !//
    !// Printout is done after calculating lightning for the 
    !// last meteorological time step of each day, giving the
    !// average scalings for days from NDAYI through NDAY.
    !//
    !// Will find scaling so that
    !//   annual_modelflashes * scale = 1.
    !// separated into land and ocean:
    !//   scaleO = obsFocean/obsFall / annual_model_Oflashes
    !//   scaleL = obsFland/obsFall  / annual_model_Lflashes
    if (mod(accDT,86400._r8) .eq. 0._r8) &
         write(*,'(A,i4,es16.7,1x,es16.7)') &
         'Lightning: Accu.scal. NDAY/ocean/land: ',NDAY, &
         obsFocean/obsFall / (accFocean/accDT * dtyear), &
         obsFland/obsFall / (accFland/accDT * dtyear)

    !// Sensiblitity check
    if (ntgyr .lt. 1._r8 .or. ntgyr .gt. 12._r8) then
       print*,'*** WARNING: lightning.f90: Low/High LNOx: '// &
              'Check scaleLand and scaleOcean!'
    end if

    !// --------------------------------------------------------------------
  end subroutine LIGHTNING_UCI2015
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine LIGHTNING_GMD2012(NDAY,NDAYI,dt_met,LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// This is the GMD2012 routine.
    !//
    !// The horizontal and temporal lightning distribution is a polynomial
    !// function of cloud-top height (Price and Rind, 1997) with additional
    !// criteria for temperature and updraft velocities. 
    !// Lightning flash rates are scaled to match the observed climatology,
    !// as described below, and a constant NOx yield per flash is used.
    !// The vertical distribution of NOx emissions are from Ott et al. (2010).
    !//
    !// Scaling factors:
    !// The global lightning flash rate is observed to be 46 flashes/s,
    !// of which 76% (35 flashes/s) occurs over land. This is derived from
    !// the space-based OTD and LIS sensors, which operated over 1998-2005. 
    !// We scale the flash rate derived from cloud top height and other
    !// criteria to match this mean rate over land and ocean. The scale
    !// factors need to be recalculated for each change of meteorology version
    !// and resolution, or whenever other parts of the algorithm change.
    !// CDH calculated factors for EC cycle 36, T42L60, over 1999-2009.
    !// --------------------------------------------------------------------
    !// Horizontal distribution is based on:
    !//   Price & Rind (1992), doi:10.1029/92JD00719
    !// Vertical distribution is based on:
    !//   Ott et al (2010), doi:10.1029/2009JD011880
    !//
    !// Modified from file I got from Chris D. Holmes at UCI.
    !//
    !// Ole Amund Sovde, August 2011
    !// --------------------------------------------------------------------
    use cmn_size, only: IPARW,IPAR,JPAR,LPAR,LPARW, LWEPAR, IDGRD
    use cmn_ctm, only: ETAA, ETAB, XDGRD, YDGRD, XDEDG, YDEDG, &
         AREAXY, PLAND, AIR
    use cmn_chem, only: NEMLIT, LITSRC, LITFAC
    use cmn_met, only: CWETE, CENTU, ZOFLE, T, PRECCNV, P, CLDFR
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input / Output
    integer, intent(in) :: NDAY, NDAYI
    logical, intent(in) :: LNEW_MONTH
    real(r8),  intent(in) :: dt_met

    !// Locals
    integer :: I, J, L
    real(r8) :: CTH      !// Cloud top height (top of convective mass flux)
    real(r8) :: UPV      !// Upward velocity [m/s]
    real(r8) :: fPR      !// Instant Price&Rind flash rate
    real(r8) :: fconv    !// Average cloud fraction

    !// Model flash rates over land and ocean. Their units do not
    !// matter, since they are scaled to match observations, but
    !// if you are curious the Price et al. equations are [fl/min].
    real(r8) :: flashL(IPAR,JPAR), flashO(IPAR,JPAR)

    !// We find the flash rates normalised to a climatological year.
    !// The sum of normFrate*dt is 1 for a year when it represents the
    !// climatological year. This is to distribute lightning emissions
    !// (e.g. 5Tg(N)) using normFrate.
    real(r8) :: normFrate(IPAR, JPAR)

    !// normFrate is the sum of flashL and flashO scaled by their
    !// respective factors, which also take care of the normalisation.
    real(r8) :: scaleOcean, scaleLand

    !// Vertical distribution of normFrate(i,j)
    real(r8) :: vertprof(LWEPAR)

    !// Meteorological variables
    real(r8) :: LNCWET(LWEPAR), AIRL(LWEPAR), TEMP(LWEPAR)
    real(r8) :: ZL(LWEPAR+1)

    !// Diagnostics to get total emissions
    real(r8) :: CNV_SUM, ntgyr
    real(r8) :: NO_SUM(IPAR, JPAR)
    real(r8) :: fLand, fOcean !// Total unscaled flash rates

    !// Counters & indices
    integer :: III, JJJ, IOS, ifnr, NLCL, NLCO !// Counters & indices
    integer :: LFREE, CBL, CTL                 !// Level indices



    !// Observed lightning flash rate (1/s), explained further down
    !// in this subroutine.
    real(r8), parameter :: obsFocean=11._r8, obsFland=35._r8
    real(r8), parameter :: obsFall = obsFocean + obsFland

    !// Scaling parameters depending on meteorological data
    !// Should be specified elsewhere at a later stage.
    !// The units of UCI calculated scalings were actually min/s,
    !// because the original flash rate equations are 1/min rather
    !// than the assumed 1/s. I decided to make the units more
    !// consistent in November 2014.
    real(r8), parameter :: z60 = 1._r8/60._r8
    !// T42L60 cycle 36r1 (1997-2010)
    !real(r8), parameter :: scaleOceanT42 = 3.5_r8 * 60._r8 / 1.45e9_r8
    !real(r8), parameter :: scaleLandT42  = 0.10_r8 * 60._r8 / 1.45e9_r8
    real(r8), parameter :: scaleOceanT42 = 1.286093e-7_r8
    real(r8), parameter :: scaleLandT42  = 3.795549e-9_r8

    !// T159L60 cycle 36r1
    real(r8), parameter :: scaleOceanT159 = 0.5439_r8 * 60._r8 / 1.45e9_r8
    real(r8), parameter :: scaleLandT159  = 0.0159_r8 * 60._r8 / 1.45e9_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'LIGHTNING_GMD2012'
    !// --------------------------------------------------------------------

    !// If no lightning, return
    if (NEMLIT .eq. 0)  return

    !// Initialize
    CNV_SUM = 0._r8
    NO_SUM(:,:) = 0._r8

    !// Get correct scaling factors
    if (IPARW .eq. 128) then
       scaleOcean = scaleOceanT42
       scaleLand  = scaleLandT42
       if (IDGRD /= 1) then
          print*,'lightning.f90: IDGRD is unknown:',idgrd
          stop 'in LIGHTNING_GMD2012'
       end if
    else if (IPARW .eq. 320) then
       scaleOcean = scaleOceanT159
       scaleLand  = scaleLandT159
       if (IDGRD /= 1) then
          print*,'lightning.f90: IDGRD is unknown:',idgrd
          stop 'in LIGHTNING_GMD2012'
       end if
    else
       print*,'* lightning.f90: scaleOcean and scaleLand not '// &
              'defined for current horizontal resolution'
       stop 'in LIGHTNING_GMD2012'
    end if

    !// Flash rates
    NormFrate(:,:) = 0._r8
    FlashL(:,:) = 0._r8
    FlashO(:,:) = 0._r8

    !// Lightning source
    LITSRC(:,:,:) = 0._r8

    !=================================================================
    ! Read lightning CDF from Ott et al [JGR, 2010]. (ltm, 1/25/11)
    !=================================================================
    if ( LNEW_MONTH .and. NDAY.eq.NDAYI ) then

       !// Initialize lightning profiles
       LITPROFILE(:,:) = 0._r8
       totlitfrq = 0._r8
       totlitn = 0._r8
       accFocean = 0._r8
       accFland  = 0._r8
       accDT     = 0._r8

       !// Echo info
       write(*,*) '   - INIT_LIGHTNING: Reading '//TRIM( INFILE_LIGHTNING )

       !// Get unused file ID
       IFNR = get_free_fileid()

       ! Open file containing lightning PDF data
       open( IFNR, file=trim( INFILE_LIGHTNING ), status='OLD', iostat=IOS)
       if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:1'
         
       ! Read 12 header lines
       do III = 1, 12
          read( IFNR, '(a)', IOSTAT=IOS ) 
          if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:2'
       end do

       ! Read NNLIGHT types of lightning profiles
       do III = 1, NNLIGHT
          read( IFNR,*,IOSTAT=IOS) (LITPROFILE(III,JJJ),JJJ=1,NLTYPE)
       end do

       ! Close file
       close( IFNR )

       !// Flash files
       if (getFlashFiles) then
          open(IFNR,file='flashrate.dta',form='unformatted')
          write(IFNR) IPAR,JPAR
          write(IFNR) XDGRD,YDGRD, AREAXY, PLAND, XDEDG,YDEDG
          close(IFNR)
       end if

    end if !// if ( LNEW_MONTH .and. NDAY.eq.NDAYI ) then


    !// Layer for free troposphere 1500m/~850hPa Oslo CTM3 uses LPAR=60
    LFREE = 12 - (LPARW-LPAR)


!$omp parallel private (I,J,L) &
!$omp         private (CTL,CTH) &
!$omp         private (LNCWET,TEMP,AIRL,ZL) &
!$omp         private (fconv,fPR,upv,vertprof) &
!$omp         shared (CWETE,CLDFR,T,ZOFLE,AIR,LITSRC) &
!$omp         shared (LFREE,PLAND,LITFAC,NormFrate,NO_SUM) &
!$omp         shared (scaleLand, scaleOcean) &
!$omp         shared (flashO, flashL) &
!$omp         default(NONE)
!$omp do
    do J = 1,JPAR
       do I = 1,IPAR
          !---------------------------------------------------------
          ! Flashrate: Horizontal and temporal distribution
          !---------------------------------------------------------

          !// Store updraft convection.
          !// I do not know the reason behind this parameterisation yet.
          !// Find convection at 850hPa. For all layers above this, we
          !// do not allow convective mass flux to be larger.
          do L = 1, LWEPAR
             LNCWET(L) = CWETE(I,J,L)
          end do

          !// Locate the top layer with convection
          !// Actually CTL is the layer above the uppermost layer having
          !// convection into it.
          CTL = 0
          do L = 1, LWEPAR
             if (LNCWET(L) .gt. 0._r8) then
                CTL = L+1 !// This is one too many
             end if
          end do
          !// Override due to error in CTL
          if (CTL .gt. LWEPAR) CTL = LWEPAR

          !// Check if there are any entraining updrafts;
          !// if none, then cycle to next box.
          if (CTL .eq. 0) cycle


          !// Retrieve temperature as 1D array
          do L = 1, LWEPAR
             TEMP(L) = T(I,J,L)
          end do


          !// Make sure the cloud top is colder than -40C
          !// and the surface is warmer than 0C. This guarantees
          !// that the cloud contains mixed liquid and ice phases.
          if ( TEMP(1) .lt. 273._r8 ) cycle
          if ( TEMP(CTL) .gt. 233._r8 ) cycle

          !// Get air mass [kg]
          do L = 1, LWEPAR
             AIRL(L) = AIR(I,J,L)
          end do
          !// Get height [km]
          do L = 1, LWEPAR+1
             ZL(L) = 1.e-3_r8 * (ZOFLE(L,I,J) - ZOFLE(1,I,J))
          end do


          !// Require updraft velocity > 0.01 m/s at layer 9,
          !// approximately 1500 m agl or 850 hPa.
          !// This eliminates instances of free tropospheric convection
          UPV = LNCWET(LFREE) / AIRL(LFREE)*(ZL(LFREE+1) - ZL(LFREE))
          if ( UPV .lt. 0.01e-3_r8 ) cycle !// unit is [km/s]


          !// Cloud top height above ground is top of level CTL, km
          CTH = ZL(CTL+1)
          !// Fraction of box area experiencing convection and lightning.
          !// use the average cloud fraction, weighted by convective flux
          fconv = sum( CLDFR(I,J,1:CTL) * LNCWET(1:CTL) ) &
                  / sum( LNCWET(1:CTL) )

          !// Use Price and Rind formulas for land and ocean
          !// flash frequency.
          !// Use additional scale factors for land and ocean to match
          !// the OTD-LIS climatology for flashes over land (35/s) 
          !// and ocean (11/s)
          !// Price and Rind equations have units 1/minutes, which we
          !// convert to 1/s.
          !// The factor fconv converts to flashrate averaged over box.
          if (PLAND(I,J) .gt. 0.25_r8) then
             fPR = 3.44e-5_r8 * CTH**4.9_r8 * z60 !// 1/min -> 1/s
             normFrate(I,J) = fPR * fconv * scaleLand !// 1/s
             flashL(I,J) = fPR * fconv !// Save unscaled oceanic
          else
             fPR = 6.4e-4_r8 * CTH**1.73_r8 * z60 !// 1/min -> 1/s
             normFrate(I,J) = fPR * fconv * scaleOcean !// 1/s
             flashO(I,J) = fPR * fconv !// Save unscaled oceanic
          end if

          !---------------------------------------------------------
          ! NOx emission, vertical distribution
          !---------------------------------------------------------

          !// If there's lightning in the column
          if ( normFrate(I,J) .gt. 0._r8 ) then

             !// Partition NOx emissions vertically using the PDFs from
             !// Ott et al. (2010) and regional definitions from
             !// Allen et al. (2010)
             call lightdist(I, J, CTL, ZL, normFrate(I,J), vertprof)

             !// Assign emissions into array for later use. 
             !// (units: fraction of annual emission per second)
             do L=1,CTL
                LITSRC(L,I,J) = vertprof(L)
                !// Sum up NO to test 5 Tg(N)/year
                NO_SUM(I,J) = NO_SUM(I,J) + vertprof(L)*LITFAC(1)
             end do

          end if

       end do !// do I = 1, IPAR
    end do !// do J = 1, JPAR
!$omp end do
!$omp end parallel

    !// Total unscaled land and ocean flashes
    fLand  = sum(flashL)
    fOcean = sum(flashO)

    !// Flash files
    if (getFlashFiles) then
       IFNR = get_free_fileid()
       open(IFNR,file='flashrate.dta',form='unformatted', position='append')
       write(IFNR) real(flashL, r4), real(flashO, r4)
       close(IFNR)
    end if

    !// How many grid boxes experienced lightning?
    NLCL = 0
    NLCO = 0
    do J = 1, JPAR
       do I = 1, IPAR
          if (flashL(I,J) .gt. 0._r8) NLCL = NLCL + 1
          if (flashO(I,J) .gt. 0._r8) NLCO = NLCO + 1
       end do
    end do
    write(*,'(A,3i7)') 'Lightning: # grid boxes w/flash (L/O/L+O): ', &
         NLCL, NLCO, NLCL+NLCO

    !// Total global fraction of annual total flashrate
    cnv_sum = sum(normFrate)
    !// Accumulate this fraction
    totlitfrq = totlitfrq + cnv_sum

    !// Sum up L-NOx until next update of LITSRC (duration dt_met)
    totlitn = totlitn + sum(no_sum) * dt_met * 1.e-9_r8 * 14._r8/30._r8
    !// Emissions this step, scaled to Tg(N)/yr
    ntgyr = sum(NO_SUM) * 1.e-9_r8 * 14._r8/30._r8 * dtyear

    !// Accumulate flashes during time step and total time.
    accFland = accFland + fLand * dt_met
    accFocean = accFocean + fOcean * dt_met
    accDT = accDT + dt_met

    !// Write total flash rate for ocean and land
    !  write(*,'(A,f12.3,1x,f12.3)')
    ! &     'Lightning: Flashrates (1/s) ocean/land:',
    ! &     fOcean*scaleOcean, fLand*scaleLand

    !// Fraction of annual total flashrates
    write(*,'(A,es23.17)') 'Lightning: Fraction of annual flashrate: ', cnv_sum

    !// Emissions during this time step (as Tg(N)/yr) and
    !// also the accumulated amount
    write(*,'(A,2f12.5)') 'Lightning: Tg(N)/yr & accumulated Tg(N):', &
         ntgyr, totlitn

    !// Accumualted flash frequency fraction of annual total, both
    !// new and old values. Their difference is cnv_sum.
    write(*,'(A,es23.17,1x,es23.17)') &
         'Lightning: Accu.flashfreq n/o: ', &
         totlitfrq,totlitfrq-cnv_sum

    !// Scaling factors. These are the values you need to modify
    !// scaleOean and scaleLand if you use other meteorological data.
    !// IMPORTANT:
    !// Remember that you need to use at least an annual mean, but
    !// preferably a mean over several years of data.
    !//
    !// Note also that for a single year of meteorological data,
    !// this printout may differ slightly from the applied scaling
    !// factors. This does not mean that the scalings are wrong,
    !// it may be that the applied scalings were calculated using
    !// a different year or several years.
    !//
    !// Printout is done after calculating lightning for the 
    !// last meteorological time step of each day, giving the
    !// average scalings for days from NDAYI through NDAY.
    !//
    !// Will find scaling so that
    !//   annual_modelflashes * scale = 1.
    !// separated into land and ocean:
    !//   scaleO = obsFocean/obsFall / annual_model_Oflashes
    !//   scaleL = obsFland/obsFall  / annual_model_Lflashes
    if (mod(accDT,86400._r8) .eq. 0._r8) &
         write(*,'(A,i4,es16.7,1x,es16.7)') &
         'Lightning: Accu.scal. NDAY/ocean/land: ',NDAY, &
         obsFocean/obsFall / (accFocean/accDT * dtyear), &
         obsFland/obsFall / (accFland/accDT * dtyear)

    !// Sensiblitity check
    if (ntgyr .lt. 1._r8 .or. ntgyr .gt. 12._r8) then
       print*,'*** WARNING: lightning.f90: Low/High LNOx: '// &
              'Check scaleLand and scaleOcean!'
    end if

    !// --------------------------------------------------------------------
  end subroutine LIGHTNING_GMD2012
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  SUBROUTINE LIGHTDIST( I, J, LTOP, ZL, TOTAL, VERTPROF )
    !// --------------------------------------------------------------------
    !// This subroutine was made by the
    !// Harvard University Atmospheric Chemistry Modeling Group
    !// and reads data used to partition the column lightning NOx
    !// into the CTM vertical layers.
    !//
    !// Code was rewritten to f90-format.
    !// Ole Amund Sovde, February 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, LWEPAR
    use cmn_ctm, only: JMON, YGRD, PLAND
    use cmn_parameters, only: ZPI180
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    INTEGER, INTENT(IN) :: I            ! Longitude index
    INTEGER, INTENT(IN) :: J            ! Latitude index 
    INTEGER, INTENT(IN) :: LTOP         ! Level of conv cloud top
    REAL(R8),  INTENT(IN) :: ZL(LWEPAR+1) ! Grid box bottom height above sfc [km]
    REAL(R8),  INTENT(IN) :: TOTAL        ! Column Total # of LNOx molec 

    !// Output
    REAL(R8),  INTENT(OUT) :: VERTPROF(LWEPAR) ! Vertical profile of LNOx
    !// --------------------------------------------------------------

    !  References:
    !  =======================================================================
    !  (1 ) Pickering et al., JGR 103, 31,203 - 31,316, 1998.
    !  (2 ) Ott et al., JGR, 2010
    !  (3 ) Allen et al., JGR, 2010
    ! 
    ! !REVISION HISTORY: 
    !  18 Sep 2002 - M. Evans - Initial version (based on Yuhang Wang's code)
    !  (1 ) Use functions IS_LAND and IS_WATER to determine if the given grid
    !        box is over land or water.  These functions work for all DAO met
    !        field data sets. (bmy, 4/2/02)
    !  (2 ) Renamed M2 to LTOP and THEIGHT to H0 for consistency w/ variable
    !        names w/in "lightning.f".  Now read the "light_dist.dat.geos3"
    !        file for  GEOS-3 directly from the DATA_DIR/lightning_NOx_200203/
    !        subdirectory.
    !        Now read the "light_dist.dat" file for GEOS-1, GEOS-STRAT directly 
    !        from the DATA_DIR/lightning_NOx_200203/ subdirectory.  Added 
    !        descriptive comment header.  Now trap I/O errors across all 
    !        platforms with subroutine "ioerror.f".  Updated comments, cosmetic 
    !        changes.  Redimension FRAC(NNLIGHT) to FRAC(LLPAR). (bmy, 4/2/02)
    !  (3 ) Deleted obsolete code from April 2002.  Now reference IU_FILE and
    !        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as the
    !        file unit number. (bmy, 6/27/02)
    !  (4 ) Now reference BXHEIGHT from "dao_mod.f" (bmy, 9/18/02)
    !  (5 ) Bug fix: add GEOS_4 to the #if block (bmy, 3/4/04)
    !  (6 ) Now bundled into "lightning.f".  CDF's are now read w/in
    !        routine INIT_LIGHTNING to allow parallelization (bmy, 4/14/04)
    !  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
    !  (8 ) Now uses near-land formulation (ltm, bmy, 5/10/06)
    !  (9 ) Added extra safety check for pathological boxes (bmy, 12/11/06)
    !  (10) Remove the near-land formulation, except for PRECON
    !       (ltm, bmy, 9/24/07)
    !  (11) Now use the Ott et al. [2010] profiles, and apply consistently with
    !        GMI model [Allen et al., 2010] (ltm, bmy, 1/25/11).
    !  10 Nov 2010 - R. Yantosca - Added ProTeX headers
    !-------------------------------------------------------------------------
    !LOCAL VARIABLES:
    INTEGER            :: M, MTYPE, L, MONTH,NNLEV
    REAL(R8)             :: ZHEIGHT, H0, YMID
    REAL(R8)             :: FRAC(LWEPAR)

    !=================================================================
    ! LIGHTDIST begins here!
    !=================================================================

    ! Initialize 
    MTYPE    = 0
    VERTPROF(:) = 0_r8
    MONTH    = JMON
    YMID     = YGRD(J) * ZPI180

    !=================================================================
    ! Test whether location (I,J) is continental, marine, or snow/ice
    !
    ! Depending on the combination of land/water and latitude, 
    ! assign a flag describing the type of lightning:
    !
    !   MTYPE = 1: ocean lightning
    !   MTYPE = 2: tropical continental lightning
    !   MTYPE = 3: midlatitude continental lightning 
    !   MTYPE = 4: subtropical lightning
    !             
    ! (ltm, bmy, 1/25/11)
    !=================================================================

    ! Assign profile kind to grid box, following Allen et al. 
    ! [JGR, 2010] (ltm, 1/25,11)
    SELECT CASE (MONTH)

    ! Southern Hemisphere Summer
    CASE ( 1,2,3,12 )

       IF ( ABS(YMID) .le. 15._r8 ) THEN
          IF ( PLAND(I,J) .gt. 0.25_r8 ) THEN
             MTYPE = 2        ! Tropical continental
          ELSE
             MTYPE = 1        ! Tropical marine
          ENDIF
       ELSE IF ( ( YMID .gt. 15._r8 ) .and. ( YMID .le. 30._r8 ) ) THEN
          MTYPE = 4           ! N. Subtropics
       ELSE IF ( ( YMID .ge. -40._r8 ) .and. ( YMID .lt. -15._r8 ) ) THEN
          MTYPE = 4           ! S. Subtropics
       ELSE
          MTYPE = 3           ! Midlatitude
       END IF

    ! Equinox months
    CASE ( 4,5,10,11 )

       IF ( ABS(YMID) .le. 15._r8 ) THEN
          IF ( PLAND(I,J) .gt. 0.25_r8 ) THEN
             MTYPE = 2        ! Tropical continental
          ELSE
             MTYPE = 1        ! Tropical marine
          ENDIF
       ELSE IF ( ABS(YMID) .le. 30._r8 ) THEN
          MTYPE = 4           ! Subtropics
       ELSE
          MTYPE = 3           ! Midlatitude
       END IF

    ! Northern Hemisphere Summer
    CASE ( 6,7,8,9 )

       IF ( ABS(YMID) .le. 15._r8 ) THEN
          IF ( PLAND(I,J) .gt. 0.25_r8 ) THEN
             MTYPE = 2        ! Tropical continental
          ELSE
             MTYPE = 1        ! Tropical marine
          ENDIF
       ELSE IF ( ( YMID .gt. 15._r8 ) .and. ( YMID .le. 40._r8 ) ) THEN
          MTYPE = 4           ! N. Subtropics
       ELSE IF ( ( YMID .ge. -30._r8 ) .and. ( YMID .lt. -15._r8 ) ) THEN
          MTYPE = 4           ! S. Subtropics
       ELSE
          MTYPE = 3           ! Midlatitude
       END IF
         
    END SELECT

    ! Extra safety check for pathological grid boxes (bmy, 11/29/06)
    IF ( MTYPE .eq. 0 ) RETURN

    !=================================================================
    ! Use the CDF for this type of lightning to partition the total
    ! column lightning into the GEOS-3, GEOS-4, or GEOS-5 layers
    !=================================================================
    !ZHEIGHT = 0._r8
    ! Height of top of uppermost box for distribution
    H0 = ZL(LTOP+1)
    FRAC(:) = 0._r8

    ! Compute the height [km] at the top of each vertical level.
    ! Look up the cumulative fraction of NOx for each vertical level
    DO L = 1, LTOP
       !ZHEIGHT = ZL(L+1)
       !NNLEV = NINT( (ZHEIGHT/H0)*real(NNLIGHT, r8))
       NNLEV = NINT( (ZL(L+1) / H0) * real(NNLIGHT, r8))
       if (NNLEV .gt. NNLIGHT) then
          print'(4I5,3F10.3,i6)', i,j,l, ltop, ZL(L+1), h0,NNLEV
       end if
       FRAC(L) = LITPROFILE( NNLEV, MTYPE ) * 0.01
    END DO

    ! Convert from cumulative fraction to fraction for each level
    DO L = LTOP, 2, - 1
       FRAC(L) = FRAC(L) - FRAC(L-1)
    END DO
      
    ! Partition lightning NOx by layer into VERTPROF
    DO L = 1, LTOP
       VERTPROF(L) = ( FRAC(L) * TOTAL )
    END DO

    !// --------------------------------------------------------------------
  END SUBROUTINE LIGHTDIST
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine filterLNC(i, j, LFREE, flux, ent, CBL, CTL)
    !// --------------------------------------------------------------------
    !// Routine to filter out convective instances that should not
    !// produce lightning.
    !//
    !// Ole Amund Sovde, March 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LWEPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)   :: I, J, LFREE
    !// Inout/output
    real(r8), intent(inout) :: flux(LWEPAR), ent(LWEPAR)
    !// Output
    integer, intent(out)  :: CBL, CTL

    !// Locals
    integer :: L, NLV, L1, L2
    integer :: intTHICK(LWEPAR) !// Thickness of intervals
    integer :: intLSTRT(LWEPAR) !// Start indices of intervals
    !// --------------------------------------------------------------------

    !// Check if there are any updrafts with entrainment
    if (sum(flux) .gt. 0._r8) then
       !// There is convective mass flux. Skip lightning if there is
       !// no entrainment (but not if metdata do not contain entrainment).
       if (LentuExist .and. sum(ent).eq.0._r8) then
          CBL = 0
          CTL = 0
          return
       end if
    else
       !// There is no convective mass flux
       CBL = 0
       CTL = 0
       return
    end if

    !// Initialize
    NLV = 0   !// # of convective intervals separated by zero convection
    L1 = 0
    L2 = 0
    intTHICK(:) = 0 !// Thickness of intervals
    intLSTRT(:) = 0 !// Start indices of intervals

    !// Then check intervals
    do L = 1, LWEPAR
       if (flux(L) .gt. 0._r8 .and. L1 .eq. 0) then
          !// New interval starts here
          L1 = L
          L2 = L !// Also initialize end here
       else if (flux(L) .gt. 0._r8 .and. L1 /= 0) then
          L2 = L !// Interval still reaches this level
       end if

       !// May have L2.eq.LWEPAR
       if (flux(L) .eq. 0._r8 .or. L2 .eq. LWEPAR) then
          if (L2 .gt. 0 .and. L1 .gt. 0) then
             NLV = NLV + 1
             !// Thickness of this interval
             intTHICK(NLV) = L2 - L1 + 1
             !// Start level of this interval
             intLSTRT(NLV) = L1
             !// Reset
             L1 = 0
             L2 = 0
          end if
       end if
    end do

    if (NLV .eq. 0) then
       !// If only one interval or if there is no convection (which
       !// should already be filtered out) we set CBL and CTL directly:
       CBL = 0
       CTL = 0
       flux(:) = 0._r8
       ent(:)  = 0._r8
    else if (NLV .eq. 1) then
       !// If only one interval or if there is no convection (which
       !// should already be filtered out) we set CBL and CTL directly:
       CBL = intLSTRT(NLV)
       CTL = intLSTRT(NLV) + intTHICK(NLV) - 1
    else
       !// More than one interval with convective mass flux.

       !// First, if any of the intervals are not thicker than
       !// 2 model level, there will be no lightning.
       if (maxval(intTHICK) .le. 2) then
          CBL = 0
          CTL = 0
          flux(:) = 0._r8
          ent(:)  = 0._r8
       else
          !// Find largest interval
          L1 = 0
          L2 = 0
          do L = 1, NLV
             !// If two intervals have the same thickness, we must
             !// chose which should be used.
             !// It is probably best to use the lowermost interval;
             !// when this occurs, we should not have deep convection.
             !// Another possibility could be to use the uppermost
             !// interval if it is thicker than some limit, e.g. 5km.
             if (intTHICK(L) .gt. L1) then
                !// Thickness and start of thickest level
                L1 = intTHICK(L)
                L2 = intLSTRT(L)
             end if
          end do
          !// Zero up to thickest level
          do L = LFREE, L2-1
             flux(L) = 0._r8
             ent(L)  = 0._r8
          end do
          !// Zero above thickest level
          do L = L2+L1, LWEPAR
             flux(L) = 0._r8
             ent(L)  = 0._r8
          end do
          CBL = L2
          CTL = L2 + L1 - 1
       end if !// if (maxval(intTHICK) .lt. 2) then
    end if !// if (NLV .eq. 0) then

    !// --------------------------------------------------------------------
  end subroutine filterLNC
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine distland(landlit)
    !// --------------------------------------------------------------------
    !// Will find LANDLIT which is all grid boxes that are land
    !// or have land in 300km proximity.
    !//
    !// Ole Amund Sovde, March 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR,JPAR
    use cmn_ctm, only: YGRD, XGRD
    use cmn_sfc, only: LSMASK
    use cmn_parameters, only: A0
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Output
    real(r8), intent(out) :: landlit(IPAR,JPAR)
    !// Locals
    real(r8) :: a, c, dist, phi1, phi2, the1, the2
    integer :: ii,jj, i,j, foundL
    real(r8) :: land(ipar,jpar)

    real(r8), parameter :: landprox = 300.e3_r8 !// [m]
    !// --------------------------------------------------------------------

    !// Get land+lakes+islands from LSMASK
    land(:,:) = 0._r8
    do j = 1, jpar
       do i = 1, ipar
          if (lsmask(i,j,2)+lsmask(i,j,3)+lsmask(i,j,4) .gt. 0._r8) &
               land(i,j) = 1._r8
       end do
    end do


    !// Use lat/lon to calculate distance to all other grid boxes
    landlit(:,:) = 0
    do j = 1, jpar
       do i = 1, ipar
          phi1 = ygrd(j)
          the1 = xgrd(i)
          if (land(i,j) .gt. 0._r8) then
             !// not ocean: treat as land
             landlit(i,j) = 1._r8
             cycle
          end if

          foundL = 0
          do jj = 1, jpar
             do ii = 1, ipar
                !// Skip the box we are checking
                if (jj .eq. j .and. ii .eq. i) then
                   cycle !// go to next ii
                end if

                phi2 = ygrd(jj)
                the2 = xgrd(ii)

                !// get distance
                a = sin((phi2-phi1)*0.5_r8)**2 &
                      + cos(phi1)*cos(phi2)*sin((the2-the1)*0.5_r8)**2
                c = 2._r8*atan2(sqrt(a),sqrt(1._r8-a))
                dist = A0 * c 
                if (dist.lt.landprox .and. land(ii,jj) .gt. 0._r8) then
                   !// i,j has land grid box within range
                   landlit(i,j) = 1._r8
                   foundL = 1
                   exit !// exit this i,j
                end if
             end do
             if (foundL .eq. 1) exit !// exit this i,j
          end do

       end do
    end do
    !// --------------------------------------------------------------------
  end subroutine distland
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine LIGHTNING_ALLEN2002(NDAY,NDAYI,dt_met,LNEW_MONTH)
    !// --------------------------------------------------------------------
    ! Subroutine lightning calculates the 3-D lightning source
    ! for use in the SOURCE subroutine. It must be called for each met
    ! field after p-wind.f
    !
    ! Algorithm:
    !   MFLUX parameterization by Allen & Pickering (2002).
    !
    !//
    !// Ole Amund Sovde, March 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: IPARW,IPAR,JPAR,LPAR,LPARW, LWEPAR, IDGRD
    use cmn_ctm, only: ETAA, ETAB, XDGRD, YDGRD, XDEDG, YDEDG, &
         AREAXY, PLAND, AIR
    use cmn_chem, only: NEMLIT, LITSRC, LITFAC
    use cmn_met, only: CWETE, CENTU, ZOFLE, T, PRECCNV, P, CLDFR
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input / Output
    integer, intent(in) :: NDAY, NDAYI
    logical, intent(in) :: LNEW_MONTH
    real(r8),  intent(in) :: dt_met

    !// Locals
    integer :: I, J, L
    real(r8) :: DZ         !// Cloud top thickness frz-top
    real(r8) :: flfcg    !// Flash rate CG
    real(r8) :: fcg      !// Fraction of CG/(CG+CC)
    !// Model flash rates over land and ocean. Their units do not
    !// matter, since they are scaled to match observations, but
    !// if you are curious the Price et al. equations are [fl/min].
    real(r8) :: flashL(IPAR,JPAR), flashO(IPAR,JPAR)

    !// We find the flash rates normalised to a climatological year.
    !// The sum of normFrate*dt is 1 for a year when it represents the
    !// climatological year. This is to distribute lightning emissions
    !// (e.g. 5Tg(N)) using normFrate.
    real(r8) :: normFrate(IPAR, JPAR)

    !// normFrate is the sum of flashL and flashO scaled by their
    !// respective factors, which also take care of the normalisation.
    real(r8) :: scaleOcean, scaleLand

    !// Vertical distribution of normFrate(i,j)
    real(r8) :: vertprof(LWEPAR)

    !// Additional factors/variables to calculate flash rates
    real(r8) :: FACT1, FACT2, fs
    real(r8) :: AY(JPAR)
    integer :: J30

    !// Meteorological variables
    real(r8) :: LNCWET(LWEPAR), LNCENTU(LWEPAR), TEMP(LWEPAR)
    real(r8) :: ZL(LWEPAR+1)
    real(r8) :: PBOT !// Pressure at bottom of grid box
    real(r8) :: MFLX !// Mass flux
    real(r8) :: MENT !// Entrainment flux
    real(r8) :: minT, maxT !// min/max temperatures in cloud

    !// Diagnostics to get total emissions
    real(r8) :: CNV_SUM, ntgyr
    real(r8) :: NO_SUM(IPAR, JPAR)
    real(r8) :: fLand, fOcean !// Total unscaled flash rates

    !// Counters & indices
    integer :: III, JJJ, IOS, ifnr, NLCL, NLCO  !// Counters & indices
    integer :: LFREE, CBL, CTL, LFLX, LFRZ !// Level indices


    !// Pressures for defining shallow or deep convection.
    real(r8), parameter :: pconvection = 440._r8



    !// Observed lightning flash rated (1/s).
    !// These are taken from  LIS/OTD data LISOTD_HRAC_V2.3.2014.hdf,
    !// as produced by Cecil et al. (2014, J. Atmos. Res.,
    !// vol 135-136, 404-414, doi:10.1016/j.atmosres.2012.06.028).
    !// With an annual global flash rate of 46fl/s, using land
    !// fractions separates this into 9.1fl/s over ocean and
    !// 36.9 over land. If we were to define all grid boxes having
    !// land within 300km reach (roughly similar to Christian et al.
    !// (2003, JGR, 108(D1),doi:10.1029/2002JD002347), these numbers
    !// would be 4.25fl/s and 41.75fl/s, respectively.
    !// Christian et al. (2003) suggested about 5fl/s over ocean.
    !// However, we will use a stricter definition of land, so
    !// we stick to the former.
    real(r8), parameter :: obsFocean=9.1_r8, obsFland=36.9_r8
    real(r8), parameter :: obsFall = obsFocean + obsFland


    !// Scaling parameters depending on meteorological data.
    !// Ole Amund Sovde, February 2015
    !//
    !// Price etal equations have units 1/minutes, while observations
    !// are given in 1/s. This conversion is taken care of by the factors
    !// scaleLand and scaleOcean. Thus, after multiplying with the
    !// scaling factors, the units of flash rate is 1/s.
    !//
    !// 1997 - 2010 cy36 average scaling factor
    real(r8), parameter :: scaleOceanT42 = 7.265455e-13_r8
    real(r8), parameter :: scaleLandT42  = 3.407279e-12_r8

    !//
    !// T159L60 cycle 36r1 (T159 2007 / T42 2007 * clim.T42)
    real(r8), parameter :: scaleOceanT159 = 1.239369e-13_r8
    real(r8), parameter :: scaleLandT159  = 5.968714e-13_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'LIGHTNING_ALLEN2002'
    !// --------------------------------------------------------------------

    !// If no lightning, return
    if (NEMLIT .eq. 0)  return

    !// Initialize
    CNV_SUM = 0._r8
    NO_SUM(:,:)  = 0._r8

    !// Get correct scaling factors
    if (IPARW .eq. 128) then
       if (IDGRD .eq. 1) then
          scaleOcean = scaleOceanT42
          scaleLand  = scaleLandT42
       else
          print*,'lightning.f90: IDGRD is unknown:',idgrd
          stop 'in LIGHTNING2'
       end if
    else if (IPARW .eq. 320) then
       if (IDGRD .eq. 1) then
          scaleOcean = scaleOceanT159
          scaleLand  = scaleLandT159
       else
          print*,'p-lit_cdh.f: IDGRD is unknown:',idgrd
          stop 'in LIGHTNING2'
       end if
    else
       print*,'* p-lit_cdh.f: scaleOcean and scaleLand not '// &
            'defined for current horizontal resolution'
    end if

    !// Flash rates
    normFrate(:,:) = 0._r8
    FlashL(:,:) = 0._r8
    FlashO(:,:) = 0._r8

    !// Lightning source
    LITSRC(:,:,:) = 0._r8

    !=================================================================
    ! Read lightning CDF from Ott et al [JGR, 2010]. (ltm, 1/25/11)
    !=================================================================
    if ( LNEW_MONTH .and. NDAY.eq.NDAYI ) then

       !// Initialize lightning profiles
       LITPROFILE(:,:) = 0._r8
       totlitfrq = 0._r8
       totlitn = 0._r8
       accFocean = 0._r8
       accFland  = 0._r8
       accDT     = 0._r8

       !// Echo info
       write(*,*) '   - INIT_LIGHTNING: Reading '//TRIM( INFILE_LIGHTNING )

       !// Get unused file ID
       IFNR = get_free_fileid()

       ! Open file containing lightning PDF data
       OPEN( IFNR, FILE=TRIM( INFILE_LIGHTNING ), STATUS='OLD', IOSTAT=IOS)
       if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:1'
         
       ! Read 12 header lines
       DO III = 1, 12
          READ( IFNR, '(a)', IOSTAT=IOS ) 
          if ( IOS /= 0 ) stop f90file//':'//subr//': lightdist not read:2'
       END DO
         
       ! Read NNLIGHT types of lightning profiles
       DO III = 1, NNLIGHT
          READ( IFNR,*,IOSTAT=IOS) (LITPROFILE(III,JJJ),JJJ=1,NLTYPE)
       END DO

       ! Close file
       CLOSE( IFNR )

       !// Flash files
       if (getFlashFiles) then
          open(IFNR,file='flashrate.dta',form='unformatted')
          write(IFNR) IPAR,JPAR
          write(IFNR) XDGRD,YDGRD, AREAXY, PLAND, XDEDG,YDEDG
          close(IFNR)
       end if

       if (maxval(CENTU) .gt. 0._r8) then
          LentuExist = .true.
       else
          LentuExist = .false.
       end if

    end if

    ! Locate layer of approximately 1500m above surface, which is
    ! about layer 12 in native vertical L60. If collapsing this
    ! will be layer 9 (as in UCI standard code).
    LFREE = 12 - (LPARW-LPAR)

    !// Area of grid boxes
    AY(:) = AREAXY(1,:)
    !// Grid box for 30S or 30N
    J30  = 1
    do while (YDEDG(J30) .lt. -30._r8)
       J30 = J30 + 1
    end do



    !=================================================================
    ! Flashrate: Horizontal and temporal distribution
    !=================================================================
!$omp parallel private (I,J,L) &
!$omp         private (CBL,CTL,LFRZ,LFLX) &
!$omp         private (DZ,LNCWET,LNCENTU,ZL,TEMP) &
!$omp         private (PBOT, MFLX, MENT) &
!$omp         private (fcg,flfcg) &
!$omp         private (FACT1,FACT2,maxT,minT,fs) &
!$omp         private (vertprof) &
!$omp         shared (CWETE,CENTU,PRECCNV,T,ZOFLE,AREAXY) &
!$omp         shared (P,ETAA,ETAB) &
!$omp         shared (LITFAC,LITSRC,normFrate,NO_SUM) &
!$omp         shared (LentuExist, LFREE, PLAND, J30, AY) &
!$omp         shared (scaleLand, scaleOcean) &
!$omp         shared (flashO, flashL) &
!$omp         default(NONE)
!$omp do
    do J = 1,JPAR
       do I = 1,IPAR

          !// Store updraft convection.
          LNCWET(1:(LFREE-1)) = 0._r8
          LNCENTU(1:(LFREE-1)) = 0._r8
          do L = LFREE, LWEPAR
             LNCWET(L) = CWETE(I,J,L)
             LNCENTU(L) = CENTU(I,J,L)
          end do

          !// Use rain as screening. There may be some convection
          !// within large scale rain, without producing convective
          !// rain falling to the ground. It is not clear to what extent
          !// that can happen for vigorous convection. Therefore I
          !// check only convective rain at the surface; it seems to
          !// be a filter that works well to get temporal distribution
          !// fairly correct.
          if (PRECCNV(I,J,1) .eq. 0._r8) cycle

          !// There may be some instances where convective mass flux
          !// occur in single layers or separated into several intervals.
          !// This is due to the nature of the IFS model physics.
          !// If e.g. level 10 and 20 are the only layers with
          !// convective mass flux, there should be no lightning.
          call filterLNC(I, J, LFREE, LNCWET, LNCENTU, CBL, CTL)


          !// Check if there are any entraining updrafts;
          !// if none, then cycle to next box.
          if (CTL .eq. 0) cycle


          !// Restrict lightning to be produced only when convection
          !// reaches a height where pressure is less than some
          !// limit (pconvection).
          LFLX = 0
          do L = CBL, CTL
             PBOT = ETAA(L) + P(I,J)*ETAB(L)
             if (LentuExist) then
                MENT = sum(LNCENTU(L:CTL))
                if (PBOT .le. pconvection .and. MENT .gt. 0._r8) then
                   LFLX = L
                   exit           !// Exit L-loop
                end if
             else
                if (PBOT .le. pconvection) then
                   LFLX = L
                   exit           !// Exit L-loop
                end if
             end if
          end do

          !// Skip lightning if convection does not reach the limit.
          if (LFLX .eq. 0) cycle



          !// Retrieve temperature as 1D array
          do L = 1, LWEPAR
             TEMP(L) = T(I,J,L)
          end do

          !// Make sure the cloud is mixed phase, containing both
          !// liquid and ice. We assume this occurs when minimum
          !// temperature in the column is colder than -40C and
          !// the maximum temperature is warmer than 0C.
          minT = minval(TEMP(LFREE:CTL))
          maxT = maxval(TEMP(LFREE:CTL))

          if (maxT .lt. 273.15_r8) cycle
          if (minT .gt. 233.15_r8) cycle


          !// Find freezing level, i.e. highest level where T>0C
          LFRZ = 0
          do L = LFREE, LWEPAR
             if (TEMP(L) .gt. 273.15_r8) then
                LFRZ = L
             else
                exit
             end if
          end do

          !// If not below zero, skip convection
          if (LFRZ .eq. 0) cycle


          !// Find layer bottom height above ground [km]
          do L = 1, LWEPAR+1
             ZL(L) = 1.e-3_r8 * (ZOFLE(L,I,J) - ZOFLE(1,I,J))
          end do

          !// Thickness cloud top vs center of LFRZ [km]
          DZ = (ZL(CTL+1) - 0.5_r8 * (ZL(LFRZ) + ZL(LFRZ+1)))
          DZ = max(5.5_r8, DZ)
          if (DZ .lt. 14._r8) then
             fcg = 1._r8 / (0.021_r8*DZ**4 -0.648_r8*DZ**3 + 7.493_r8*DZ**2 &
                           -36.54_r8*DZ + 63.09_r8 + 1._r8)
          else
             fcg = 0.02_r8
          end if


          !// MFLX [kg/m2/min]
          MFLX = min(LNCWET(LFLX) / AY(J) * 60._r8, 9.6_r8)


          fs = 1._r8 !FACT1 * FACT2

          if (PLAND(I,J) .gt. 0.25_r8) then
             flfcg = AY(J)/AY(J30) &
                  * (-2.78e-2_r8 + 3.33e-1_r8*MFLX - 6.93e-1_r8*MFLX**2 &
                     + 5.22e-1_r8*MFLX**3 - 3.74e-2_r8*MFLX**4)


             normFrate(I,J) = flfcg / fcg * fs * scaleLand
             flashL(I,J) = flfcg / fcg * fs !// Save unscaled oceanic
          else
             flfcg = AY(J)/AY(J30) &
                  * (6.41e-3_r8 + 5.88e-2_r8*MFLX - 3.26e-1_r8*MFLX**2 &
                     + 2.85e-1_r8*MFLX**3 - 7.41e-3_r8*MFLX**4)

             normFrate(I,J) = flfcg / fcg * fs * scaleOcean
             flashO(I,J) = flfcg / fcg * fs !// Save unscaled oceanic
          end if



          !---------------------------------------------------------
          ! NOx emission, vertical distribution
          !---------------------------------------------------------

          ! If there's lightning in the column
          if ( normFrate(I,J) .gt. 0._r8 ) then

             ! Partition NOx emissions vertically using the PDFs from
             ! Ott et al. (2010) and regional definitions from
             ! Allen et al. (2010)
             call lightdist( I, J, CTL, ZL, normFrate(I,J), vertprof )

             ! Assign emissions into array for later use. 
             ! (units: fraction of annual emission per second)
             do L=1,CTL
                LITSRC(L,I,J) = vertprof(L)
                !// Sum up NO to test 5 Tg(N)/year
                NO_SUM(I,J) = NO_SUM(I,J) + (vertprof(L)*LITFAC(1))
             end do

          end if

       end do
    end do !// do J = 1,JPAR
!$omp end do
!$omp end parallel

    !=================================================================
    ! Diagnostics for lightning production
    !=================================================================
    !// Total unscaled land and ocean flashes
    fLand  = sum(flashL)
    fOcean = sum(flashO)

    !// Flash files
    if (getFlashFiles) then
       IFNR = get_free_fileid()
       open(IFNR,file='flashrate.dta',form='unformatted', position='append')
       write(IFNR) real(flashL, r4), real(flashO, r4)
       close(IFNR)
    end if

    !// How many grid boxes experienced lightning?
    NLCL = 0
    NLCO = 0
    do J = 1, JPAR
       do I = 1, IPAR
          if (flashL(I,J) .gt. 0._r8) NLCL = NLCL + 1
          if (flashO(I,J) .gt. 0._r8) NLCO = NLCO + 1
       end do
    end do
    write(*,'(A,3i7)') &
         'Lightning: # grid boxes w/flash (L/O/tot): ',&
         NLCL,NLCO,NLCL+NLCO

    !// Total global fraction of annual total flashrate. Must
    !// multiply by sum up to make annual total 1
    cnv_sum = sum(normFrate) * dt_met

    !// Accumulate this fraction
    totlitfrq = totlitfrq + cnv_sum

    !// Sum up L-NOx until next update of LITSRC (duration dt_met)
    totlitn = totlitn + sum(no_sum) * dt_met * 1.e-9_r8 * 14._r8/30._r8
    !// Emissions this step, scaled to Tg(N)/yr
    ntgyr = sum(NO_SUM) * 1.e-9_r8 * 14._r8/30._r8 * dtyear

    !// Accumulate flashes during time step and total time.
    accFland = accFland + fLand * dt_met
    accFocean = accFocean + fOcean * dt_met
    accDT = accDT + dt_met

    !// Write total flash rate for ocean and land
    !  write(*,'(A,f12.3,1x,f12.3)')
    ! &     'Lightning: Flashrates (1/s) ocean/land:',
    ! &     fOcean*scaleOcean, fLand*scaleLand

    !// Fraction of annual total flashrates
    write(*,'(A,es23.17)') 'Lightning: Fraction of annual flashes: ', cnv_sum

    !// Emissions during this time step (as Tg(N)/yr) and
    !// also the accumulated amount
    write(*,'(A,2f12.5)') 'Lightning: Tg(N)/yr & accumulated Tg(N):', &
         ntgyr, totlitn

    !// Accumualted flash frequency fraction of annual total, both
    !// new and old values. Their difference is cnv_sum.
    write(*,'(A,es23.17,1x,es23.17)') &
         'Lightning: Accu.flashfreq n/o: ', &
         totlitfrq, totlitfrq - cnv_sum

    !// Scaling factors. These are the values you need to modify
    !// scaleOean and scaleLand if you use other meteorological data.
    !// IMPORTANT:
    !// Remember that you need to use at least an annual mean, but
    !// preferably a mean over several years of data.
    !//
    !// Note also that for a single year of meteorological data,
    !// this printout may differ slightly from the applied scaling
    !// factors. This does not mean that the scalings are wrong,
    !// it may be that the applied scalings were calculated using
    !// a different year or several years.
    !//
    !// Printout is done after calculating lightning for the 
    !// last meteorological time step of each day, giving the
    !// average scalings for days from NDAYI through NDAY.
    !//
    !// Will find scaling so that
    !//   annual_modelflashes * scale = 1.
    !// separated into land and ocean:
    !//   scaleO = obsFocean/obsFall / annual_model_Oflashes
    !//   scaleL = obsFland/obsFall  / annual_model_Lflashes
    if (mod(accDT,86400._r8) .eq. 0._r8) &
         write(*,'(A,i4,es16.7,1x,es16.7)') &
         'Lightning: Accu.scal. NDAY/ocean/land: ',NDAY, &
         obsFocean/obsFall / (accFocean/accDT * dtyear), &
         obsFland/obsFall / (accFland/accDT * dtyear)


    !// Sensiblitity check
    if (ntgyr .lt. 1._r8 .or. ntgyr .gt. 12._r8) then
       print*,'*** WARNING: lightning.f90: Low/High LNOx: '// &
              'Check scaleLand and scaleOcean!'
    end if


    !// --------------------------------------------------------------------
  end subroutine LIGHTNING_ALLEN2002
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
end module lightning
!//=========================================================================
