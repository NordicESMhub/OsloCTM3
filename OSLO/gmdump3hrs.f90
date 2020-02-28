!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Dump instantaneous 3D data every 3 hours.
!//=========================================================================
module gmdump3hrs
  !// ----------------------------------------------------------------------
  !// MODULE: gmdump3hrs
  !// DESCRIPTION:  Puts out selected tracer components in 3D each hour.
  !//               Based on old dump_3hrs.F from CTM2.
  !//
  !// Ole Amund Sovde, June 2012
  !// ----------------------------------------------------------------------
  use cmn_size, only: NPAR, NOTRPAR
  use cmn_ctm, only: LDUMP3HRS
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// List of tracers to put out
  integer, parameter :: trp_nr = 7, sul_nr = 3, slt_nr = 8, min_nr = 8, &
       nit_nr = 5, bio_nr = 4, moa_nr=2, ffc_nr = 4, bfc_nr = 4, &
       soa_ant_nr = 17, soa_nat_nr = 4, ntr_nr = 1
  integer, parameter,dimension(trp_nr) :: &
       trp_list = (/ 1, 6, 46, 41, 42, 43, 44 /)
  integer, parameter,dimension(sul_nr) :: &
       sul_list = (/ 71, 72, 73 /)
  integer, parameter,dimension(slt_nr) :: &
       slt_list = (/ 201, 202, 203, 204, 205, 206, 207, 208 /)
  integer, parameter,dimension(min_nr) :: &
       min_list = (/ 211, 212, 213, 214, 215, 216, 217, 218 /)
  integer, parameter,dimension(nit_nr) :: &
       snn_list = (/ 62,  63,  64,  65,  73 /)
  integer, parameter,dimension(bio_nr) :: &
       bio_list = (/ 230, 231, 240, 241 /)
  integer, parameter,dimension(moa_nr) :: &
       moa_list = (/ 236, 237 /)
  integer, parameter,dimension(ffc_nr) :: &
       ffc_list = (/ 232, 233, 242, 243 /)
  integer, parameter,dimension(bfc_nr) :: &
       bfc_list = (/ 234, 235, 244, 245 /)
  ! NOTRPAR is of same size as the length of list in tracer_list
  integer, parameter,dimension(ntr_nr) :: &
       ntr_list = (/ 40 /)

  !// Not included yet (usually put out as the sum, not every species)
  integer, parameter,dimension(soa_ant_nr) :: &
       soa_ant_list = (/ 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, &
                     175, 176, 177, 178, 179, 182, 183 /)
  integer, parameter,dimension(soa_nat_nr) :: &
       soa_nat_list = (/ 188, 189, 190, 191 /)

  character(len=10),dimension(NPAR) :: gmname
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'gmdump3hrs.f90'
  !// ----------------------------------------------------------------------
  private
  public dump3hrs
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine dump3hrs(NDAYI,NDAY,NMET)
    !// --------------------------------------------------------------------
    !// Dump data to files.
    !// The routine gm_dump_nc will divide by volume, so whatever is sent
    !// to that routine will be per volume.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR,JPAR, LPAR, LOSLOCTROP, LSULPHUR, LSALT, &
         LDUST, LNITRATE, LBCOC, LSOA, TRACER_ID_MAX
    use cmn_ctm, only: STT, AIR, JYEAR, JMON, JDATE
    use cmn_chem, only: TNAME
    use cmn_oslo, only: trsp_idx, Xtrsp_idx, XSTT, XTNAME
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in) :: NDAYI, NDAY, NMET

    !// Locals
    real(r8) :: TMPFIELD(IPAR,JPAR,LPAR,2)
    character(len=3)  :: cday                            !Char of day of year
    character(len=4)  :: cyear                           !Char of year
    character(len=14) :: ncfilename                      !Output filename
    character(len=80) :: time_label                      !Time label or unit
    character(len=10) :: tmpname
    integer :: N,sl

    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'dump3hrs'
    !// --------------------------------------------------------------------

    !// Check whether to do dump or not
    if (.not.LDUMP3HRS(1)) return

    !// Initialize
    If ((NDAY.eq.NDAYI).AND.(NMET.eq.1)) Then
       !// Initialize GMNAME (only 8 characters)
       GMNAME(:) = TNAME(:)
       !do N = 1, NPAR
       !   tmpname = GMNAME(N)
       !   sl = strlen(tmpname)
       !   tmpname = '__________' !// Trailings: '_'
       !   tmpname(1:sl) = GMNAME(N)
       !   GMNAME(N) = tmpname(1:10)
       !end do

       !// Override names
       if (trsp_idx(230).gt.0) GMNAME(trsp_idx(230)) = 'OMPHOB  '
       if (trsp_idx(231).gt.0) GMNAME(trsp_idx(231)) = 'OMPHIL  '
       if (trsp_idx(232).gt.0) GMNAME(trsp_idx(232)) = 'OMPHOB  '
       if (trsp_idx(233).gt.0) GMNAME(trsp_idx(233)) = 'OMPHIL  '
       if (trsp_idx(234).gt.0) GMNAME(trsp_idx(234)) = 'OMPHOB  '
       if (trsp_idx(235).gt.0) GMNAME(trsp_idx(235)) = 'OMPHIL  '
       if (trsp_idx(236).gt.0) GMNAME(trsp_idx(236)) = 'OMPHOB  '
       if (trsp_idx(237).gt.0) GMNAME(trsp_idx(237)) = 'OMPHIL  '

       if (trsp_idx(240).gt.0) GMNAME(trsp_idx(240)) = 'BCPHOB  '
       if (trsp_idx(241).gt.0) GMNAME(trsp_idx(241)) = 'BCPHIL  '
       if (trsp_idx(242).gt.0) GMNAME(trsp_idx(242)) = 'BCPHOB  '
       if (trsp_idx(243).gt.0) GMNAME(trsp_idx(243)) = 'BCPHIL  '
       if (trsp_idx(244).gt.0) GMNAME(trsp_idx(244)) = 'BCPHOB  '
       if (trsp_idx(245).gt.0) GMNAME(trsp_idx(245)) = 'BCPHIL  '
       
       if (trsp_idx( 73).gt.0) GMNAME(trsp_idx( 73)) = 'SULFUR  '
       if (trsp_idx( 62).gt.0) GMNAME(trsp_idx( 62)) = 'NH4FIN  '
       if (trsp_idx( 63).gt.0) GMNAME(trsp_idx( 63)) = 'NH4CRS  '
       if (trsp_idx( 64).gt.0) GMNAME(trsp_idx( 64)) = 'NO3FIN  '
       if (trsp_idx( 65).gt.0) GMNAME(trsp_idx( 65)) = 'NO3CRS  '

       write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
       write(6,'(a)') f90file//':'//subr// &
            ' Will write data to netCDF every 3rd hour'
    end if

    !// Strings for filenames
    write(cyear,'(i4.4)') JYEAR
    write(cday,'(i3.3)') NDAY

    !// String for time label
    time_label='hours since yyyy-mm-dd 00:00:00'
       write(time_label(13:23),'(i4.4,a1,i2.2,a1,i2.2)') &
            JYEAR,'-',JMON,'-',JDATE

    !// Output of air density
    !// All variables will be divided by volume in gm_dump_nc, so we
    !// send in air mass (AIR).
    if(LDUMP3HRS(2)) then
       ncfilename = 'air'//cyear//'_'//cday//'.nc'
       TMPFIELD(:,:,:,1) = AIR(:,:,:)
       call gm_dump_nc(1, 'air_density', TMPFIELD(:,:,:,1), 1, (/1/), &
            NMET, 1, (/1/), time_label, ncfilename)
    end if

    if (LOSLOCTROP .and. LDUMP3HRS(3)) then
       !// Output of selected tropospheric tracers + air (netCDF)
       ncfilename = 'trp'//cyear//'_'//cday//'.nc'
       call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
            NMET, trp_nr, TRP_LIST, time_label, ncfilename)
    end if

    if (LOSLOCTROP .and. LDUMP3HRS(13)) then
       !// Output of selected non transported tracers + air (netCDF)
       ncfilename = 'ntr'//cyear//'_'//cday//'.nc'
       call gm_dump_nc(NOTRPAR, XTNAME, XSTT, TRACER_ID_MAX, Xtrsp_idx, &
            NMET, ntr_nr, NTR_LIST, time_label, ncfilename)
    end if

    if (LSULPHUR .and. LDUMP3HRS(4)) then
       !// Output of sulphur
       ncfilename = 'sul'//cyear//'_'//cday//'.nc'
       call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
            NMET, sul_nr, SUL_LIST, time_label, ncfilename)
    end if

    if (LSALT .and. LDUMP3HRS(5)) then
       !// Output of salt
       ncfilename = 'slt'//cyear//'_'//cday//'.nc'
       call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
            NMET, slt_nr, SLT_LIST, time_label, ncfilename)
    end if

    if (LDUST .and. LDUMP3HRS(6)) then
       !// Output of dust
       ncfilename = 'min'//cyear//'_'//cday//'.nc'
       call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
            NMET, min_nr, MIN_LIST, time_label, ncfilename)
    end if

    if (LNITRATE .and. LDUMP3HRS(7)) then
       !// Output of nitrate
       ncfilename = 'snn'//cyear//'_'//cday//'.nc'
       call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
            NMET, nit_nr, SNN_LIST, time_label, ncfilename)
    end if

    if (LBCOC) then
       !// Output of bio BC/OC
       if (LDUMP3HRS(8)) then
          ncfilename = 'bio'//cyear//'_'//cday//'.nc'
          call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
               NMET, bio_nr, BIO_LIST, time_label, ncfilename)
       end if
       if (LDUMP3HRS(9)) then
          ncfilename = 'moa'//cyear//'_'//cday//'.nc'
          call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
               NMET, moa_nr, MOA_LIST, time_label, ncfilename)
       end if
       !// Output of fossil fuel BC/OC
       if (LDUMP3HRS(10)) then
          ncfilename = 'ffc'//cyear//'_'//cday//'.nc'
          call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
               NMET, ffc_nr, FFC_LIST, time_label, ncfilename)
       end if
       !// Output of fossil fuel BC/OC
       if (LDUMP3HRS(11)) then
          ncfilename = 'bfc'//cyear//'_'//cday//'.nc'
          call gm_dump_nc(NPAR, GMNAME, STT, TRACER_ID_MAX, trsp_idx, &
               NMET, bfc_nr, BFC_LIST, time_label, ncfilename)
       end if
    end if

    !// Output of SOA (only total)
    if (LSOA .and. LDUMP3HRS(12)) then
       !// SOA - ANT
       TMPFIELD(:,:,:,1) = STT(:,:,:,trsp_idx(soa_ant_list(1)))
       do N = 2, soa_ant_nr
          TMPFIELD(:,:,:,1) = TMPFIELD(:,:,:,1) &
               + STT(:,:,:,trsp_idx(soa_ant_list(N)))
       end do

       !// SOA - NAT
       TMPFIELD(:,:,:,2) = STT(:,:,:,trsp_idx(soa_nat_list(1)))
       do N = 2, soa_nat_nr
          TMPFIELD(:,:,:,2) = TMPFIELD(:,:,:,2) &
               + STT(:,:,:,trsp_idx(soa_nat_list(N)))
       end do

       ncfilename = 'soa'//cyear//'_'//cday//'.nc'
       call gm_dump_nc(2, (/'SOAANT    ','SOANAT    '/), TMPFIELD, 2, (/1,2/), &
            NMET, 2, (/1,2/), time_label, ncfilename)

    end if

    !// --------------------------------------------------------------------
  end subroutine dump3hrs
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine gm_dump_nc( &
       NPAR_IN,     & !I Number of components in STT_IN (not necessarily STT)
       TNAME_IN,    & !I Name of all the tracers sent in
       STT_IN,      & !I 4-D array with components
       IREPMX,      & !I Number of chemical components in MTC_IN
       MTC_IN,      & !I Transport order
       nbr_steps,   & !I Counter for netCDF timestep (i.e. NMET)
       NCOMP,       & !I Number of compenents to put out (size of ICOMP_LIST)
       ICOMP_LIST,  & !I Array with component names to dump
       time_label,  & !I Correct time label (hours since yyyy-mm-dd 00:00:00)
       filename)      !I Name of file (minYYYY_DDD.nc)
    !// --------------------------------------------------------------------
    !// Purpose: Write data to netcdf output file.
    !// 
    !//
    !// Author: Alf Grini, 2002
    !// To CTM3: Amund Sovde, 2012
    !// --------------------------------------------------------------------
    use netcdf
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only:  AREAXY, XDGRD, YDGRD, ZGRD, NRMETD
    use cmn_met, only:  ZOFLE
    use ncutils, only:  handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    !// IM,JM,LM are fixed to global
    !// MTC_IN may differ from MTC
    !// STT_IN may differ from STT
    !// 
    integer, intent(in)           :: NPAR_IN                         !Number of components
    character(len=10),intent(in)  :: TNAME_IN(NPAR_IN)               !Tracer names
    real(r8), intent(in)          :: STT_IN(IPAR,JPAR,LPAR,NPAR_IN)  !4-D fields
    integer, intent(in)           :: NCOMP                           !Number of components in list
    integer, intent(in)           :: ICOMP_LIST(NCOMP)               !Numbers of comps to put out
    integer, intent(in)           :: IREPMX                          !Number of species in MTC_IN
    integer,intent(in)            :: MTC_IN(IREPMX)                  !Comp.order (not necessarily MTC)
    character(len=*),intent(in)   :: filename                        !Name of file (minYYYY_DDD.nc)
    integer, intent(in)           :: nbr_steps                       !Step number
    character(len=80), intent(in) :: time_label                      !Time label or unit
    !// Output
    !// No output from this file
  
    !// Locals
    character(len=3)        :: cday                            !Name of start day in character
    character(len=4)        :: cyear                           !Name of start year in character
    character(len=10)       :: field_name                      !Name of field 
    integer                 :: cnt_lon_lat_lev_time(4)         !Count memory for added tracer field
    integer                 :: dim_lon_lat_lev_time(4)         !Dimension id for field
    integer                 :: lat_dim_id                      !Dimension id for latitude
    integer                 :: lat_id                          !Variable id for latitude
    integer                 :: lev_dim_id                      !Dimension id for level
    integer                 :: lev_id                          !Variable id for level
    integer                 :: lon_dim_id                      !Dimension id for longitude
    integer                 :: lon_id                          !Variable id for longitude
    integer                 :: nlats                           !Number of latitudes from nc-file
    integer                 :: nlevs                           !Number of levels from nc-file
    integer                 :: nlons                           !Number of longitudes from nc-file
    integer                 :: ncid                            !fileno for output netcdf file
    integer                 :: nsteps                          !Number of timesteps already in nc-file
    real(r8)                :: rtcdmp(IPAR,JPAR,LPAR)          !3D field to dump to file
    character(len=10)       :: TCNAME_AVG(NPAR_IN)             !Tracer name in the order not following MTC
    real(r8)                :: time                            !Time steps
    integer                 :: time_dim_id                     !Dimension id for time
    integer                 :: time_id                         !Variable id for time
    integer                 :: TCINFO(NPAR_IN,2)               !Tracer info (compno,dim_id,comp_id)
    integer                 :: srt_lon_lat_lev_time(4)         !Start value for tracer
    integer                 :: srt_time(1)                     !Start value for time
    integer                 :: status                          !Status for netcdf file 0 =OK, other value = error
    integer :: I, J, L, N
    real(r8) :: xdgrd_loc(ipar) !// To put out xdgrd inreasing monotonically

    integer, parameter :: nc4deflate=9 !1:141597536, 5:131449992 9:128082865
    integer, parameter :: nc4shuffle=1
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gm_dump_nc'
    !// --------------------------------------------------------------------
    
    !// Get the names and components in MTC order
    DO N = 1, NCOMP
       !write(6,*)'N,IC,MT,NP',N,ICOMP_LIST(N),MTC_IN(ICOMP_LIST(N)),NPAR
       TCNAME_AVG(N) = TNAME_IN(MTC_IN(ICOMP_LIST(N))) !Component name
       TCINFO(N,1)   = MTC_IN(ICOMP_LIST(N))        !Transport number
       if ( TCINFO(N,1) .le. 0) then 
          write(6,'(a,i5,a)') f90file//':'//subr//': COMPONENT ', &
               ICOMP_LIST(N), 'is not included.'
          stop
       end if
    END DO
  
    !// longitudes 0:360
    do i = 1, ipar
       if (XDGRD(i) .lt. 0._r8) then
          xdgrd_loc(i) = XDGRD(i) + 360._r8
       else
          xdgrd_loc(i) = XDGRD(i)
       end if
    end do

    !// No errors so far
    status=nf90_noerr 
       
    if (nbr_steps .eq. 1) then
       !First time this is done. Need to define variables

       !// Create file
       write(6,*)'creating file: ',filename,status
       !// Clobber means you can overwrite existing data (netcdf3 only):
       !// Use nf90_netcdf4 to create netcdf4:
       status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
       if (status/=nf90_noerr) call handle_error(status,'in creating file')

       !//File headers
       status=nf90_put_att(ncid,nf90_global,'title','3h output fields from Oslo CTM3')
       if (status/=nf90_noerr) call handle_error(status,'file header')
       status=nf90_put_att(ncid,nf90_global,'Modelinfo1','Oslo CTM3 is a 3D Chemical Transport Model')
       if (status/=nf90_noerr) call handle_error(status,'modelifo1')
       status=nf90_put_att(ncid,nf90_global,'Modelinfo2','Oslo CTM3 is driven by ECMWF meteorological data')
       if (status/=nf90_noerr) call handle_error(status,'modelinfo2')
       status=nf90_put_att(ncid,nf90_global,'contactinfo','For errors, contact CICERO')
       if (status/=nf90_noerr) call handle_error(status,'contactinfo')

       !// Define dimensions (JM, IM, LM, time)
       status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define lat dim')
       status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define lon dim')
       status = nf90_def_dim(ncid,"lev",LPAR,lev_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define lev dim')
       !// Defining dimension time of length unlimited
       status = nf90_def_dim(ncid,"time",nf90_unlimited,time_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define time dim')
     
       !// Defining the combined id for a field (lon / lev /lat /time)
       dim_lon_lat_lev_time(1)=lon_dim_id
       dim_lon_lat_lev_time(2)=lat_dim_id
       dim_lon_lat_lev_time(3)=lev_dim_id
       dim_lon_lat_lev_time(4)=time_dim_id
     
       !// Defining the lon/lat/lev-variable
       status = nf90_def_var(ncid,"lon",nf90_float,lon_dim_id,lon_id)
       if (status/=nf90_noerr) call handle_error(status,'define lon variable')
       status = nf90_def_var(ncid,"lat",nf90_float,lat_dim_id,lat_id)
       if (status/=nf90_noerr) call handle_error(status,'define lat variable')
       status = nf90_def_var(ncid,"lev",nf90_float,lev_dim_id,lev_id)
       if (status/=nf90_noerr) call handle_error(status,'define lev variable')
       status = nf90_def_var(ncid,"time",nf90_float,time_dim_id,time_id)
       if (status/=nf90_noerr) call handle_error(status,'define time variable')

       !// Putting attributes to /lon/lat/lev variables
       status = nf90_put_att(ncid,lon_id,'units','degree_east')
       if (status/=nf90_noerr) call handle_error(status,'attribute units lon')
       status = nf90_put_att(ncid,lat_id,'units','degree_north')
       if (status/=nf90_noerr) call handle_error(status,'attribute units lat')
       status = nf90_put_att(ncid,lev_id,'units','mb')
       if (status/=nf90_noerr) call handle_error(status,'attribute units lev')

       !// Putting attributes to time variable
       status = nf90_put_att(ncid,time_id,'units',time_label)
       if (status/=nf90_noerr) call handle_error(status,&
            f90file//':'//subr//':attribute units time')
       
       !// Define the tracer field variables
       DO N = 1, NCOMP
        
          field_name = TCNAME_AVG(N)
          status = nf90_def_var(ncid,field_name,nf90_float, &
               dim_lon_lat_lev_time, TCINFO(N,2))      !Defining TCINFO(N,2)
                                                       !As the variable_id for field
          if (status/=nf90_noerr) call handle_error(status,'define tracer variable')
          !// Deflate netcdf4
          status = nf90_def_var_deflate(ncid,TCINFO(N,2),nc4shuffle,1,nc4deflate)

          if (status/=nf90_noerr) call handle_error(status,'define tracer variable deflate')
        
          !// Add text descriptions and units
          status = nf90_put_att(ncid,TCINFO(N,2),'units','g/m3')        
          if (status/=nf90_noerr) call handle_error(status,'attribute units tracer')
        
          !// Put more attributes here if you want to . Ex: long_name
        
       END DO  !Loop on averageable variables
     
       !// End definition mode
       status = nf90_enddef(ncid)
       if (status/=nf90_noerr) call handle_error(status,'end defmode')
     
       !// Putting the lon/lat/lev variables
       !write(6,*)'Putting variables for lon/lev/lat. They do not change'
       status = nf90_put_var(ncid,lon_id,xdgrd_loc)
       if (status/=nf90_noerr) call handle_error(status,'putting lon')
       status = nf90_put_var(ncid,lat_id,ydgrd)
       if (status/=nf90_noerr) call handle_error(status,'putting lat')
       status = nf90_put_var(ncid,lev_id,ZGRD)
       if (status/=nf90_noerr) call handle_error(status,'putting lev')
       !write(6,*)'done writing lon/lev/lat dimensions'
     
    else                      !// THE FILE HAS BEEN USED BEFORE

       !// Open the existing file
       status = nf90_open(filename, nf90_write, ncid)
       if (status/=nf90_noerr) call handle_error(status,'open existing file')
     
       !// Inquire dimension ids
       status = nf90_inq_dimid(ncid,"lat",lat_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'getting lat')
       status = nf90_inq_dimid(ncid,"lon",lon_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'getting lon')
       status = nf90_inq_dimid(ncid,"lev",lev_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'getting lev')
       status = nf90_inq_dimid(ncid,"time",time_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'getting time')
     
       !// Inquire dimensions
       status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
       if (status/=nf90_noerr) call handle_error(status,'inq lat dim')
       if (nlats/=JPAR) then
          write(6,*)'Month avg file reports JM = ',nlats, JPAR
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
       if (status/=nf90_noerr) call handle_error(status,'inq lon dim')
       if (nlons/=IPAR) then
          write(6,*)'Month avg file reports IM = ',nlons,IPAR
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,lev_dim_id,len=nlevs)
       if (status/=nf90_noerr) call handle_error(status,'inq lev dim')
       if (nlevs/=LPAR) then
          write(6,*)'Month avg file reports LM = ',nlevs,LPAR
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
       if (status/=nf90_noerr) call handle_error(status,'inq time dim')
       if (nsteps.lt.0) then
          write(6,*)'Month avg file reports already added nsteps = ',nsteps
          stop
       end if
          
       do N = 1, NCOMP           
          !// Check variable id
          field_name = TCNAME_AVG(N)
          !// TCINFO(N,1) is MTC_IN(ICOMP_LIST(N))
          !// TCINFO(N,2) is variable_id for N in nc-file
          status = nf90_inq_varid(ncid,field_name,TCINFO(N,2))
          if (status/=nf90_noerr) call handle_error(status,'inq tracer varid')
       end do

       !// Get variable id for time
       status = nf90_inq_varid(ncid,"time",time_id) 
       if (status/=nf90_noerr) call handle_error(status,'inq time varid')
     
    end if  !// if (nbr_steps .eq. 1) then
  

    !// PUT THE VARIABLES TO FILE
    !// --------------------------------------------------------------------

    !// For the tracer fields:
    !// Defining how far to count for each time a data set is added
    !write(6,*)'setting count vecor'
    cnt_lon_lat_lev_time = (/IPAR , JPAR , LPAR , 1/)
    !// Defining where to start adding the new time step
    !write(6,*)'setting start vector'
    srt_lon_lat_lev_time = (/1, 1, 1, nbr_steps/)

    !// Start value for new time step
    srt_time(1) = nbr_steps
    !write(6,*)'srt_time is set',nbr_steps
    time = real(nbr_steps-1, r8)*24._r8 / real(NRMETD, r8) !3.0
    !write(6,*)'time = ',time

    status = nf90_put_var(ncid,time_id,time,start=srt_time)
    if (status/=nf90_noerr) call handle_error(status,'putting time')

    write(6,'(a)') '* Putting tracer data to nc-file'
    DO N = 1, NCOMP
    
       !//Go through the input array and pick the right IDs, fields and names
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR

                !!// Field for component ICOMPLIST(N)
                RTCDMP(I,J,L) = 1.e3_r8   &                         !kg-->g
                     *STT_IN(I,J,L,TCINFO(N,1))  &                  !kg
                     /(AREAXY(I,J) * (ZOFLE(L+1,I,J)-ZOFLE(L,I,J))) !m^3
             end do
          end do
       end do

       !Convert to real 4 to not get very small numbers
       !This is silly; let nf90 take care of it.
       !!RTCDMP(:,:,:) = real(RTCDMP(:,:,:),4)

       write(6,'(a,i3,1x,2(i4,1x),a8,a,es20.12)') '  Tracer/ID/trsp_idx ', &
            N,ICOMP_LIST(N),TCINFO(N,1),TCNAME_AVG(N), &
            ' total [kg]:',sum(STT_IN(:,:,:,TCINFO(N,1)))

       !// Put variable tracer field
       status = nf90_put_var(ncid,  &           !File id
            tcinfo(N,2), &                      !field_id for netCDF file (should match id set in def_var) 
            RTCDMP,                          &  !Tracer field (REAL4)
            start=srt_lon_lat_lev_time,      &  !starting point for writing
            count=cnt_lon_lat_lev_time )        !Counts how many bytes written
       if (status/=nf90_noerr) call handle_error(status,'putting tracer data')

    END DO
  
    !// close netcdf file
    status = nf90_close(ncid)
    if (status/=nf90_noerr) call handle_error(status,'close file')

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// --------------------------------------------------------------------
  end subroutine gm_dump_nc
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
end module gmdump3hrs
!//=========================================================================
