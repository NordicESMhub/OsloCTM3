!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Initialize the CTM.
!//=========================================================================
module initialize
  !//-----------------------------------------------------------------------
  !// MODULE: initialize
  !// DESCRIPTION: Initializes CTM.
  !//
  !// Contains:
  !//   subroutine INPUT
  !//   subroutine SETUP_SPECIES
  !//   subroutine SETUP_UNF_OUTPUT
  !//   subroutine report_zeroinit
  !//-----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'initialize.f90'
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine INPUT(NDAYI,NDAYE)
    !//---------------------------------------------------------------------
    !// Read std input(5) file that contains basic data for CTM run.
    !// Read file names for additional data sets and call the readers:
    !// tracer specifications and supporting data,
    !// met field data and file/directory info.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: IPAR, IPARW, JPAR, JPARW, LPAR, LPARW, IDBLK, &
         JDBLK, MPBLK, MPIPAR, MPJPAR, LOSLOCSTRAT, LE90, IMDIV, &
         NTDPAR, NTBPAR, NSBPAR, NSTPAR, NPAR, NRMETD
    use cmn_chem, only: INFILE_T, INFILE_WET, INFILE_DRY, INFILE_EMIS, &
         INFILE_POLAR_O3LOSS, INFILE_RES, INFILE_MEGAN, INFILE_LIGHTNING
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW, &
         IYEAR, JDAY, JYEAR, JMON, JDATE, GMTAU, TMON, TMET, &
         JDAY_NEXT, JDATE_NEXT, JMON_NEXT, JYEAR_NEXT, &
         LMTSOM, CFLLIM, NBLX, NSCX, NDPX, NROPSM, NRCHEM, &
         LLPYR, LYEAR, LFIXMET, LCONT, START_AVG, XDEDG, YDEDG, &
         ETAAW, ETABW, ETAA, ETAB, LMMAP, IMEPZ, LNCR, LDUMP3HRS, &
         SOLDEC, SOLDIS, IDTLN, modelTimeIntegrated, LSTOM1HRS, &
         LDLYSCAV
         !LDLYTOT, LDLY2D, Lbrd, Lls, Lc,n Ldry, Lsto 
    use cmn_diag, only: RUNTITLE, LFLXDG, JDO_C, JDO_T, JDO_A, JDO_X, &
         JDO_S, NTND, TLDIAG, &
         NTND_SOURCE, NTND_BNDLYR, NTND_DRYDEP, NTND_CHEM, &
         NTND_LSSCAV, NTND_CNSCAV, NTND_WADV, NTND_UVADV, &
         LBGA1, LBGA2, LBGTA, &
         NBOXD, TBOXD, IBOXD, JBOXD, KBOXD, NBOXLT, &
         NRGBLT, NRGLTD, &
         IBOXDMP, JBOXDMP, LBOXDMP, &
         NBOXS, TSTAX, STLAT, STLNG, LTSTNSV, &
         TLTAX, TLTRST, LTLAT, TLTRGL, LTLNG, LTSTN1, LTSTN2, ILTX, JLTX, &
         LTGBL1, LTGBL2, LTGBLSV, &
         LBGT1, LBGT2, LBGTS, ISTA, JSTA, GM0_LT, &
         metTypeInfo, resolutionInfo
    use cmn_met, only: MPATH1, MPATH2, MFILE3, metTYPE, metCYCLE, metREVNR, &
         MET_ROOT, MYEAR, PMEANW, HnativeRES, VRES, VRESW, PPFDPATH, PPFDFILE
    use cmn_sfc, only: LANDUSE_IDX, LANDUSE_YEAR, fileLandSurfTypeFrac, &
         LAI_YEAR, fileLAI, ZOI_YEAR, fileZOI, LDDEPmOSaic, DDEP_PAR, &
         LGSMAP, fileGSMAP, NLCAT
    use cmn_fjx, only: INFILE_FJX_SPEC, INFILE_FJX_SCAT, INFILE_FJX_AERO, &
         INFILE_FJX_JS, INFILE_FJX_O1D, INFILE_FJX_CLIM
    use cloudjx, only: CLDFLAG, NRANDO, RANSEED, &
         LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ, LCLDAVG
    use grid, only: SET_GRID, SET_GRID_VERT, SET_MEAN_PSFC, &
         DIAG_LTSTN, SURF_IN, DIAGBLK, DIAG_LTGL, &
         read_landSurfTypeFrac, read_LAI, read_ZOI, read_ZOI_LAI, &
         read_growing_season
    use utilities, only: calendar, get_soldecdis, get_free_fileid
    use cmn_oslo, only: LJCCYC
    use chem_oslo_rates, only: inSAD
    use input_oslo, only: clear_oslo_variables, read_tracer_list
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Output
    integer, intent(out) ::  NDAYI, NDAYE

    !// Locals
    character(len=80) :: INFILE1,INFILE2,INFILE3, TITLE
    real(r8) :: XLNG(2), YLAT(2), GM00Z
    integer :: I,J,L,M,N, IMX,JMX,LMX, NOPSTL, K, I1,I2,J1,J2,IFNR,IOS
    integer :: POLAVG(25)
    integer :: II, JMPOLR
    logical :: LISLSCP2, LCLDFLAG
    character(len=160) :: fileDDEPpar, fileGS
    character(len=70) :: TITCLD(8)
    character(len=4) :: CCNR,CRNR,NRES
    !// --------------------------------------------------------------------
    !// Temporare ascii array
    character(len=8), dimension(29,NLCAT) :: tmpDDEP_PAR
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'input'
    !//---------------------------------------------------------------------

    !// Initialise seconds counter
    modelTimeIntegrated = 0._r8


    !// Check capability for 'degridding' with 2x, 3x, ... box size in I or J
    !// params to reduce horizontal resolution of read-in met fields.
    !// integer, parameter :: NDEGM=2, IDGRD=2, JDGRD=2, NDGRD=IDGRD*JDGRD
    !// should set IM = IPARW/IDGRD, JM = (JPARW-2)/JDGRD+2
    !// if IM < IPARW then degrade resolution, must have mod (IPARW,IPAR)=0
    !// if JM < JPARW then degrade resolution, must have mod (JPARW,JPAR)=0

    if (mod(IPARW,IPAR) .ne. 0 .or. (mod(JPARW,JPAR) .ne. 0)) then
       write(6,'(a)') f90file//':'//subr// &
         ': Check degridding: IPARW IPAR JPARW JPAR'
       stop 'STOP in '//subr
    end if
    if (IDBLK*MPIPAR .lt. IPAR)  then
       write(6,'(a)') f90file//':'//subr//': Mismatch MPIPAR vs IPAR'
       write(6,'(a,i5)') '  Number of OMP blocks in longitude: ',MPIPAR
       write(6,'(a,i5)') '  Longitude dimension in OMP block:  ',IDBLK
       write(6,'(2(a,i5))') '  Total longitudnal grids:',IDBLK*MPIPAR, &
            ' is less than IPAR',IPAR
       stop 'STOP in '//subr
    end if

    if (JDBLK*MPJPAR .lt. JPAR)  then
       write(6,'(a)') f90file//':'//subr//': Mismatch MPJPAR vs JPAR'
       write(6,'(A,I5)') '  Number of OMP blocks in latitude: ',MPJPAR
       write(6,'(A,I5)') '  Latitude dimension in OMP block:  ',JDBLK
       write(6,'(2(A,I5))') '  Total latitudinal grids: ',JDBLK*MPJPAR, &
            ' is less than JPAR',JPAR
       stop 'STOP in '//subr
    end if

    !// setup OpenMp blocks
    do J = 1,MPJPAR
       do I = 1,MPIPAR
          N         = (J-1)*MPIPAR + I
          MPBLKIB(N) = (I-1)*IDBLK + 1
          MPBLKIE(N) = min(IPAR,I*IDBLK)
          MPBLKJB(N) = (J-1)*JDBLK + 1
          MPBLKJE(N) = min(JPAR,J*JDBLK)
       end do
    end do
    if (mod(MPBLKIE(1)-MPBLKIB(1)+1, IMDIV) .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Check OMP-block size (IMDIV) for DYN2W'
       write(6,'(A,I5)') 'Longitude length of block: ',MPIPAR
       write(6,'(A,I5)') 'Columns to stack (IMDIV) ',IMDIV
       stop 'STOP in '//subr
    end if

    !// Resolution strings
    if (IPARW .eq. 64) then
       HnativeRES = 'T21'
    else if (IPARW .eq. 128) then
       HnativeRES = 'T42'
    else if (IPARW .eq. 196) then
       HnativeRES = 'T63'
    else if (IPARW .eq. 320) then
       HnativeRES = 'T159'
    else if (IPARW .eq. 640) then
       HnativeRES = 'T319'
    else if (IPARW .eq. 360) then
       HnativeRES = '1x1'
    end if
    write(VRESW(1:3),'(a1,i2.2)') 'L',LPARW !// Metdata resolution
    write(VRES(1:3),'(a1,i2.2)') 'L',LPAR   !// Evt. degraded resolution


    !// Clear tracer and air mass
    AIR(:,:,:)   = 0._r8
    STT(:,:,:,:) = 1.e-20_r8
    SUT(:,:,:,:) = 0._rMom
    SVT(:,:,:,:) = 0._rMom
    SWT(:,:,:,:) = 0._rMom
    SUU(:,:,:,:) = 0._rMom
    SVV(:,:,:,:) = 0._rMom
    SWW(:,:,:,:) = 0._rMom
    SUV(:,:,:,:) = 0._rMom
    SUW(:,:,:,:) = 0._rMom
    SVW(:,:,:,:) = 0._rMom

    !// Clear additional arrays for chemistry, because they are
    !// used below.
    !// Note that additional initialisation is done in input_oslo,
    !// called from pmain after subroutine INPUT is done.
    call clear_oslo_variables()

    !// Read PRIMARY CTM DATA (unit=5)
    !//---------------------------------------------------------------------

    read(5,'(a80)') RUNTITLE
      write(6,'(a10,a)') 'CTM run>>>', trim(RUNTITLE)
    read(5,*)
    read(5,'(i5)') IYEAR
    read(5,'(2i5)') NDAYI
    read(5,'(2i5)') NDAYE
      write(6,'(a,3i8)') ' base yr, day begin/end:',IYEAR,NDAYI,NDAYE
    GMTAU  = 0._r8
    read(5,*)
    !// *** LCLDAVG and LCLDRAN combinations: AVG:T/F; QUAD:F/F; RAN:F/T ***
    read(5,'(l5)') LCLDQMD
    read(5,'(l5)') LCLDQMN
    read(5,'(l5)') LCLDRANA
    read(5,'(l5)') LCLDRANQ
    !// New cloud parameters cloudjx
    !CLDFLAG = 0
    !do I = 1,8
    !   read(5,'(L5,I3,1X,A70)') LCLDFLAG,II,TITCLD(I)
    !   if (LCLDFLAG)  then
    !      CLDFLAG = I
    !      NRANDO  = II
    !   end if
    !end do
    read(5,'(i5)') RANSEED
      write(6,'(a,i5)') ' RANSEED',RANSEED
    read(5,'(3i5)') NROPSM,NRCHEM
      write(6,'(a,3i5)') ' NRMETD/NROPSM/NRCHEM',NRMETD,NROPSM,NRCHEM


    !// CTM3: J-values only for first NRCHEM (once every NOPS)?
    read(5,'(l5)') LJCCYC
    !// SOM limiter and CFL criteria
    read(5,'(i5)') LMTSOM
      LMTSOM = min(3,max(1, LMTSOM))
    read(5,'(f5.3)') CFLLIM

    !// Meteorolgy info (should match paths below!)
    read(5,*)
    read(5,'(a)') metTYPE
    read(5,'(2i5)') metCYCLE, metREVNR
      write(6,'(A,A)')  ' Meteorology type:  ', METTYPE
      write(6,'(A,i3)') '   Cycle number:    ', metCYCLE
      write(6,'(A,i3)') '   Revision number: ', metREVNR
    !// Meteorology paths and names (must be given by '/path/')
    read(5,*) MET_ROOT

    !// The rest of the metdata path is set from metTYPE/metCYCLE/metREVNR
    write(CCNR,'(i4)') metCYCLE
    CCNR = adjustl(CCNR)
    write(CRNR,'(i4)') metREVNR
    CRNR = adjustl(CRNR)
    write(NRES,'(i4)') JPARW/2
    NRES = 'N'//trim(adjustl(NRES))

    !// Set up paths and file names.
    !// Important: Must match positions hardcoded in metdata_xxx.f90
    if (trim(metTYPE) .eq. 'ECMWF_IFS') then
       !// ECMWF IFS - old type
       write(6,'(A)') '   File type: EC binary'
       MPATH1 = trim(MET_ROOT)//'cy'//trim(CCNR)//'_' !'r'//trim(CRNR)//'_'
       !// Second part of directory path
       MPATH2 = 'YYYY/'//trim(HnativeRES)//'_data/'
       MFILE3 = 'ECDDD'//trim(HnativeRES)
       !// Meteorology type info to put on diagnostic files
       metTypeInfo = 'ECMWF IFS cycle '//trim(CCNR)

    else if (trim(metTYPE) .eq. 'ECMWF_oIFS') then
       !// ECMWF Open IFS standard EC binary format
       write(6,'(A)') '   File type: EC binary'
       MPATH1 = trim(MET_ROOT)//'cy'//trim(CCNR)//'r'//trim(CRNR)//'_'
       !// Second part of directory path
       MPATH2 = 'YYYY/'//trim(HnativeRES)//trim(VRESW)//trim(NRES)//'/'
       MFILE3 = 'ECopenYYYYDDD'//trim(HnativeRES)//trim(VRESW)// &
            'C'//trim(CCNR)//'a'
       !// Meteorology type info to put on diagnostic files
       metTypeInfo = 'ECMWF OpenIFS cycle '//trim(CCNR)//' revision '//trim(CRNR)

    else if (trim(metTYPE) .eq. 'ECMWF_oIFSnc4') then
       !// ECMWF Open IFS - netcdf4 files
       write(6,'(A)') '   File type: netcdf4'
       MPATH1 = trim(MET_ROOT)//'cy'//trim(CCNR)//'r'//trim(CRNR)//'nc4_'
       !// Second part of directory path
       MPATH2 = 'YYYY/'//trim(HnativeRES)//trim(NRES)//trim(VRESW)
       !// File name: ECopenIFSc38r1_yYYYYmMMdDDhHH_T159N80L60.nc
       MFILE3 = 'ECopenIFSc'//trim(CCNR)//'r'//trim(CRNR)// &
            '_yYYYYmMMdDDhHH_'//trim(HnativeRES)//trim(NRES)//trim(VRESW)
       !// Meteorology type info to put on diagnostic files
       metTypeInfo = 'ECMWF OpenIFS cycle '//trim(CCNR)//' revision '//trim(CRNR)

    else
       !// Unknown type
       write(6,'(a)') f90file//':'//subr//': Unknown metTYPE'
       write(6,'(a,a)')  '  metTYPE: ', metTYPE
       write(6,'(a,i5)') '  metCYCLE: ', metCYCLE
       write(6,'(a,i5)') '  metREVNRycle: ', metREVNR
       write(6,'(a,a)')  '  MET_ROOT:       ', trim(MET_ROOT)
       metTypeInfo = 'NOT DEFINED!!!'
       stop 'STOP in '//subr
    end if

      write(6,'(A,A)') '   MET_ROOT:       ',trim(MET_ROOT)
      write(6,'(A,A)') '   MPATH1:         ',trim(MPATH1)
      write(6,'(A,A)') '   MPATH2:         ',trim(MPATH2)
      write(6,'(A,A)') '   MFILE3:         ',trim(MFILE3)

    !// Set metdata type label
    if (MFILE3(1:2) .ne. 'EC') then
       write(6,'(a)') f90file//':'//subr// &
            ': Unknown metTYPE file: '//trim(MFILE3)
       write(6,'(a)') '  metTYPE is: '//trim(metTYPE)
       stop 'STOP in '//subr
    end if


    !// Resolution info
    if (IPAR .eq. IPARW .and. JPAR .eq. JPARW .and. LPAR .eq. LPARW) then
       resolutionInfo = 'Model resolution is the same as meteorology data resolution'
    else
       resolutionInfo = 'Model resolution differs from meteorology data resolution, see variabes IPARW, JPARW and LPARW.'
    end if

    read(5,'(l5)') LLPYR
    read(5,'(l5)') LFIXMET
      write(6,'(A,L2)') ' allow for leap year (LLPYR): ',LLPYR
      write(6,'(A,L2)') ' recycle met fields (LFIXMET):',LFIXMET



    !// Check cloud cover parameters
    !if (CLDFLAG.lt.1 .or. CLDFLAG.gt.8) then
    !   write(6,'(a,i5)') f90file//':'//subr// &
    !      ': Error in select cloud cover type: CLDFLAG: ',CLDFLAG
    !   stop 'STOP in '//subr
    !else
    !   write(6,'(A,I3)') 'Cloud '//TITCLD(CLDFLAG),NRANDO
    !end if
    LCLDAVG = .false.
    if (LCLDQMD)  then
       if (LCLDQMN .or. LCLDRANA .or. LCLDRANQ) then
          write(6,'(a)') f90file//':'//subr// &
            ': Only one of LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ can be true'
          stop 'STOP in '//subr
       else
          write(6,'(A,L2)') ' Mid point of quadrature cloud cover'//&
               ' ICAs (LCLDQMD):',LCLDQMD
       end if
    else if (LCLDQMN) then
       if (LCLDRANQ .or. LCLDRANA) then
          write(6,'(a)') f90file//':'//subr// &
            ': Only one of LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ can be true'
          stop 'STOP in '//subr
       else
          write(6,'(A,L2)') ' Mean quadrature cloud cover ICAs'//&
               ' (LCLDQMN):',LCLDQMN
       end if
    else if (LCLDRANA) then
       if (LCLDRANQ) then
          write(6,'(a)') f90file//':'//subr// &
            ': Only one of LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ can be true'
          stop 'STOP in '//subr
       else
          write(6,'(A,L2,A,I6)') ' Random selected from all cloud'//&
               ' cover ICAs (LCLDRANA):',LCLDRANA,', RANSEED:',RANSEED
       end if
    else if (LCLDRANQ) then
       write(6,'(A,L2,A,I6)') ' Random selected from 4 Mean '// &
            'quadrature cloud cover ICAs (LCLDRANQ):',LCLDRANQ, &
            ', RANSEED:',RANSEED
    else
       LCLDAVG = .true.
       write(6,'(A,L2)') ' use averaged cloud cover (LCLDAVG):',LCLDAVG
    end if


      write(6,'(a)') '--- OpenMP Blocks: no. IBeg/IEnd JBeg/JEnd'
      write(6,'(5i8)') (M,MPBLKIB(M),MPBLKIE(M),MPBLKJB(M),MPBLKJE(M),M=1,MPBLK)

    !// update the calendar to the new, upcoming day = 1st day of CTM run
    !call CALENDR(IYEAR,NDAYI, LLPYR,LFIXMET,MYEAR, &
    !     JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)
    call calendar(IYEAR,NDAYI, LLPYR,LFIXMET,MYEAR, &
         JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET, &
         JYEAR_NEXT,JDAY_NEXT,JMON_NEXT,JDATE_NEXT)
    call get_soldecdis(JDAY,SOLDEC,SOLDIS)

    !// No need to read resolution
    read(5,'(i5)') JMPOLR
      write(6,'(A,i5)') ' JMPOLR: ', JMPOLR
    read(5,'(f5.2)') GM00Z
      write(6,'(A,f5.1)')  ' GM0000: ', GM00Z

    !// Read in eta levels of met-field grid (collapse later)
    !// if LM < LPARW (met field layers) use LMMAP to collapse
    !// met-layers on CTM L
    read(5,*)
    read(5,*)
    do L = 1, LPARW+1
       read(5,'(5x,i5,f20.10,f17.10)') LMMAP(L),ETAAW(L),ETABW(L)
    end do

    !// efine polar filter to average U,V,T,Q,P over EW boxes (check
    !// for multiple)
    read(5,'(5x,25i3)') POLAVG
    IMEPZ(:) = 1
    do J = 1, 25
       if (J .lt. JPAR/2) then
          if (POLAVG(J) .gt. 1) then
             if (mod(IPAR,POLAVG(J)) .eq. 0) then
                IMEPZ(J) = POLAVG(J)
                IMEPZ(JPAR+1-J) = POLAVG(J)
             end if
          end if
       end if
    end do
    write(6,'(a9,30I3)') ' IMEPZ: ',(IMEPZ(J), J=1,30)

    !// Define type of boundary-layer mixing, dry deposition,
    !// washout/rainout code
    read(5,'(i5)') NBLX
      write(6,'(a,i5)') ' BL type:     ',NBLX
    read(5,'(i5)') NDPX           
      write(6,'(a,i5)') ' dry dep type:',NDPX
    read(5,'(i5)') NSCX
      write(6,'(a,i5)') ' scavng  type:',NSCX

    !// SET UP HORIZONTAL GRID (I,J)
    call SET_GRID(GM00Z, JMPOLR)

    !// Annual mean pressure field in CTM 
    read(5,*)
    read(5,*)  INFILE1  ! annual mean pressure field
    call set_mean_psfc(trim(infile1))

    !// SET UP VERTICAL GRID (I,J,L)
    call SET_GRID_VERT()

    call inSAD('Indata_CTM3/SAD_12m_1x1.nc')

    !// Land surface type fractions
    read(5,*)
    read(5,'(I5)')  LANDUSE_IDX
    if (LANDUSE_IDX .eq. 2) then
       write(6,'(a)') f90file//':'//subr//': Land surface type fractions from MODIS'
    else if (LANDUSE_IDX .eq. 3) then
       write(6,'(a)') f90file//':'//subr//': Land surface type fractions from CLM4-input'
    else
       write(6,'(a,i5)') f90file//':'//subr//': No such LANDUSE_IDX',LANDUSE_IDX
       stop
    end if
    read(5,'(I5)')  LANDUSE_YEAR ! year to get (9999 for metdata)
    read(5,*)  fileLandSurfTypeFrac  ! land surface type fraction

    !// Roughness lenght
    read(5,*)
    read(5,'(I5)')  ZOI_YEAR ! year to get (9999 for metdata, 0000 for climatology)
    read(5,*)  fileZOI       ! surface roughness length data set

    !// Leaf area index
    read(5,*)
    read(5,'(I5)')  LAI_YEAR ! year to get (9999 for metdata, 0000 for climatology)
    read(5,*)  fileLAI       ! leaf area index

    !// Get surface data (vegetation map, etc.)
    !call SURF_IN(INFILE1,INFILE2,INFILE3,LANDUSE_IDX.eq.2)

    !// Must also be read each year if LANDUSE_YEAR = 9999
    call read_landSurfTypeFrac()
    call read_LAI()
    call read_ZOI()

    !call read_ZOI_LAI(INFILE2,INFILE3)

    !// Dry deposition scheme based on Simpson et al. (2012) and parameter list therein
    read (5,*)
    read(5,'(l5)') LDDEPmOSaic ! Switch between new and old dry deposition scheme
    if (LDDEPmOSaic) then
       write(6,'(a)') f90file//':'//subr// ': DRYDEP: Using mOSaic scheme.'
    else
       write(6,'(a)') f90file//':'//subr// ': DRYDEP: Using old scheme.'
    end if
    !// Read the location of the dry deposition parameter list taken from Simpson et al. (2012)
    read(5,*)  fileDDEPpar
    !// Read the location of the PPFD files
    read(5,*)
    read(5,*)  PPFDPATH
    read(5,*)  PPFDFILE
    !// Initialize dry deposition parameters
    if (LDDEPmOSaic) then
       IFNR = get_free_fileid()
       open(IFNR,file=fileDDEPpar,Status='OLD',action='read',IOSTAT=IOS)
       if (IOS .eq. 0) then
          write(6,'(a)') '** Reading dry deposition parameters from '//trim(fileDDEPpar)
       else
          write(6,'(a)') f90file//':'//subr//': File not found: '//trim(fileDDEPpar)
          stop
       end if
       ! Read the table header
       read(IFNR, *) 
       ! Read the whole table as ascii (none floating point value in column one)
       read(IFNR, *) tmpDDEP_PAR
       ! Split the temporary table and save the data
       read(tmpDDEP_PAR(2:,:),'(f10.0)') DDEP_PAR
       write(6,*) tmpDDEP_PAR
       close(unit=ifnr)
    end if
    !// Growing season
    read(5,*)
    read(5,'(l5)') LGSMAP ! Switch between preprocessed GDAY/GLEN and fixed latitude based callculation
    read(5,*) fileGSMAP   ! Growing season related maps
    if (LGSMAP) then
       write(6,'(a)') f90file//':'//subr// ': GROWING_SEASON: Using preprocessed data.'
       write(6,'(a)') f90file//':'//subr// ': GROWING_SEASON WARNING: Only called during initialization. For runs exceeding one year a restart and change of file in *.inp is currently needed!'
       call read_growing_season()
    else
       write(6,'(a)') f90file//':'//subr// ': GROWING_SEASON: Using fixed scheme.'
    end if
    

    !// CTM3: Tracer specific input - previously in Ltracer.inp
    read(5, '(a)') TITLE
      write(6,'(a)') TITLE
    read(5,'(l2,1x,I2)')  LCONT, START_AVG
    write(6,'(a,l1)') 'Continuation run: ', LCONT
    if (LCONT) then
       write(6,'(a,l1)') '  Initializing tracers with restart file'
    else
       if (START_AVG .eq. 0) then
          write(6,'(a)') '  Initializing tracers with zero'
       else if (START_AVG .eq. 1) then
          write(6,'(a)') '  Initializing tracerswith CTM3 month average'
       else
          write(6,'(a,i5)') '  Unknown START_AVG; ',START_AVG
       end if
    end if

    ! Switch between old (.sav) and new (.nc) restart files
    read(5,'(l5)') LNCR 

    read (5,*)
    read(5,*)  INFILE_T  ! tracer list
    write(6,'(a)') 'tracer list: '//trim(INFILE_T)

    !// Will read tracer list now
    call read_tracer_list(infile_T)
    write(6,'(a)') '-- Resume reading input file'

    read (5,*)
    read(5,*)  INFILE_WET  ! wet scav list
    write(6,'(a)') 'wet scavenging list: '//trim(INFILE_WET)
    read (5,*)
    read(5,*)  INFILE_DRY  ! dry scav list
    write(6,'(a)') 'UCI dry scavenging list: '//trim(INFILE_DRY)
    read (5,*)
    read(5,*)  INFILE_EMIS  ! tracer emission list
    write(6,'(a)') 'tracer emission list: '//trim(INFILE_EMIS)

    !// Read FASTJ related tables
    read (5,*)
    read(5,*)  INFILE_FJX_SPEC  ! spectral data
    write(6,'(a)') 'FASTJ spectral data: '//trim(INFILE_FJX_SPEC)
    read(5,*)  INFILE_FJX_SCAT  ! mie scattering data
    write(6,'(a)') 'FASTJ mie scattering data: '//trim(INFILE_FJX_SCAT)
    read(5,*)  INFILE_FJX_AERO  ! aerosol optical data
    write(6,'(a)') 'FASTJ aerosol optical data: '//trim(INFILE_FJX_AERO)
    read(5,*)  INFILE_FJX_JS  ! FASTJ->Oslo CTM3 species mapping
    write(6,'(a)') 'FASTJ->Oslo CTM3 species mapping: '//trim(INFILE_FJX_JS)
    read(5,*)  INFILE_FJX_O1D  ! FASTJ->Oslo CTM3 species mapping O1D
    write(6,'(a)') 'FASTJ->Oslo CTM3 species mapping O1D: '//trim(INFILE_FJX_O1D)
    read(5,*)  INFILE_FJX_CLIM  ! T and O3 climatology
    write(6,'(a)') 'FASTJ temperature and ozone climatology: '//trim(INFILE_FJX_CLIM)
    
    !// Read Polar ozone loss data
    read (5,*)
    read(5,*)  INFILE_POLAR_O3LOSS  ! polar ozone loss data
    write(6,'(a)') 'Polar ozone loss data: '//trim(INFILE_POLAR_O3LOSS)

    !// Read chemical rates in specified resolution
    read (5,*)
    read(5,*)  INFILE_RES  ! chemical rates in given resolution
    write(6,'(a)') 'Chemical rates: '//trim(INFILE_RES)

    !// Read MEGAN tables
    read (5,*)
    read(5,*)  INFILE_MEGAN  ! MEGAN input tables
    write(6,'(a)') 'MEGAN tables: '//trim(INFILE_MEGAN)

    !// Read Lightning distributions
    read (5,*)
    read(5,*)  INFILE_LIGHTNING  ! lightning distributions
    write(6,'(a)') 'Lightning distributions: '//trim(INFILE_LIGHTNING)

    read (5,*)
    !// DIAGNOSTICS
    !//---------------------------------------------------------------------
    read (5,*)
    read (5,*)
    read (5,*)
    !// Switch on 3 hourly output  
    read(5,*) LDUMP3HRS
    if (LDUMP3HRS(1)) then
       write(6,*) f90file//':'//subr//': Diagnostic output: Dump3hrs (all/trp/sul/slt/min/nit/bio/moa/ffc/bfc/ntr)', LDUMP3HRS
    end if
    read (5,*)
    ! Switch on 1 hourly output of GstO3 and FstO3
    read(5,'(l5)') LSTOM1HRS 
    if (LSTOM1HRS) then
       write(6,'(a)') f90file//':'//subr//': Diagnostic output: Stom1hrs'
    end if
    read (5,*)
    read (5,*)
    !// Switch to scavenging daily output (tot/2d-brd/ls/cn/dry/sto)
    read(5,*) LDLYSCAV
    if (LDLYSCAV(1) .or. LDLYSCAV(2)) then
       write(6,*) f90file//':'//subr//': Diagnostic output: DlyScav (tot/2d-brd/ls/cn/dry/sto) ', LDLYSCAV
    end if

    !// CTM3 override; calculate when stratospheric chemistry is on
    !// and when E90-tracer is used (need LPAUZ)
    LFLXDG = LOSLOCSTRAT .and. LE90

    !// restart / continuation file calendar
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_C

    write(6,'(a)') ' <<<<<restart/cont save calendar:'
    write(6,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_C

    !// tendency budgets: labels, calendar for write & reset
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_T

    write(6,'(a)') ' <<<<<tendency budget calendar:'
    write(6,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_T

    read(5,'(i5)') NTND
    if (NTND .gt. NTDPAR) then
       write(6,'(a)') f90file//':'//subr// &
            ': # tendency items > NTDPAR: ',NTND,NTND,NTDPAR
       stop 'STOP in '//subr
    end if
    !// So far there are 8 processes to diagnose. In addition we will
    !// print out the sum of all tendencies '>SUMS<', and instant
    !// tracer field 'INST-M'.
      write(6,'(a,i2,a)') 'Will do tendencies for ',NTND,' processes'

    TLDIAG(:) = ' '
    NTND_SOURCE = 0
    NTND_BNDLYR = 0
    NTND_DRYDEP = 0
    NTND_LSSCAV = 0
    NTND_CNSCAV = 0
    NTND_WADV   = 0
    NTND_UVADV  = 0
    do I = 1, NTND
       read(5,'(i5,3x,a6)') II, TLDIAG(II)
       if (trim(TLDIAG(II)) .eq. 'SOURCE') then
          NTND_SOURCE = II
       else if (trim(TLDIAG(II)) .eq. 'BNDLYR') then
          NTND_BNDLYR = II
       else if (trim(TLDIAG(II)) .eq. 'DRYDEP') then
          NTND_DRYDEP = II
       else if (trim(TLDIAG(II)) .eq. 'CHEMIS') then
          NTND_CHEM = II
       else if (trim(TLDIAG(II)) .eq. 'LSSCAV') then
          NTND_LSSCAV = II
       else if (trim(TLDIAG(II)) .eq. 'CNSCAV') then
          NTND_CNSCAV = II
       else if (trim(TLDIAG(II)) .eq. 'W_ADV_') then
          NTND_WADV = II
       else if (trim(TLDIAG(II)) .eq. 'UV_ADV') then
          NTND_UVADV = II
       else
          write(6,'(a,3i5)') f90file//':'//subr// &
            'Unknown process diag. name '//TLDIAG(II)
          stop 'STOP in '//subr
       end if
       write(6,'(a,i3)') ' Assiged process diag. '//TLDIAG(II)// &
            ' to process # ',II
    end do

    !// Need this section to be species, not transport numbers
    read(5,*)
    read(5,'(18x,L1,18(1x,10L1))') LBGA1,(LBGT1(I),I=1,NPAR)
    read(5,'(18x,L1,18(1x,10L1))') LBGA2,(LBGT2(I),I=1,NPAR)

      write(6,'(a9,11a7)')' tends: ', (TLDIAG(I), I=1,NTND)
      write(6,'(a9,2x,L1,18(1x,10L1))') ' 1D budgt',LBGA1,LBGT1
      write(6,'(a9,2x,L1,18(1x,10L1))') ' 2D budgt',LBGA2,LBGT2


    !// 3D tendency blocks: box(1) always = global
    read(5,'(5x,i3)') NBOXD
      write(6,'(a,i5)')' <<<<<budget tendencies boxes:', NBOXD
    if (NBOXD .gt. NTBPAR) then
       write(6,'(a,3i5)') f90file//':'//subr// &
            ': # boxes > NTBPAR: ',NBOXD,NBOXD,NTBPAR
       stop 'STOP in '//subr
    end if

    TBOXD(1) = ' GLOBAL '
    IBOXD(1,1)  = 1
    JBOXD(1,1)  = 1
    KBOXD(1,1)  = 1
    IBOXD(2,1)  = IPAR
    JBOXD(2,1)  = JPAR
    KBOXD(2,1)  = LPAR
    NBOXD = NBOXD+1
    do N = 2,NBOXD
       read(5,'(10x,a15,4f8.2,2i5)') TBOXD(N), &
            (XLNG(M),M=1,2),(YLAT(M),M=1,2), KBOXD(1,N),KBOXD(2,N)
       KBOXD(1,N)  = max(1,KBOXD(1,N))
       KBOXD(2,N)  = min(LPAR,KBOXD(2,N))
       call DIAGBLK(XLNG,YLAT,XDEDG,YDEDG, &
            IBOXD(1,N),JBOXD(1,N),IPAR,JPAR,IDTLN)
    end do

    do M = 1,MPBLK
       do K = 1,NBOXD
          J1 = max(JBOXD(1,K),MPBLKJB(M))
          J2 = min(JBOXD(2,K),MPBLKJE(M))
          I1 = max(IBOXD(1,K),MPBLKIB(M))
          I2 = min(IBOXD(2,K),MPBLKIE(M))
          II = mod(IBOXD(2,K)-1,IPAR) + 1
          if (J1 .gt. J2)  then
             LBOXDMP(K,M) = .false.
          else
             JBOXDMP(1,K,M) = J1 - MPBLKJB(M) + 1
             JBOXDMP(2,K,M) = J2 - MPBLKJB(M) + 1
             if (II .ge. IBOXD(1,K)) then
                if (I1 .gt. I2)  then
                   LBOXDMP(K,M) = .false.
                else
                   LBOXDMP(K,M) = .true.
                   IBOXDMP(1,K,M) = I1 - MPBLKIB(M) + 1
                   IBOXDMP(2,K,M) = I2 - MPBLKIB(M) + 1
                end if
             else
                if ( (MPBLKIE(M)-MPBLKIB(M)+1) .eq. IPAR) then
                   LBOXDMP(K,M) = .true.
                   IBOXDMP(1,K,M) = IBOXD(1,K)
                   IBOXDMP(2,K,M) = IBOXD(2,K)
                else
                   if (MPBLKIE(M) .ge. IBOXD(1,K))  then
                      I2 = min(IPAR,MPBLKIE(M))
                   end if
                   if (II .ge. MPBLKIB(M)) then
                      I1 = max(1,MPBLKIB(M))
                      I2 = min(II,MPBLKIE(M))
                   end if
                   if (I1 .gt. I2)  then
                      LBOXDMP(K,M) = .false.
                   else
                      LBOXDMP(K,M) = .true.
                      IBOXDMP(1,K,M) = I1 - MPBLKIB(M) + 1
                      IBOXDMP(2,K,M) = I2 - MPBLKIB(M) + 1
                   end if
                end if
             end if
          end if
 
       end do
    end do

    do N = 1,NBOXD
       write(6,'(a5,i4,4x,a15,3(a7,2i5))')' Tbox',N,TBOXD(N), &
            '  I1/I2',IBOXD(1,N),IBOXD(2,N),'  J1/J2',JBOXD(1,N),JBOXD(2,N), &
            '  K1/K2',KBOXD(1,N),KBOXD(2,N)
       do M = 1,MPBLK
          if (LBOXDMP(N,M))  then
             write(6,'(A8,I4,7H  IB/IE,2I4,7H  JB/JE,2I4, 2X,11H  IB_B/IE_B,2I4,11H  JB_B/JE_B,2I4)') &
                  'IJ-BLK',M,MPBLKIB(M),MPBLKIE(M),MPBLKJB(M),MPBLKJE(M), &
                  (IBOXDMP(I,N,M),I=1,2),(JBOXDMP(J,N,M),J=1,2)
          end if
       end do
    end do

    !// 3-D average tracer (m/m):  calendar for write & reset
    !// ---------------------------------------------------------------------
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_A
      write(6,'(a)') ' <<<<<averages calendar:'
      write(6,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_A


    !// STEFLUX: STE output (and instantaneous mixing ratios, not
    !// implemented in CTM3)
    !// ---------------------------------------------------------------------
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_X
    read(5,*)
      write(6,'(a)') ' <<<<<STE output calendar:'
      write(6,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_X

    !// 1-D average tracer profile (m/m):  calendar above, stdout, uses
    !// tendency boxes
    !// ---------------------------------------------------------------------
    read(5,*)
    read(5,'(19x,18(1x,10L1))') (LBGTA(I),I=1,NPAR)
      write(6,'(a,18(1x,10L1))') ' 1D tracer m/m:',LBGTA

    !// time series at specified stations, full L-profiles
    read(5,*)
    read(5,*)
    read(5,*)
    read(5,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_S

    read(5,'(i5,14x,18(1x,10l1))') NBOXS,(LBGTS(I),I=1,NPAR)
    if (NBOXS .gt. NSBPAR) then
       write(6,'(a,3i5)') f90file//':'//subr// &
            ': # boxes > NSBPAR: ',NBOXS,NBOXS,NSBPAR
       stop 'STOP in '//subr
    end if
    if (NRMETD*NROPSM .gt. NSTPAR) then
       write(6,'(a,3i5)') f90file//':'//subr// &
            ': # time steps > NSTPAR: ',NRMETD,NROPSM,NSTPAR
       stop 'STOP in '//subr
    end if

    do N = 1, NBOXS
       read(5,'(10x,a10,2x,2f8.2)') TSTAX(N),STLAT(N),STLNG(N)
       !// check for sign change at dateline (eg, box[1] from +175. to -175.
       if (STLNG(N) .gt. XDEDG(IDTLN)) STLNG(N) = -180._r8
       II = IDTLN + 1
       I  = mod(II-1,IPAR) + 1
       do while (STLNG(N) .gt. XDEDG(I))
          II = II+1
          I  = mod(II-1,IPAR) + 1
       end do
       II = II-1
       ISTA(N) = mod(II-1,IPAR) + 1
       STLAT(N) = min(YDEDG(JPAR+1), max(YDEDG(1), STLAT(N)))
       J = 2
       do while  (STLAT(N) .gt. YDEDG(J))
          J = J+1
       end do
       JSTA(N) = J-1
    end do

    write(6,*) ' <<<<<time series profile stations:',NBOXS

    write(6,'(a)') ' <<<<<station-time series calendar:'
    write(6,'(5x,31i1/5x,29i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1/5x,31i1/5x,30i1/5x,31i1/5x,30i1/5x,31i1)') JDO_S

    do N = 1, NBOXS
       write(6,'(i5,2x,a10,f10.2,i4,f10.2,i4)') &
            N,TSTAX(N),STLAT(N),JSTA(N),STLNG(N),ISTA(N)
    end do

    !// Local hours
    NOPSTL = NROPSM * NRMETD
    do I = 1,IPAR+1
       GM0_LT(I) = dmod((XDEDG(I)/15._r8)+24._r8,24._r8)
    end do
    LTSTNSV(:,:) = .false.

    !// Stations to put out local time diags
    read(5,*)
    read(5,'(i5)') NBOXLT
    if (NBOXLT .gt. NSBPAR) then
       write(6,'(a,3i5)') f90file//':'//subr// &
            ': # boxes > NSBPAR: ',NBOXLT,NBOXLT,NSBPAR
       stop 'STOP in '//subr
    end if

    do N = 1, NBOXLT
       read(5,'(3x,2(a10,2x),2f8.2,2x,2f6.2)') TLTRST(N),TLTAX(N), &
            LTLAT(N),LTLNG(N),LTSTN1(N),LTSTN2(N)
       call DIAG_LTSTN(LTLNG(N),LTLAT(N),XDEDG,YDEDG,GM0_LT,LTSTN1(N), &
            LTSTN2(N),N,ILTX(N),JLTX(N),IPAR,JPAR,IDTLN,NOPSTL,LTSTNSV(1,N))
    end do

    write(6,*) ' <<<<<local time diag. stations:',NBOXLT
    do N = 1, NBOXLT
       write(6,'(i2,2(1x,a10),2(f10.2,i4),2x,2f6.2)') &
            N,TLTRST(N),TLTAX(N),LTLAT(N),JLTX(N),LTLNG(N),ILTX(N),&
            LTSTN1(N),LTSTN2(N)
    end do

    !// Tracers to put out at the selected stations
    read(5,*)
    read(5,'(i5)') NRGBLT
    write(6,'(i5)') NRGBLT
    do N = 1, NRGBLT
       read(5,'(3x,a10,2x,2f6.2)') TLTRGL(N), LTGBL1(N), LTGBL2(N)
       write(6,'(2A,2X,2F6.2)') &
            ' Read in 3-D local time diag. tracer and time interval ', &
            TLTRGL(N), LTGBL1(N), LTGBL2(N)
       call DIAG_LTGL(GM0_LT,LTGBL1(N),LTGBL2(N),IPAR,NOPSTL, &
            NRGLTD(N),LTGBLSV(1,1,N))
    end do

    close(5)

    !//---------------------------------------------------------------------
  end subroutine INPUT
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine SETUP_SPECIES(NDAY,NDAYI)
    !//---------------------------------------------------------------------
    !// Sets up the tracer arrays.
    !//
    !// Ole Amund Sovde, October 2015, January 2009
    !//---------------------------------------------------------------------
    use cmn_size, only: IPAR,JPAR,LPAR,NPAR, IPARW,JPARW,LPARW, NOTRPAR, &
         LSULPHUR, LBCOC, Le90
    use cmn_ctm, only: START_AVG, LCONT, IYEAR, LLPYR, LYEAR, LFIXMET, &
         JYEAR, JDAY, JMON, TMON, JDATE, TMET, SOLDEC, SOLDIS, &
         JDAY_NEXT, JDATE_NEXT, JMON_NEXT, JYEAR_NEXT, STT, AIR, &
         LNCR
    use cmn_diag, only: RUNTITLE
    use cmn_met, only: MYEAR
    use utilities, only: calendar, get_soldecdis
    use cmn_oslo, only: trsp_idx, ZEROINIT, XZEROINIT
    use bcoc_oslo, only: bcsnow_init
    use stt_save_load, only: load_restart_file, &
         oslo_con_run, oslo_con_run42, oslo_con_runxx, &
         restart_from_CTM3avg, oslo_restartfile_info
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NDAYI

    !// Locals
    integer           :: VERSION, IDW,JDW,LDW,ID,JD,LD,ND,XND,readInType
    character(len=80) :: FN_CON, TITLE
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'SETUP_SPECIES'
    !// --------------------------------------------------------------------

    write(6,'(a)') '-------------------------------------------------'// &
         '----------------------'


    if (LCONT) then

       !// Reload a continuation/restart run
       if (LNCR) then
          FN_CON = 'restart.nc'
          call load_restart_file(FN_CON, 0)

          !// You can add read-in of several files:
          !// call load_restart_file('additional_file.nc', 1)
       else
          !// Old read-ins
          FN_CON = 'restart.sav'
          !// Start from the same resolution as current run:
          !call OSLO_CON_RUN(FN_CON, 0)
          !// Start from different hor.res. (sav-file version > 2):
          !call OSLO_CON_RUNxx(FN_CON, 0)
          !// Start from T42 horizontal resolution (sav-file version 1):
          call OSLO_CON_RUN42(FN_CON, 0)
       end if

       !// Restart fields for BCsnow
       if (LBCOC) call bcsnow_init(NDAYI)

       !// The dates are not overwritten by restart file dates. This is
       !// to allow restarting a simulation from a differeny year.
       !// I keep the call to calendar anyway.
       !call CALENDR (IYEAR,NDAYI, LLPYR,LFIXMET,MYEAR, &
       !     JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)
       call calendar(IYEAR,NDAYI, LLPYR,LFIXMET,MYEAR, &
            JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET, &
            JYEAR_NEXT,JDAY_NEXT,JMON_NEXT,JDATE_NEXT)
       call get_soldecdis(JDAY,SOLDEC,SOLDIS)

    else

       !// Initialize tracer fields when LCONT is FALSE, i.e. when
       !// a restart field is not used. Use a monthly average.
       if (START_AVG .eq. 1) then
          !// CTM3 average
          write(6,'(a,i3,a)'), 'START_AVG .eq. ',START_AVG,' not defined'
          stop
          !call restart_from_CTM3avg(zeroinitLOCAL)
          !// If you need to restart from a T42-field, here is how to do it:
          !call restart_from_CTM3avg_T42(zeroinitLOCAL)
       else if (START_AVG .ne. 0) then
          !// Not defined
          write(6,'(a,i3,a)'), 'START_AVG .eq. ',START_AVG,' not defined'
          stop
       end if

       !// Restart fields for BCsnow
       if (LBCOC) call bcsnow_init(NDAYI)

    end if
    !// --------------------------------------------------------------------
  end subroutine SETUP_SPECIES
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine SETUP_UNF_OUTPUT(NDAY,NDAYI)
    !//---------------------------------------------------------------------
    !// Open the unformatted output files.
    !//
    !// Ole Amund Sovde, October 2015, January 2009
    !//---------------------------------------------------------------------
    use cmn_size, only: IPAR,JPAR,LPAR,NPAR, IPARW, NOTRPAR, &
         LSULPHUR, LBCOC, Le90
    use cmn_ctm, only: START_AVG, LCONT, IYEAR, LLPYR, LYEAR, LFIXMET, &
         JYEAR, JDAY, JMON, TMON, JDATE, TMET, SOLDEC, SOLDIS
    use cmn_diag, only: RUNTITLE
    use cmn_met, only: MYEAR
    use cmn_oslo, only: trsp_idx, ZEROINIT, XZEROINIT
    use bcoc_oslo, only: bcsnow_init
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NDAYI

    !// Locals
    logical           :: ex
    !//---------------------------------------------------------------------

    write(6,'(a)') '-------------------------------------------------'// &
         '----------------------'

    if (.not.LCONT)  then


       !// unit=23 = S.unf = unformatted time series profiles at stations
       !// done every op-split time (eg, 8, 24, 48 times per day)
       !// accumulates over the CTM run, daily writes set in calendar JDO_T
       open (23,file='S.unf',form='unformatted')
       write(23)  RUNTITLE
       call TSER_0(23)


       !// unit=24 = T.unf = unformatted 1D & 2D budget tendencies
       !// accumulates over the CTM run, between dates set in calendar JDO_T
       open (24,file='T.unf',form='unformatted')
       write(24)  RUNTITLE

       !// unit=25 = D.unf =instant dump of 3D (IJL) tracers & air (set by flags)
       !// done every op-split time (eg, 8, 24, 48 times per day)
       !// accumulates over the CTM run, per tracer flag, until K3DMAX reached
       open (25,file='D.unf',form='unformatted')
       write(25)  RUNTITLE
       write(25)  IPAR,JPAR,LPAR,NPAR

    else


       !// SERIES = 23 = unformatted times series
       inquire (File='S.unf',exist=ex)
       if (ex) then
          open (23,file='S.unf',status='old',form='unformatted',position='append')
       else
          open (23,file='S.unf',status='new',form='unformatted')
          write(23)  RUNTITLE
       end if

       !// TEND = unit=24 = unformatted budget tendencies
       inquire (File='T.unf',exist=ex)
       if (ex) then
          open (24,file='T.unf',status='old',form='unformatted',position='append')
       else
          open (24,file='T.unf',status='new',form='unformatted')
          write(24)  RUNTITLE
       end if
       !// 3-D = unit=25 = unformatted real*4 3-D dump of specified tracers & air
       !// done every op-split time, AIR given in kg (dry), STT in m/m
       inquire (File='D.unf',exist=ex)
       if (ex) then
          open (25,file='D.unf',status='old',form='unformatted',position='append')
       else
          open (25,file='D.unf',status='new',form='unformatted')
          write(25)  IPAR,JPAR,LPAR,NPAR
       end if

       LCONT = .false.
    end if

    !// --------------------------------------------------------------------
  end subroutine SETUP_UNF_OUTPUT
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine report_zeroinit(zeroinit)
    !// --------------------------------------------------------------------
    !// Report which tracers are set to zero.
    !//
    !// Ole Amund Sovde, November 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR, NOTRPAR
    use cmn_ctm, only: STT
    use cmn_chem, only: TNAME
    use cmn_oslo, only: trsp_idx, chem_idx, Xtrsp_idx, Xchem_idx, XTNAME, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: zeroinit(NPAR+NOTRPAR)
    !// Local variables
    integer :: N, TRACER_ID, K, TRSPNR, N1, N2
    !// --------------------------------------------------------------------

    K = 0   !// Total number of tracers with zeroinit
    N1 = 0  !// Number of zeroinitialized transported tracers
    N2 = 0  !// Number of zeroinitialized non-transported tracers

    !// Find N1 and N2
    do N = 1, NPAR + NOTRPAR
       if (N.le.NPAR) then
          if (zeroinit(N) .gt. 0) N1 = N1+1
       end if
       if (N.gt.NPAR) then
          if (zeroinit(N) .gt. 0) N2 = N2+1
       end if
    end do

    !// Report zero initializations - transported species
    if (N1 .gt. 0) then
       write(6,'(a)') 'ZERO-initialized transported components'
       write(6,'(2x,A2,1x,A6,4x,2x,A2,1x,A8,1x,A9,1x,A15)') &
            'NR','Tracer','ID','trsp_idx','Xtrsp_idx','TOTAL MASS (kg)'
       do N = 1, NPAR
          if (zeroinit(N) .eq. 1) then
             K = K + 1             !// Counting
             TRACER_ID = chem_idx(N)  !// Tracer ID
             write(6,'(1x,i3,1x,a10,1x,i3,1x,i8,11x,es12.6)') &
                  K,TNAME(N),TRACER_ID,N, sum(STT(:,:,:,N))
          end if
       end do
    end if

    !// Report zero initializations - non-transported species
    if (N2 .gt. 0) then
       write(6,'(a)') 'ZERO-initialized non-transported components'
       write(6,'(2x,A2,1x,A6,4x,2x,A2,1x,A8,1x,A9,1x,A15)') &
             'NR','Tracer','ID','trsp_idx','Xtrsp_idx','TOTAL MASS (kg)'
       do N = NPAR+1, NPAR + NOTRPAR
          if (zeroinit(N) .eq. 1) then
             K = K + 1             !// Counting
             TRSPNR = N-NPAR       !// Transport number
             if (TRSPNR .gt. 0) then
                TRACER_ID = Xchem_idx(TRSPNR) !// Tracer ID
                write(6,'(1x,i3,1x,a10,1x,i3,12x,i8,es12.6)') &
                     K,XTNAME(TRSPNR),TRACER_ID,TRSPNR,sum(XSTT(:,TRSPNR,:,:))
             end if
          end if
       end do
    end if

    !// Are all tracers initialized?
    if ((N1 + N2) .eq. 0) then
       write(6,'(a)') 'All tracers are initialized from restart file'
       write(6,'(a71)') '--------------------------------------------'// &
            '---------------------------'
       return
    else
       write(6,'(a)') 'Not all tracers are initialized from restart file'
       write(6,'(a,i3)') 'Transported species not initialized:     ',N1
       write(6,'(a,i3)') 'Non-transported species not initialized: ',N2

       !// Print out info about what "zero" is
       write(6,'(a)') 'NOTE: Zero means 1.d-20 per grid box!'
       write(6,'(a71)') '--------------------------------------------'// &
            '---------------------------'
    end if

    !// --------------------------------------------------------------------
  end subroutine report_zeroinit
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
end module initialize
!//-------------------------------------------------------------------------
