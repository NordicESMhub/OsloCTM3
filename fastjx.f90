!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, October 2015
!//=========================================================================
!// FastJX routines.
!//=========================================================================
module fastjx
  !//-----------------------------------------------------------------------
  !// MODULE: fastjx
  !// DESCRIPTION: fast-JX parameters, variables and routines.
  !//
  !//              Routines are from p-setc_oc.f, and thus fast-JX 6.7c.
  !//              WILL BE UPDATED for version 7.3 eventually.
  !//
  !// Contains:
  !//   subroutine PHOT_IN
  !//   subroutine RD_XXX
  !//   subroutine RD_MIE
  !//   subroutine RD_UM
  !//   subroutine RD_JS
  !//   subroutine RD_O1D
  !//   subroutine RD_PROF
  !//   subroutine SET_ATM
  !//
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8, r4
  !//-----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'fastjx.f90'
  !//-----------------------------------------------------------------------
  private
  public PHOT_IN, SET_ATM
  save
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine PHOT_IN
    !//---------------------------------------------------------------------
    !// Routine to initialise photolysis rate data.
    !// Will be updated for fastJX7
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_fjx, only: AER1N, AER2N, AER1P, AER2P, LTROP, &
         INFILE_FJX_SPEC, INFILE_FJX_SCAT, INFILE_FJX_AERO, &
         INFILE_FJX_JS, INFILE_FJX_O1D, INFILE_FJX_CLIM
    use utilities, only: get_free_fileid
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer  IPH, I, J, K
    !//---------------------------------------------------------------------

    !// Find a free file id
    IPH = get_free_fileid()

    !// Read in fast-J X-sections (spectral data) <<<<<< new fast-JX
    call RD_XXX(IPH,INFILE_FJX_SPEC)

    !// Read in aerosol/cloud scattering data <<<<<<<<<< new fast-JX
    call RD_MIE(IPH,INFILE_FJX_SCAT)

    !// Read in UMich aerosol optical data   <<<<<<<<<<< new fast-JX 6.1
    call RD_UM (IPH,INFILE_FJX_AERO)

    !// Read in labels of photolysis rates required
    !// >>>>> keyed to users chem code
    !// this is a tranfer map from the J_s automatically calculated in
    !// fast-JX onto the names and order in the users chemistry code
    !// CTM3: skip TSPECI
    call RD_JS(IPH,INFILE_FJX_JS) !,TSPECI)
    call RD_O1D(IPH,INFILE_FJX_O1D)

    !// Read in T & O3 climatology            >>>> general backup clim.
    call RD_PROF(IPH,INFILE_FJX_CLIM)

    !// Aerosol paths
    AER1N(:,:,:) = 0
    AER2N(:,:,:) = 0
    AER1P(:,:,:) = 0._r8
    AER2P(:,:,:) = 0._r8

    LTROP  = .false.       ! if .T. eliminates wavelengths

    !//---------------------------------------------------------------------
  end subroutine PHOT_IN
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine RD_XXX(NJ1,NAMFIL)
    !//---------------------------------------------------------------------
    !// Version FJX6.7 (checked December 2012)
    !// To f90 October 2015.
    !//---------------------------------------------------------------------
    !// Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections. 
    !//
    !// NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
    !//     WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengths
    !//     W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).
    !//  if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
    !//     W_ = 8, reverts to quick fix: fast-J (12-18) plus bin (5) scaled
    !//
    !//---------------------------------------------------------------------
    !// NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
    !// NJ1      Channel number for reading data file
    !//
    !// NJVAL    Number of species to calculate J-values for
    !// NWWW     Number of wavelength bins, from 1:NWWW
    !// WBIN     Boundaries of wavelength bins
    !// WL       Centres of wavelength bins - 'effective wavelength'
    !// FL       Solar flux incident on top of atmosphere (cm-2.s-1)
    !// QRAYL    Rayleigh parameters (effective cross-section) (cm2)
    !// QO2      O2 cross-sections
    !// QO3      O3 cross-sections
    !// Q1D      O3 => O(1D) quantum yield
    !// TQQ      Temperature for supplied cross sections
    !// QQQ      Supplied cross sections in each wavelength bin (cm2)
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_fjx, only: AER1N, AER2N, AER1P, AER2P, LTROP, &
         TQQ, NW1, NW2, X_, W_, WX_, NJVAL, &
         WL, FL, QRAYL, TITLEJ, QO2, QO3, Q1D, QQQ
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) :: NJ1
    character(len=*), intent(in) ::  NAMFIL

    integer :: I, J, JJ, K, IW, NQQQ, NWWW, NQRD
    character(len=7) :: TITLEJ2, TITLEJ3
    character(len=78) :: TITLE0
    character(len=20) :: FMT
    !//---------------------------------------------------------------------

    TQQ(:,:) = 0._r8

    !// ----------spectral data----set for new format data------------------
    !//     note that NQRD = # Xsects read in (must be  .le. X_) and
    !//               NJVAL = # fast-JX J-values derived from this
    !//     note NQQQ is not used outside this subroutine!

    !// SOLAR: Amund Sovde Haslerud
    !// If you want to interpolate FL between solar min and max, you need
    !// to read in min and max. FLmin and FLmax have to be defined in
    !// cmn_fjx.f90. See also routines at the bottom here, namely
    !// solar_flget and solar_flscale_read.
    !// I only used this for testing (ala study of Joanna Haigh, 
    !// doi:10.1038/nature09426) but it never went anywhere.
    !// The whole effect of solar min/max cannot be studied,
    !// since the direct effect on temperature is already part of
    !// meteorological data.
    !open (NJ1,FILE='FJX_solminmax.dat',status='old',form='formatted')
    !read (NJ1,100) TITLE0
    !read (NJ1,102) (FLmin(IW),IW=1,18)
    !read (NJ1,100) TITLE0
    !read (NJ1,102) (FLmax(IW),IW=1,18)
    !close(NJ1)
    !// Cannot adjust to 100% for all bins; leave at 80%.


    !// W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

    open (NJ1,FILE=NAMFIL,status='old',form='formatted')
    read (NJ1,'(a)') TITLE0
    read (NJ1,101) NQRD, NWWW
101 format(10x,5i5)
    NW1 = 1
    NW2 = NWWW
    if (NQRD.gt.X_) then
       write(6,201) NQRD,X_
201 format(' Number of x-sections supplied to Fast-JX: ',i3,/, &
            ' Maximum number allowed (X_) only set to: ',i3, &
            ' - increase in cmn_jv.f')
       stop
    end if
    write(6,'(1X,A)') TITLE0
    write(6,'(2i8)') NQRD,NWWW

    !// ----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
    read (NJ1,102) (WL(IW),IW=1,NWWW)
    read (NJ1,102) (FL(IW),IW=1,NWWW)
    read (NJ1,102) (QRAYL(IW),IW=1,NWWW)
102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))

    !// ---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields
    !//    (each at 3 temps)
    read (NJ1,103) TITLEJ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
    read (NJ1,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
    read (NJ1,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)
103 format(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))

    read (NJ1,103) TITLEJ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
    read (NJ1,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
    read (NJ1,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

    read (NJ1,103) TITLEJ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
    read (NJ1,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
    read (NJ1,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

    do J = 1,3
       write(6,200) J,TITLEJ(J),(TQQ(I,J),I=1,3)
    end do
200 format(1x,' x-sect:',i3,a10,3(3x,f6.2))

    !// ---Read remaining species:  X-sections at 2 T_s
    JJ = 4
    do J = 4,NQRD
       read (NJ1,103) TITLEJ(JJ),TQQ(1,JJ),(QQQ(IW,1,JJ),IW=1,NWWW)
       read (NJ1,103) TITLEJ2,  TQQ(2,JJ),(QQQ(IW,2,JJ),IW=1,NWWW)

       if (W_.eq.18 .or. TITLEJ2(7:7).ne.'x') then
          !// include stratospheric J's (this also includes Cl
          !// and Br compounds!)
          write(6,200) JJ,TITLEJ(JJ),(TQQ(I,JJ),I=1,2)
          JJ = JJ+1
       end if

    end do

    NJVAL = JJ - 1
    NQQQ = NJVAL

    !// ---truncate number of wavelengths to do troposphere-only
    if (W_ .ne. WX_) then
       !// ---TROP-ONLY
       if (W_ .eq. 12) then
          write(6,'(a)') &
               ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
          NW2 = 12
          do IW = 1, 4
             WL(IW) = WL(IW+4)
             FL(IW) = FL(IW+4)
             QRAYL(IW) = QRAYL(IW+4)
             do K = 1, 3
                QO2(IW,K) = QO2(IW+4,K)
                QO3(IW,K) = QO3(IW+4,K)
                Q1D(IW,K) = Q1D(IW+4,K)
             end do
             do J = 4, NQQQ
                QQQ(IW,1,J) = QQQ(IW+4,1,J)
                QQQ(IW,2,J) = QQQ(IW+4,2,J)
             end do
          end do
          do IW = 5, 12
             WL(IW) = WL(IW+6)
             FL(IW) = FL(IW+6)
             QRAYL(IW) = QRAYL(IW+6)
             do K = 1, 3
                QO2(IW,K) = QO2(IW+6,K)
                QO3(IW,K) = QO3(IW+6,K)
                Q1D(IW,K) = Q1D(IW+6,K)
             end do
             do J = 4, NQQQ
                QQQ(IW,1,J) = QQQ(IW+6,1,J)
                QQQ(IW,2,J) = QQQ(IW+6,2,J)
             end do
          end do
          !// ---TROP-QUICK  (must scale solar flux for W=5)
       else if (W_ .eq. 8) then
          write(6,'(a)') &
               ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
          NW2 = 8
          do IW = 1, 1
             WL(IW) = WL(IW+4)
             FL(IW) = FL(IW+4)  * 2._r8
             QRAYL(IW) = QRAYL(IW+4)
             do K = 1, 3
                QO2(IW,K) = QO2(IW+4,K)
                QO3(IW,K) = QO3(IW+4,K)
                Q1D(IW,K) = Q1D(IW+4,K)
             end do
             do J = 4, NQQQ
                QQQ(IW,1,J) = QQQ(IW+4,1,J)
                QQQ(IW,2,J) = QQQ(IW+4,2,J)
             end do
          end do
          do IW = 2, 8
             WL(IW) = WL(IW+10)
             FL(IW) = FL(IW+10)
             QRAYL(IW) = QRAYL(IW+10)
             do K = 1, 3
                QO2(IW,K) = QO2(IW+10,K)
                QO3(IW,K) = QO3(IW+10,K)
                Q1D(IW,K) = Q1D(IW+10,K)
             end do
             do J = 4, NQQQ
                QQQ(IW,1,J) = QQQ(IW+10,1,J)
                QQQ(IW,2,J) = QQQ(IW+10,2,J)
             end do
          end do

       else
          write(6,*) ' number of used wavelengths wrong:',W_
          stop
       end if
    end if

    close(NJ1)


    !//---------------------------------------------------------------------
  end subroutine RD_XXX
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine RD_MIE(NJ1,NAMFIL)
    !//---------------------------------------------------------------------
    !// Version FJX6.7 (checked December 2012)
    !//---------------------------------------------------------------------
    !// aerosols/cloud scattering data set for fast-JX (ver 5.7r)
    !// >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
    !//---------------------------------------------------------------------
    !// NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
    !// NJ1      Channel number for reading data file
    !// NAA      Number of categories for scattering phase functions
    !// QAA      Aerosol scattering phase functions
    !// NK       Number of wavelengths at which functions supplied (set as 4)
    !// WAA      Wavelengths for the NK supplied phase functions
    !// PAA      Phase function: first 8 terms of expansion
    !// RAA      Effective radius associated with aerosol type
    !// SAA      Single scattering albedo
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_fjx, only: A_, JTAUMX, ATAU, ATAU0, &
         NAA, RAA, DAA, WAA, SAA, PAA, QAA
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) :: NJ1
    character(len=*), intent(in) ::  NAMFIL

    integer :: I, J, K
    character(len=78) :: TITLE0
    character(len=20) :: TITLAA(A_)
    !//---------------------------------------------------------------------

    open (NJ1,FILE=NAMFIL,status='old',form='formatted')

    read (NJ1,'(i2,a78)') NAA,TITLE0
    if (NAA .gt. A_) then
       write(6,*) ' too many scat-data sets:', NAA, A_
       stop
    end if
    read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0
    write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX
    read (NJ1,*)
      
    do J = 1, NAA
       read (NJ1,'(3x,a20,32x,f5.3,15x,f5.3)') TITLAA(J),RAA(J),DAA(J)
       do K = 1, 5  ! ver 6.0 extend to 5 ref wavelengths for mie-scat data
          read (NJ1,'(f4.0,f7.4,f7.4,7f6.3,1x,f7.3,f8.4)') &
               WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1._r8
       end do
    end do

    close(NJ1)

    write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):' &
                  ,(WAA(K,1),K=1,5)
    write(6,*) TITLE0
    do J = 1, NAA
       write(6,'(i3,1x,a8,7f8.3)') &
            J,TITLAA(J),RAA(J),DAA(J),(QAA(K,J),K=1,5)
    end do

    !//---------------------------------------------------------------------
  end subroutine RD_MIE
  !//-----------------------------------------------------------------------




  !//-----------------------------------------------------------------------
  subroutine RD_UM(NJ1,NAMFIL)
    !//---------------------------------------------------------------------
    !// Version FJX6.7 (checked December 2012)
    !//---------------------------------------------------------------------
    !// UMich aerosol optical data for fast-JX (ver 6.1+)
    !//---------------------------------------------------------------------
    !// NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
    !// NJ1      Channel number for reading data file
    !//---------------------------------------------------------------------
    use cmn_fjx, only: n_umset, UMAER
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) :: NJ1
    character(len=*), intent(in) ::  NAMFIL

    integer :: I, J, K, L
    character(len=78) :: TITLE0
    character(len=20) :: TITLUM(n_umset)
    !//---------------------------------------------------------------------

    open (NJ1,FILE=NAMFIL,status='old',form='formatted')

    read (NJ1,'(a78)') TITLE0
    write(6,*) 'UMichigan Aerosol optical data'
    write(6,*) TITLE0

    !// 33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
    !//       FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
    do L = 1, n_umset
       read(NJ1,'(a8)') TITLUM(L)
       !// 21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
       do K = 1, 21
          !// 6 wavelengths:
          !//   J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 4=600nm, 5=1000nm
          !// 3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
          read(NJ1,'(9f9.5,27x,6f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,5)
       end do
    end do

    close(NJ1)

    write(6,'(5(i5,1x,a8))') (L,TITLUM(L), L=1,n_umset)

    !//---------------------------------------------------------------------
  end subroutine RD_UM
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine RD_JS(NJ1,NAMFIL) !,TSPECI)
    !//---------------------------------------------------------------------
    !// Version FJX6.7 (checked December 2012)
    !//---------------------------------------------------------------------
    !// Reread the ratj.d file to map photolysis rate to reaction
    !// Read in quantum yield 'jfacta' and fastj2 label 'jlabel'
    !//---------------------------------------------------------------------
    !// JFACTA    Quantum yield (or multiplication factor) for photolysis
    !// JLABEL    label of J-value used in the main chem model
    !// JMAP      label of J-value used to match with fast-JX J's
    !// NRATJ     number of Photolysis reaction counter - must be .le. 'JVN_'
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_fjx, only: JVN_, JFACTA, JLABEL, NRATJ, JIND, NJVAL, TITLEJ
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) ::  NJ1
    character(len=*), intent(in) :: NAMFIL
    !//character*10, intent(in) :: TSPECI(JPSPEC)

    integer :: IPR, JR, JP, J, K, ioerr
    character(len=120) :: CLINE
    character(len=7) :: T_JX, T_CHEM
    !// for Oslo Chemistry instead of asad method.
    integer, parameter :: jpspj = 5
    character(len=10) acspj(jpspj)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'RD_JS'
    !//---------------------------------------------------------------------

    !// Reread the ratj.d file to map photolysis rate to reaction
    !// Read in quantum yield jfacta and fastj2 label jlabel
    open (NJ1,file=NAMFIL,status='old',form='formatted',iostat=ioerr)
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr//': Error opening ratefile: '//NAMFIL
       stop 'STOP in '//subr
    end if

    !// Read off header, assume no longer than 1000
    do J = 1, 1000
       read (NJ1,'(A)',iostat=ioerr)  CLINE
       if (ioerr .ne. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Error opening ratefile: '//NAMFIL//' header'
          stop 'STOP in '//subr
       end if
       if (CLINE(1:1).eq.'#') then
          !// Go to next line
          cycle
       else
          !// Finished headers
          exit
       end if
    end do
    !// Get one step back to get first reaction
    backspace(NJ1)

    !// Read through JVN_ entries.
    IPR = 0
    do J = 1, JVN_
       read (NJ1,'(A)',iostat=ioerr)  CLINE
       if (ioerr .ne. 0) then
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Error opening ratefile: '//NAMFIL//', line ',J
          stop 'STOP in '//subr
       end if

       if (CLINE(2:5).eq.'9999') then
          !// should in principle not happen: exit loop
          exit
       else if (CLINE(1:1).eq.'#') then
          !// Go to next line
          cycle
       else if (CLINE(5:5).eq.'$') then
          !// Go to next line
          cycle
       end if

       !// Store values
       IPR = IPR + 1
       read (CLINE,9902) (acspj(jp), JP=1,JPSPJ), JFACTA(IPR), JLABEL(IPR)
 9902 format(6x,4(1a10,1x),a10,18x,f5.1,2x,a7)
       write(6,'(4(1a10,1x),a10,f5.1,2x,a7)') &
            (acspj(JP), JP=1,JPSPJ),JFACTA(IPR), JLABEL(IPR)

       JFACTA(IPR) = JFACTA(IPR) / 100._r8
    end do

    !// Check if 9999 is next line
    read (NJ1,'(A)',iostat=ioerr)  CLINE
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Error reading 9999 '//NAMFIL
       stop 'STOP in '//subr
    end if
    if (CLINE(2:5).ne.'9999') then
       write(6,'(a)') f90file//':'//subr// &
            ': 9999 is missing '//NAMFIL
       write(6,'(a,2i5)') '  Number of entries read IPR/JVN_: ',IPR,JVN_
       stop 'STOP in '//subr
    end if

    close(NJ1)

    NRATJ = IPR

    !//---------------------------------------------------------------------
    !// compare Xsections titles with J-values listed in chem code (jratd.dat)
    !// map the J-values needed for chemistry (ratj.d) onto the fast-JX rates
    !// >>>>>>>>>>>>>>>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<<<<<
    !//      >>>this must now follow the read in of Xsects, etc<<<
    !//---------------------------------------------------------------------

    !// ---Zero / Set index arrays that map Jvalue(j) onto rates
    do K = 1, NRATJ
       JIND(K) = 0
       T_CHEM = JLABEL(K)
       do J = 1, NJVAL
          T_JX = TITLEJ(J)
          if (T_CHEM(1:6) .eq. T_JX(1:6)) then
             JIND(K) = J
          end if
       end do
    end do

    write(6,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
    do K = 1, NRATJ
       J = JIND(K)
       if (J.eq.0) then
          write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K), &
               ' has no mapping onto fast-JX'
       else
          write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K), &
               ' mapped onto fast-JX:',J,TITLEJ(J)
       end if
    end do

    !//---------------------------------------------------------------------
  end subroutine RD_JS
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine RD_O1D(NJ1,NAMFIL)
    !//---------------------------------------------------------------------
    !// Read the ratO1D.d file to map O3 + PHOTON --> O2 + O3(1D)
    !// photolysis rate to reaction
    !// Read in quantum yield 'jfactao' and fastj2 label 'jlabelo'
    !//
    !// jfactao    Quantum yield (or multiplication factor) for photolysis
    !// jlabelo    Reference label identifying appropriate J-value to use
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_fjx, only: JFACTAO, JLABELO, TITLEJ, JINDO, NJVAL
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) ::  NJ1
    character(len=*), intent(in) :: NAMFIL

    integer :: J, ioerr
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'RD_O1D'
    !//---------------------------------------------------------------------

    !// Read the ratO1D.d file to map photolysis rate to reaction
    !// Read in quantum yield jfactao and fastj2 label jlabelo
    open (NJ1,file=NAMFIL,status='old',form='formatted',iostat=ioerr)
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr//': Error opening ratefile: '//NAMFIL
       stop 'STOP in '//subr
    end if

    read (NJ1,*)
    read (NJ1,*)
    read (NJ1,'(77X,f5.1,2x,a7)')   JFACTAO, JLABELO
    write(6,'(f5.1,2x,a7)')     JFACTAO, JLABELO

    close(NJ1)

    JFACTAO = JFACTAO / 100._r8

    do J = 1, NJVAL
       if (JLABELO .eq. TITLEJ(J)) JINDO = J
    end do
    !//---------------------------------------------------------------------
  end subroutine RD_O1D
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine RD_PROF(NJ2,NAMFIL)
    !//---------------------------------------------------------------------
    !// Version FJX6.7 (checked December 2012)
    !//---------------------------------------------------------------------
    !// Routine to input T and O3 reference profiles
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_fjx, only: TREF, OREF
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) ::  NJ2
    character(len=*), intent(in) ::  NAMFIL

    integer :: IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
    real(r8) :: OFAC, OFAK
    character(len=78) :: TITLE0
    !//---------------------------------------------------------------------

    open (NJ2,file=NAMFIL,status='old',form='formatted')
    read (NJ2,'(A)') TITLE0
    read (NJ2,'(2I5)') NTLATS,NTMONS
    write(6,*) 'RD_PROF: Reading '//NAMFIL
    write(6,1000) NTLATS,NTMONS
1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')
    N216  = min(216, NTLATS*NTMONS)
    do IA = 1, N216
       read (NJ2,'(1X,I3,3X,I2)') LAT, MON
       M = min(12, max(1, MON))
       L = min(18, max(1, (LAT+95)/10))
       read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
       read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
    end do
    close (NJ2)

    !// Extend climatology to 100 km
    OFAC = exp(-2.e5_r8 / 5.e5_r8)
    do I = 32, 51
       OFAK = OFAC**(I-31)
       do M = 1, NTMONS
          do L = 1, NTLATS
             OREF(I,L,M) = OREF(31,L,M) * OFAK
          end do
       end do
    end do
    do L = 1, NTLATS
       do M = 1, NTMONS
          do I = 42, 51
             TREF(I,L,M) = TREF(41,L,M)
          end do
       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine RD_PROF
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine SET_ATM
    !//---------------------------------------------------------------------
    !// Sets up atmosphere for fast-JX, i.e. O3 from STT and also
    !// O3 above the model domain (from climatology).
    !//
    !// Oslo CTM3 with only tropospheric chemistry uses O3 climatology
    !// in stratosphere.
    !//
    !// Ole Amund Sovde, 28 April 2010
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, IPAR, JPAR, NPAR, LOSLOCSTRAT
    use cmn_ctm, only: ETAA, ETAB, YDGRD, AREAXY, STT, JMON
    use cmn_chem, only: TNAME
    use cmn_fjx, only: TREF, OREF, TOPT, TOPM, TOP3, &
         LPART, LPARM, LPAR3, DMS, DO3, MASFAC
    use cmn_met, only: P
    use cmn_parameters, only: AVOGNR
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Local arrays
    integer :: I, J, K, L, M, N, LMAX
    real(r8) :: PSTD(52), OREF2(51), TREF2(51), PJ(LPAR+1)
    real(r8) :: F0,T0,PB,PC,XC,DLOGP,PTOPCTM,PTOPATM,TOPMASS
    real(r8) :: PLMM1
    !//---------------------------------------------------------------------


    !//---code block to calculate mass, O3, T in atmosphere above top of CTM
    !//---  assumes that PTOP is same for all (I,J) = ETAA(LPAR+1)
    !//---  pressure levels for O3/T climatology (at 2 km in z*)
    PSTD(1) = 1000._r8
    PSTD(2) = 1000._r8 * 10._r8**(-1._r8 / 16._r8)
    DLOGP   = 10._r8**(-2._r8 / 16._r8)
    do K = 3, 51
       PSTD(K) = PSTD(K-1) * DLOGP
    end do
    PSTD(52)  = 0._r8
    PLMM1   = ETAA(LPAR)     !// pressure bottom of LPAR
    PTOPCTM = ETAA(LPAR+1)   !// pressure top of LPAR
    PTOPATM = 0._r8
    TOPMASS = (PTOPCTM - PTOPATM) * MASFAC

    M = max(1, min(12,JMON))
    do J = 1, JPAR
       N = max(1, min(18, (int(YDGRD(J)) + 99)/10 ))

       do K = 1, 51
          OREF2(K) = OREF(K,N,M)
          TREF2(K) = TREF(K,N,M)
       end do
       F0 = 0._r8
       T0 = 0._r8
       do K = 1, 51
          PC   = min(PTOPCTM,PSTD(K))
          PB   = max(PTOPATM,PSTD(K+1))
          if (PC .gt. PB) then
             XC = (PC - PB) / (PTOPCTM - PTOPATM)
             F0 = F0 + OREF2(K) * XC
             T0 = T0 + TREF2(K) * XC
          end if
       end do

       do I = 1, IPAR
          TOPT(I,J) = T0
          TOPM(I,J) = TOPMASS
          TOP3(I,J) = 1.e-6_r8*F0*TOPMASS
       end do

       !// Similar values for layer LPAR; between ETAA(LPAR) and PTOPCTM
       !// ONLY do this when chemistry is not calculated in LPAR.
       F0 = 0._r8
       T0 = 0._r8
       do K = 1,51
          PC   = min(PLMM1,PSTD(K))
          PB   = max(PTOPCTM,PSTD(K+1))
          if (PC .gt. PB) then
             XC = (PC - PB) / (PLMM1 - PTOPCTM) !Frac. of clim box inside CTM box
             F0 = F0 + OREF2(K) * XC   !Sum up contribution for O3
             T0 = T0 + TREF2(K) * XC   !Sum up contribution for T
          end if
       end do
       do I = 1, IPAR
          !// T in layer LPAR
          LPART(I,J) = T0
          !// Mass from bottom to top of LPAR
          LPARM(I,J) = (PLMM1 - PTOPCTM) * MASFAC
          !// O3 in layer LPAR
          LPAR3(I,J) = 1.e-6_r8 * F0 * (PLMM1 - PTOPCTM) * MASFAC
       end do

    end do


    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LPAR+1
             PJ(L) = ETAA(L) + ETAB(L) * P(I,J)
          end do
          do L = 1, LPAR
             DMS(L,I,J)  = (PJ(L) - PJ(L+1)) * MASFAC
          end do
       end do
    end do

    !//---load O3 from CTM is being calculated:
    do N = 1, NPAR
       if (TNAME(N) .eq. 'O3') then
          do J = 1, JPAR
             if (.not. LOSLOCSTRAT) then
                !// We do not calculate DO3 from climatology above.
                !// This can be changed back to CTM2 style, but for now we
                !// use LPAR, since O3 is set from CTM2-climatology in
                !// stratosphere.
                LMAX = LPAR !maxval(LMTROP(:,J))
             else
                LMAX = LPAR
             end if
             do I = 1, IPAR
                do L = 1, LMAX
                   ! unit of DO3:  # molecules/cm^2
                   !// Value in LPAR may be overwritten in p-phot_oc.f, by 
                   !// LPARO3.
                   DO3(L,I,J) = AVOGNR * 1.e-1_r8 * STT(I,J,L,N) &
                                / (48._r8 * AREAXY(I,J))
                end do
             end do
          end do
       end if
    end do

    !//---------------------------------------------------------------------
  end subroutine SET_ATM
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  !// For monthly solar variation, include the following in cmn_fjx.f90 and
  !// use these routines.
  !real(r8) :: FLmin(WX_), FLmax(WX_)
  !real(r8) :: FLscale(12,16) !// 12 months for 16 years (1997-2012)
  !//-----------------------------------------------------------------------
  !subroutine solar_flget()
  !  !//---------------------------------------------------------------------
  !  !// Call this every month to update FL.
  !  !// You need to make the proper changes.
  !  !// Amund Sovde Haslerud, December 2017
  !  !//---------------------------------------------------------------------
  !  use cmn_fjx, only: FL, FLmin, FLmax, FLscale, WX_
  !  !//---------------------------------------------------------------------
  !  implicit none
  !  !//---------------------------------------------------------------------
  !  real(r8) :: fracmax
  !  integer :: YY, K
  !  !//---------------------------------------------------------------------
  !
  !  !// YY=1 is 1997
  !  YY = JYEAR - 1996
  !  fracmax = FLscale(JMON,YY)
  !
  !  !// Interpolate between FLmax and FLmin
  !  do K = 1, WX_
  !     FL(K) = FLmin(K) * (1.d0 - fracmax) + FLmax(K) * fracmax
  !  end do
  !  write(*,'(a,i5,i3,f8.3)') 'FRACMAX',YY,JMON,fracmax
  !  !// Print the value
  !  write(*,'(a3,18es10.3)') 'FL:',FL
  !
  !  !//---------------------------------------------------------------------
  !end subroutine solar_flget
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
  !subroutine solar_flscale_read()
  !  !//---------------------------------------------------------------------
  !  !// Call this at model start to read the FLscale, i.e. fraction between
  !  !// FLmin and FLmax for each month.
  !  !// You need to make the proper changes.
  !  !// Amund Sovde Haslerud, December 2017
  !  !//---------------------------------------------------------------------
  !  use cmn_fjx, only: WX_, FLscale
  !  !//---------------------------------------------------------------------
  !  implicit none
  !  !//---------------------------------------------------------------------
  !  integer :: Y,M
  !  real(r8) :: frac
  !  !//---------------------------------------------------------------------
  !
  !  open(1,file='tables/solflux_7month_running_mean_1997-2012.dat', &
  !        form='formatted',status='old')
  !  print*,'solflux_7month_running_mean_1997-2012'
  !  do Y=1,16
  !     do M=1,12
  !        read(1,'(4x,1x,2x,1x,f12.10)') frac
  !        FLscale(M,Y) = frac
  !        print*,Y+1996,M,frac
  !     end do
  !  end do
  !  close(1)
  !
  !  !//---------------------------------------------------------------------
  !end subroutine solar_flscale_read
  !//-----------------------------------------------------------------------



  


  !//-----------------------------------------------------------------------
end module fastjx
!//=========================================================================
