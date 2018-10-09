!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Secondary organic aerosols (SOA) - Oslo treatment.
!//=========================================================================
module soa_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: soa_oslo
  !// DESCRIPTION: Module for SOA package.
  !//
  !// Contains SOA variables and routines.
  !//   subroutine soa_init
  !//   subroutine soa_setdrydep
  !//   subroutine SOA_v9_separate
  !//   subroutine soa_diag_soagas
  !//   subroutine FINDM0
  !//   subroutine FM0
  !//   subroutine soa_diag_drydep
  !//   subroutine soa_diag_separate
  !//   subroutine soa_nopsdiag
  !//   subroutine soa_diag_lsscav
  !//
  !// This module also contains some comprehensive diagnostics for SOA,
  !// to keep track of what is produced and lost of all SOA.
  !//
  !// Reference often used for SOA:
  !//   Chung and Seinfeld, JGR 2002, doi:10.1029/2001JD001397
  !//
  !// Amund Sovde, June 2012
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: LSOA, IPAR, JPAR, LPAR, NPAR, NPAR_SOA, &
       IDBLK, JDBLK, MPBLK
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// SOA component indices
  !// SOA gas phase:  150-164
  !// SOA (aerosols): 165-191
  !// SOA precursors: 192-193 and 280-291

  !// Variables
  !// List of all transport numbers used for SOA components
  integer :: soa_trsp_idx(NPAR_SOA)

  !// List of all SOAs with drydeps (i.e. all the aerosols)
  integer,parameter :: ndep_soa = 21
  integer,dimension(ndep_soa),parameter   :: soa_deps = &
       (/ 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, &
          175, 176, 177, 178, 179, 182, 183, 188, 189, 190, &
          191 /)
  real(r8) :: VSOA(ndep_soa)
  !// List of all SOA gases (not currently used)
  integer,parameter :: ngas_soa = 21
  integer,dimension(ngas_soa),parameter   :: soa_gas = &
       (/ 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, &
          160, 161, 162, 163, 164, 180, 181, 184, 185, 186, &
          187 /)

  !// SOA accounting/budgets (in units of grams [g])
  !// The following processes are assessed
  !//   - Production through the routine SOA_v9_separate (PROD)
  !//   - Evaporation through the routine SOA_v9_separate (EVAP)
  !//   - Dry deposition (DDEP)
  !//   - Large scale wet scavenging (LS)
  !//   - Convective wet scavenging (CNV)
  !// The last two are done by the model core diagnosics, and
  !// do not need separate arrays. Here we define the necessary
  !// arrays to diagnose the rest:
  real(r8), dimension(IPAR,JPAR,NDEP_SOA) :: TOTAL_SOA_DDEP
  real(r8), dimension(IPAR,JPAR,LPAR,NDEP_SOA) :: TOTAL_SOA_PROD
  real(r8), dimension(IPAR,JPAR,LPAR,NDEP_SOA) :: TOTAL_SOA_EVAP
  real(r8), dimension(IPAR,JPAR,LPAR,NDEP_SOA) :: TOTAL_SOA_LSSCAV
  real(r8), dimension(NDEP_SOA) :: oldProd,oldEvap,oldDDep,oldLscv,oldCscv

  character(len=80),parameter :: filesoabudget = 'soa_budgets.dta'
  character(len=80),parameter :: filesoagbudget = 'soag_budgets.dta'
  !// Also save daily totals throughout the year
  real(r8), dimension(NDEP_SOA,366) :: &
       dailyProd,dailyEvap,dailyDDep,dailyLscv,dailyCscv,dailyBrdn

  !// Diagnose SOAGAS
  real(r8), dimension(IPAR,JPAR,LPAR,NGAS_SOA) :: TOTAL_SOAG_PROD
  real(r8), dimension(IPAR,JPAR,LPAR,NGAS_SOA) :: TOTAL_SOAG_SEPPROD
  real(r8), dimension(IPAR,JPAR,LPAR,NGAS_SOA) :: TOTAL_SOAG_SEPCOND
  real(r8), dimension(IPAR,JPAR,LPAR,NGAS_SOA) :: TOTAL_SOAG_LSSCAV
  real(r8), dimension(NGAS_SOA) :: oldProdGas,oldLscvGas,oldCscvGas, &
       oldSepProdGas, oldSepCondGas

  real(r8), dimension(NGAS_SOA,366) :: &
       dailyProdGas, dailyLscvGas, dailyCscvGas, dailySepProdGas, &
       dailySepCondGas, dailyBrdnGas
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file='soa_oslo.f90'
  !// ----------------------------------------------------------------------
  save !// All variables are to be saved.
  private
  public soa_init, soa_setdrydep, SOA_v9_separate, soa_diag_drydep, &
       soa_diag_separate, soa_diag2file, soa_nopsdiag, &
       ndep_soa, soa_deps, soa_diag_lsscav, soa_diag_soagas
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine soa_init()
    !// --------------------------------------------------------------------
    !// Initialize SOA simulations.
    !//
    !// Amund Sovde, June 2012
    !// --------------------------------------------------------------------
    use cmn_chem, only: TNAME
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// For file stuff
    integer :: ifnr
    logical :: fnr_ok 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'soa_init'
    !// --------------------------------------------------------------------

    !// Check for name of component id 12
    if (LSOA .and. TNAME(trsp_idx(12)).ne.'C6HXR_SOA') then
       write(6,'(a)') f90file//':'//subr//': SOA Secondary organic aerosols:'
       write(6,'(a)') '  Tracer 12 is wrongly named: '//trim(TNAME(trsp_idx(12)))
       write(6,'(a)') '  It must be named C6HXR_SOA'
       write(6,'(a)') '  This is because the original C6HXR is split up in SOA'
       stop
    end if


    !// Dry deposition of SOA
    !// Assuming that all aerosol species sediment with velocity 0.1cm/s,
    !// as in Chung and Seinfeld JGR 2002.
    VSOA(:) = 0.1_r8

    !// Diagnoses
    TOTAL_SOA_DDEP(:,:,:) = 0._r8
    TOTAL_SOA_PROD(:,:,:,:) = 0._r8
    TOTAL_SOA_EVAP(:,:,:,:) = 0._r8
    TOTAL_SOA_LSSCAV(:,:,:,:) = 0._r8

    TOTAL_SOAG_PROD(:,:,:,:) = 0._r8
    TOTAL_SOAG_SEPPROD(:,:,:,:) = 0._r8
    TOTAL_SOAG_SEPCOND(:,:,:,:) = 0._r8
    TOTAL_SOAG_LSSCAV(:,:,:,:) = 0._r8

    !// For NOPS diags
    oldProd(:) = 0._r8
    oldEvap(:) = 0._r8
    oldDDep(:) = 0._r8
    oldLscv(:) = 0._r8
    oldCscv(:) = 0._r8

    oldProdGas(:) = 0._r8
    oldLscvGas(:) = 0._r8
    oldCscvGas(:) = 0._r8
    oldSepProdGas(:) = 0._r8
    oldSepCondGas(:) = 0._r8

    !// For daily diags
    dailyProd(:,:) = 0._r8
    dailyEvap(:,:) = 0._r8
    dailyDDep(:,:) = 0._r8
    dailyLscv(:,:) = 0._r8
    dailyCscv(:,:) = 0._r8
    dailyBrdn(:,:) = 0._r8

    dailyProdGas(:,:) = 0._r8
    dailyLscvGas(:,:) = 0._r8
    dailyCscvGas(:,:) = 0._r8
    dailySepProdGas(:,:) = 0._r8
    dailySepCondGas(:,:) = 0._r8
    dailyBrdnGas(:,:) = 0._r8

    !// Open result file
    !// Find non-used file number for input file
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do
    !// Initialize output file
    write(6,'(a)') f90file//':'//subr// &
         ': Initializing SOA budget file '//trim(filesoabudget)
    open(ifnr,file=filesoabudget,form='unformatted')
    write(ifnr) IPAR,JPAR,LPAR,NDEP_SOA
    write(ifnr) soa_deps
    close(ifnr)

    write(6,'(a)') f90file//':'//subr// &
         ': Initializing SOAG budget file '//trim(filesoagbudget)
    open(ifnr,file=filesoagbudget,form='unformatted')
    write(ifnr) IPAR,JPAR,LPAR,NGAS_SOA
    write(ifnr) soa_gas
    close(ifnr)

    !// --------------------------------------------------------------------
  end subroutine soa_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine soa_setdrydep(VDEP,MP)
    !// --------------------------------------------------------------------
    !// Set dry deposition for SOA, in a given IJ-block.
    !// Called from subroutine setdrydep, which also modifies deposition
    !// rate by stability.
    !//
    !// Amund Sovde, June 2012
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    !// Input/Output
    real(r8), intent(inout) :: VDEP(NPAR,IPAR,JPAR)

    !// Locals
    integer :: I,II,J,JJ,N, K
    !// --------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          do K = 1, ndep_soa

            !// Get transport number
            N = trsp_idx(soa_deps(K))
            if (N .gt. 0) then
               !// Convert from cm/s to m/s
               VDEP(N,I,J) = VSOA(K)*1.e-2_r8
            end if

          end do !// do K = 1, ndep_soa
       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine soa_setdrydep
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine SOA_v9_separate( &
       ZC_MASS, & !// the mass of each species
       TEMP, &    !// temperature in column
       DV, &      !// volume in column
       ILOC, &    !// I-box
       JLOC, &    !// J-box
       LSTART, &  !// Start level for calculations
       LEND)      !// End level for calculations
    !// --------------------------------------------------------------------
    use cmn_size, only: TRACER_ID_MAX
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    !// Restructured for CTM3
    !//   Treats a column, initially called after oslo_chem.
    !//   Should consider calling it in a separate aerosol master routine.
    !//
    !// Amund Sovde, June 2012
    !// --------------------------------------------------------------------
    !//
    !// --------------------------------------------------------------------
    !// The separate bit stands for the fact that in this version the precursor
    !// hydrocarbons are treated separately before oxidation.
    !// --------------------------------------------------------------------
    !// History:
    !//  CRH began developing this code on Friday, 13 Jan 2006.
    !//  V3 came from Terje, and was not changed at all by me.
    !//  V4 was the first test version of my developed code
    !//  V5 is a version to be adapted to be used by the model. (2/5/06)
    !//     to set up to run with real data for each class look for "CRH ms"
    !//     to get rid of any testing stuff, look for "!CRH test"
    !//  V6 Im shifting all the hydrocarbon oxidation to the tropospheric
    !//     chemistry solver, it dosn't seem to work here. From now on all
    !//     this subroutine will do is partition between gas and aerosol phase.
    !//     Another way of doing it would have been to record the amount of
    !//     hydrocarbon lost via oxidation from each oxidant in the chemistry
    !//     code and send it here as deltaHC.  
    !//  V7 (15/6/06) I added isoprene oxidation products today. look for
    !//     "v7" All coefficients and stuff are from Henze and Seinfeld GRL2006
    !//  V8 (7/7/2006) I added the oxidation products of aromatics today.
    !//     (9/8/2006) I didnt think it was worth a new version, but I am going 
    !//     to remove the aerosol phase oxidation products from AR1 and ISOR1.
    !//     If aerosol evaporates, I'll put them back. This will lead to loss
    !//     of AR1 and ISOR1 via processes affecting aerosol, which makes sense,
    !//     since we are assuming that a fraction of AR1 and ISOR1 is capable
    !//     of forming aerosol. THAT GOT REMOVED AGAIN!!
    !//  V9 (24.4.2007) Will allow aerosol to partition on ammonium sulphate 
    !//     aerosol, and will prevent evaporation of aerosol. Benzene may
    !//     be added too. Actually I'll just lump the benzene in with toluene.
    !//     changes marked with "V9"
    !//
    !//
    !// Hm, it seems that nowhere in this code is any story about a bear.
    !// In particular not about a cheese bear. That is odd! How can CRH 
    !// expect that a program works if it does not contain a bear story? 
    !// Bears don't know much about modelling, but lots about bears.
    !// And cooking. And baking. And all these things are essential to a model. 
    !// So therefore, there will soon be some cooking to keep CHRs brain in
    !// proper function. Good bye, the bear.
    !//
    !// More history from Amund Sovde, developing CTM3:
    !// - Not sure what CRH meant by preventing evaporation in V9.
    !// - Isoprene + OH gives ISOR1 and SOAGAS61+SOAGAS62. The fraction
    !//   going to SOAGAS is about 25%, so it does not make sense to keep
    !//   the full path to ISOR1. It will in turn give too much AR1.
    !//   Therefore, the production of ISOR1 is reduced accordingly.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in) :: ILOC, JLOC, LSTART, LEND
    real(r8), intent(in) :: TEMP(LPAR), DV(LPAR)

    !// In/Out
    real(r8), intent(inout) :: ZC_MASS(TRACER_ID_MAX,LPAR)


    !// Number of precursor hydrocarbon classes
    integer, parameter :: NHC = 8
    !// Number of oxidants
    integer, parameter :: NOX = 3
    !// Number of products
    integer, parameter :: NP = 2


    !// species concentrations at the grid box in question:
    REAL(r8) :: MUMASS_CONC(TRACER_ID_MAX)

  
    !// Define all the different variables:
    integer :: L,ID,N1,N2,N3, TR1, TR2, TR3
    real(r8)  :: X1,SG, TEM

    !// Mass concentrations
    real(r8) :: &
         GC(NHC,NOX,NP), &  !// Gas-phase conc. of HC oxidation products
         GC0(NHC,NOX,NP), & !// Gas-phase conc. at the beginning of timestep
         A(NHC,NOX,NP), &   !// Aerosol-phase conc. of HC oxidation products
         A0(NHC,NOX,NP)     !// Aerosol-phase conc. at beginning of timestep
    real(r8) :: SOAM, POA, SOA0, M0

    !// Temperature dependence factor, for product classes
    real(r8) :: TFAC(NHC)

    !// Equlibirum array
    real(r8) :: KOM(NHC,NOX,NP)

    !// Equilibrium data from Chung and Seinfeld (JGR, 2002)
    real(r8), parameter :: &
         !// Products of O3 oxidation
         K0_111 = 0.184_r8, K0_211 = 0.055_r8, &
         K0_311 = 0.133_r8, K0_411 = 0.224_r8, &
         K0_511 = 0.0459_r8, &
         K0_112 = 0.0043_r8, K0_212 = 0.0053_r8, &
         K0_312 = 0.0035_r8, K0_412 = 0.0082_r8, &
         K0_512 = 0._r8, &
         !// Products of OH oxidation (equal to O3)
         K0_121 = 0.184_r8, K0_221 = 0.055_r8, &
         K0_321 = 0.133_r8, K0_421 = 0.224_r8, &
         K0_521 = 0.0459_r8, &
         K0_122 = 0.0043_r8, K0_222 = 0.0053_r8, &
         K0_322 = 0.0035_r8, K0_422 = 0.0082_r8, &
         K0_522 = 0._r8, &
         !// Products of NO3 reactions
         K0_131 = 0.0163_r8, K0_231 = 0.0163_r8, &
         K0_331 = 0.0163_r8, K0_431 = 0.0163_r8, &
         K0_531 = 0.0163_r8, &
         K0_132 = 0._r8,    K0_232 = 0._r8, &
         K0_332 = 0._r8,    K0_432 = 0._r8, &
         K0_532 = 0._r8, &
         !// The rest (no NO3 reaction)
         !// partitioning coeff from Henz and seinfeld GRL 2006 
         K0_611 = 0._r8,      K0_612 = 0._r8, &
         K0_621 = 0.00862_r8, K0_622 = 1.62_r8, &
         K0_63x = 0._r8, &
         !// partitioning coeff from Odum et al Science (1997?) v8 :
         K0_711 = 0.042_r8, K0_712 = 0.0014_r8, &
         K0_721 = 0.042_r8, K0_722 = 0.0014_r8, &
         K0_73x = 0._r8, &
         K0_811 = 0.053_r8, K0_812 = 0.0019_r8, &
         K0_821 = 0.053_r8, K0_822 = 0.0019_r8, &
         K0_83x = 0._r8


   
    !// dH/R, to calculate temp. dep. of the partion coefficient
    !// (chung and seinfeld JGR 2002)
    real(r8), parameter :: DHR = 5.0e3_r8   !// (K)
    logical, parameter :: LDEBUG_SEP = .false.
    !// --------------------------------------------------------------------

    !// Initialize
    MUMASS_CONC(:) = 0._r8

    !// Loop over L
    do L = LSTART, LEND

       TEM = TEMP(L)

       !// Temperature dependence factor, for classes 1-5, 6 and 7-8
       TFAC(1:5) = TEM/298._r8 * EXP(DHR*(1._r8/TEM - 1._r8/298._r8))
       TFAC(6)   = TEM/295._r8 * EXP(DHR*(1._r8/TEM - 1._r8/295._r8))
       TFAC(7:8) = TEM/298._r8 * EXP(DHR*(1._r8/TEM - 1._r8/298._r8))

       !// Set up array with reaction rates
       KOM(:,1,1) = (/K0_111,K0_211,K0_311,K0_411,K0_511,K0_611,K0_711,K0_811/)
       KOM(:,1,2) = (/K0_112,K0_212,K0_312,K0_412,K0_512,K0_612,K0_712,K0_812/)
       KOM(:,2,1) = (/K0_121,K0_221,K0_321,K0_421,K0_521,K0_621,K0_721,K0_821/)
       KOM(:,2,2) = (/K0_122,K0_222,K0_322,K0_422,K0_522,K0_622,K0_722,K0_822/)
       KOM(:,3,1) = (/K0_131,K0_231,K0_331,K0_431,K0_531,K0_63x,K0_73x,K0_83x/)
       KOM(:,3,2) = (/K0_132,K0_232,K0_332,K0_432,K0_532,K0_63x,K0_73x,K0_83x/)

       !// Modify by TFAC
       do N1 = 1, NHC
          KOM(N1,:,:) = KOM(N1,:,:) * TFAC(N1)
       end do

       !// According to Zhang et al., PNAS, doi:10.1073/pnas.1404727111
       !// (2014), the conversion to aerosols should be increased by
       !// up to 4, so we implement this.
       KOM(:,:,:) = KOM(:,:,:) * 4._r8

       !// Debug check
       if (LDEBUG_SEP) then
         do ID = 150, 191
           if ((ZC_MASS(ID,L) .ne. ZC_MASS(ID,L)) .or. &
               (ZC_MASS(ID,L) .lt. 0._r8) ) then
             if (abs(ZC_MASS(ID,L)) .lt. 1.e-30_r8) then
               print*,'small value in SOA, setting to 0 ',ZC_MASS(ID,L),ID,L
               ZC_MASS(ID,L) = 0._r8
             else
               print*,'dodgy mass at first point,'
               print*,'soa_oslo.f90:',ZC_MASS(ID,L),ID,L
               stop
             end if
           end if
         end do
       end if !// if (LDEBUG_SEP) then


       !// Calculations are done in ug/m3, need to convert from kg/gridbox
       do ID = 150, 191
          !// kg -> ug: 1.d9
          MUMASS_CONC(ID) = ZC_MASS(ID,L) * 1.e9_r8 / DV(L)
       end do

       !// Also for some additional species OC and SO4:
       !// POA Primary organic aerosols
       !// Hydrophobic and hydrophilic, from biomass burning
       !// and fossil fuel contribution.
       !//
       !// Component 230-237 are already organic matter. Previously, they were
       !// OC, and had to be multiplied by 2.6 for biomass burning species
       !// and 1.6 for fossil fuel species. That is no longer necessary.
       !// For SO4, use mass as (NH4)2SO4 132.14/96
       POA = ( ZC_MASS(230,L)   & !// omBB1fob
               + ZC_MASS(231,L) & !// omBB1fil
               + ZC_MASS(232,L) & !// omFF1fob
               + ZC_MASS(233,L) & !// omFF1fil
               + ZC_MASS(234,L) & !// omBF1fob
               + ZC_MASS(235,L) & !// omBF1fil
               + ZC_MASS(236,L) & !// omOCNfob
               + ZC_MASS(237,L) & !// omOCNfil
               + ZC_MASS( 73,L) * 132.14_r8/96._r8 &  !// (NH4)2SO4
             ) * 1.e9_r8 / DV(L) !// -> ug/m3



       !// Gas phase products
       !// The products behave (are?) the same for oxidation by O3 and OH.
       !// Since we dont want to transport identical tracers, we transport
       !// the products lumped together. Then we split them up into the
       !// two different compounds here again, half each (maybe not necessary,
       !// but probably more flexible).
       !//
       !// The transported tracers are one for each class, and two for the
       !// products of oxidation by both O3 and OH (two identical products
       !// from each oxidation reaction), i.e.
       !//   O3 + HC -> product1 + product2
       !//   OH + HC -> product1 + product2 (the same product1 and product2)
       !// There is one tracer for the products of oxidation by NO3
       !// and only one product from its oxidation reaction.

       !// Originally, the tracers were distributed in a double loop,
       !// more specifically for gas phase as:
       !//ID = 150
       !//do N2 = 1, 3
       !//   !// Loop over classes
       !//   do N1 = 1,5 !v7 changed NHC to 5
       !//      if (N2 .eq. 1) then
       !//         !// Distribute tracer on oxidants O3 and OH
       !//         GC0(N1,1,1) = 0.5d0 * MUMASS_CONC(ID)
       !//         GC0(N1,2,1) = MUMASS_CONC(ID) - GC0(N1,1,1)
       !//      else if (N2 .eq. 2) then
       !//         GC0(N1,1,2) = 0.5d0 * MUMASS_CONC(ID)
       !//         GC0(N1,2,2) = MUMASS_CONC(ID) - GC0(N1,1,2)
       !//      else if (N2 .eq. 3) then
       !//         GC0(N1,N2,2) = 0.d0  ! only one product for NO3 oxidn.
       !//         GC0(N1,N2,1) = MUMASS_CONC(ID)
       !//      end if
       !//      ID = ID + 1
       !//   end do
       !//end do
       !//
       !// But this loop is a bit obscure. What is done is that
       !// the gas phase species of each class are divided into two
       !// species. A simpler loop, although also a bit difficult to
       !// grasp, is:
       do N1 = 1, 5
          !// Tracer ID for 5 species
          TR1 = 150 + N1 - 1 !// Tracer for oxidants 1&2, product 1
          TR2 = 155 + N1 - 1 !// Tracer for oxidants 1&2, product 2
          TR3 = 160 + N1 - 1 !// Oxidant 3, product 1&2

          !// Separate TR1 into O3 and OH oxidants, product 1
          GC0(N1,1,1) = 0.5_r8 * MUMASS_CONC(TR1)
          GC0(N1,2,1) = MUMASS_CONC(TR1) - GC0(N1,1,1)

          !// Separate TR2 into O3 and OH oxidants, product 2
          GC0(N1,1,2) = 0.5_r8 * MUMASS_CONC(TR2)
          GC0(N1,2,2) = MUMASS_CONC(TR2) - GC0(N1,1,2)

          !// Separate TR3 into NO3 oxidant, product 1 only
          GC0(N1,3,1) = MUMASS_CONC(TR3)
          GC0(N1,3,2) = 0._r8
       end do
       !// Continue with the rest of the classes:
       !// isoprene OH products in the gas phase:
       GC0(6,1,:) = 0._r8
       GC0(6,2,1) = MUMASS_CONC(180)
       GC0(6,2,2) = MUMASS_CONC(181)
       GC0(6,3,:) = 0._r8
       !// aromatic oh and o3 products in the gas phase:
       GC0(7,1,1) = MUMASS_CONC(184) * 0.5_r8
       GC0(7,2,1) = MUMASS_CONC(184) - GC0(7,1,1)
       GC0(7,1,2) = MUMASS_CONC(185) * 0.5_r8
       GC0(7,2,2) = MUMASS_CONC(185) - GC0(7,1,2)
       GC0(7,3,:) = 0._r8
       GC0(8,1,1) = MUMASS_CONC(186) * 0.5_r8
       GC0(8,2,1) = MUMASS_CONC(186) - GC0(8,1,1)
       GC0(8,1,2) = MUMASS_CONC(187) * 0.5_r8
       GC0(8,2,2) = MUMASS_CONC(187) - GC0(8,1,2)
       GC0(8,3,:) = 0._r8

       !// Old aerosol loop
       !//ID = 165
       !//do N2 = 1, 3
       !//   !// Loop over classes
       !//   do N1 = 1, 5
       !//      if (N2 .eq. 1) then
       !//         A0(N1,1,1) = 0.5d0 * MUMASS_CONC(ID)
       !//         A0(N1,2,1) = MUMASS_CONC(ID) - A0(N1,1,1)
       !//      else if (N2 .eq. 2) then
       !//         A0(N1,1,2) = 0.5d0 * MUMASS_CONC(ID)
       !//         A0(N1,2,2) = MUMASS_CONC(ID) - A0(N1,1,2)
       !//      else if (N2 .eq. 3) then 
       !//         A0(N1,N2,2) = 0.d0  ! only one product for NO3 oxidn.
       !//         A0(N1,N2,1) = MUMASS_CONC(ID)
       !//      end if
       !//      ID = ID + 1
       !//   end do
       !//end do

       do N1 = 1, 5
          !// Tracer ID for 5 species
          TR1 = 165 + N1 - 1 !// Aerosol from oxidants 1&2, product 1
          TR2 = 170 + N1 - 1 !// Aerosol from oxidants 1&2, product 2
          TR3 = 175 + N1 - 1 !// Aerosol from oxidant 3, product 1

          A0(N1,1,1) = 0.5_r8 * MUMASS_CONC(TR1)
          A0(N1,2,1) = MUMASS_CONC(TR1) - A0(N1,1,1)

          A0(N1,1,2) = 0.5_r8 * MUMASS_CONC(TR2)
          A0(N1,2,2) = MUMASS_CONC(TR2) - A0(N1,1,2)

          A0(N1,3,1) = MUMASS_CONC(TR3)
          A0(N1,3,2) = 0._r8
       end do

       !// isoprene OH products in the aerosol phase:
       A0(6,1,:) = 0._r8
       A0(6,2,1) = MUMASS_CONC(182)
       A0(6,2,2) = MUMASS_CONC(183)
       A0(6,3,:) = 0._r8

       !// aromatic oh and o3 products in the aerosol phase:
       A0(7,1,1) = MUMASS_CONC(188) * 0.5_r8
       A0(7,2,1) = MUMASS_CONC(188) - A0(7,1,1)
       A0(7,1,2) = MUMASS_CONC(189) * 0.5_r8
       A0(7,2,2) = MUMASS_CONC(189) - A0(7,1,2)
       A0(7,3,:) = 0._r8

       A0(8,1,1) = MUMASS_CONC(190) * 0.5_r8
       A0(8,2,1) = MUMASS_CONC(190) - A0(8,1,1)
       A0(8,1,2) = MUMASS_CONC(191) * 0.5_r8
       A0(8,2,2) = MUMASS_CONC(191) - A0(8,1,2)
       A0(8,3,:) = 0._r8

       !// Total A0 mass concentration
       SOA0 = sum(A0)
       !// Total aerosol mass concentration
       M0 = POA + SOA0

       if (M0 .lt. 0._r8) then
          print*,'soa_oslo.f90: m0 less than 0, at begining of SOA code'
          stop
       end if

       !// Initialize new mass concentrations
       GC(:,:,:) = 0._r8
       A(:,:,:)  = 0._r8


       !// New total aerosol concentration (M0), Iterative procedure to solve
       !// M0 from eq.5
       !// -----------------------------------------------------------------
       call findM0(NHC,NOX,NP,KOM,A0,GC0,POA,SOA0,M0)

       do N3 = 1, NP
          do N2 = 1, NOX
             do N1 = 1, NHC
                !// CRH I removed the AST*dhc term below, as the newly reacted
                !// HCs are already accounted for in GC0, - they were added
                !// in OSLO_CHEM:
                X1 = KOM(N1,N2,N3) * M0 * (A0(N1,N2,N3) + GC0(N1,N2,N3))

                A(N1,N2,N3) = X1 / (1._r8 + KOM(N1,N2,N3) * M0)

                if ((KOM(N1,N2,N3) * M0) .gt. 0._r8) then
                   GC(N1,N2,N3) = A(N1,N2,N3) / (KOM(N1,N2,N3) * M0)
                else if (M0 .eq. 0._r8) then
                   !// if M0=0.d0, then A=0.d0 above, and all should be
                   !// in gas phase
                   !CRH hope this is the right thing to do!!
                   GC(N1,N2,N3) = A0(N1,N2,N3) + GC0(N1,N2,N3)
                else
                   !// KOM is zero (and possibly M0 also). Do nothing,
                   !// as A and GC are already initialized as zero.
                end if

             end do
          end do
       end do

       !// Useful diagnostics (?), but not used at the moment:
       SOAM = sum(A)
       SG = sum(GC)


       !// Putting back into the transport variables
       !// -----------------------------------------
       !// Gas phase
       do N1 = 1, 5
          !// Tracer ID for 5 species
          TR1 = 150 + N1 - 1 !// Tracer for oxidants 1&2, product 1
          TR2 = 155 + N1 - 1 !// Tracer for oxidants 1&2, product 2
          TR3 = 160 + N1 - 1 !// Oxidant 3, product 1&2
          MUMASS_CONC(TR1) = GC(N1,1,1) + GC(N1,2,1)
          MUMASS_CONC(TR2) = GC(N1,1,2) + GC(N1,2,2)
          MUMASS_CONC(TR3) = GC(N1,3,1)
       end do

       !// isoprene OH products in the gas phase:
       MUMASS_CONC(180) = GC(6,2,1)
       MUMASS_CONC(181) = GC(6,2,2)

       !// aromatic oh and o3 products in the gas phase:
       MUMASS_CONC(184) = GC(7,2,1) + GC(7,1,1)
       MUMASS_CONC(185) = GC(7,2,2) + GC(7,1,2)
       MUMASS_CONC(186) = GC(8,2,1) + GC(8,1,1)
       MUMASS_CONC(187) = GC(8,2,2) + GC(8,1,2)


       !// Aerosol phase
       do N1 = 1, 5
          TR1 = 165 + N1 - 1 !// Aerosol from oxidants 1&2, product 1
          TR2 = 170 + N1 - 1 !// Aerosol from oxidants 1&2, product 2
          TR3 = 175 + N1 - 1 !// Aerosol from oxidant 3, product 1
          MUMASS_CONC(TR1) = A(N1,1,1) + A(N1,2,1)
          MUMASS_CONC(TR2) = A(N1,1,2) + A(N1,2,2)
          MUMASS_CONC(TR3) = A(N1,3,1)
       end do

       !// isoprene OH products in the aerosol phase:
       MUMASS_CONC(182) = A(6,2,1)
       MUMASS_CONC(183) = A(6,2,2)

       !// aromatic oh and o3 products in the aerosol phase:
       MUMASS_CONC(188) = A(7,2,1) + A(7,1,1)
       MUMASS_CONC(189) = A(7,2,2) + A(7,1,2)
       MUMASS_CONC(190) = A(8,2,1) + A(8,1,1)
       MUMASS_CONC(191) = A(8,2,2) + A(8,1,2)


       !// Convert from ug/m3 to kg/gridbox
       do ID = 150, 191
          ZC_MASS(ID,L) = max(MUMASS_CONC(ID) * DV(L) * 1.e-9_r8, 1.e-20_r8)
       end do
    end do !// do L = LSTART, LEND

    !// --------------------------------------------------------------------
  end subroutine SOA_v9_separate
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine FINDM0(NHC,NOX,NP,KOM,A0,GC0,POA,SOA0,M0)
    !// --------------------------------------------------------------------
    !// Subroutine to find M0 (Concentration of total organic aerosol) from
    !// equation 5 in Chung and Seinfeld (2002). Equation 5 is an
    !// NHC*NOX*NP order equation in M0, and can thus not be solved
    !// analytically. This subroutine solves this equation by and iterative
    !// method.
    !// TKB March 2005
    !//
    !// To CTM3
    !// Amund Sovde, June 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NHC, NOX, NP
    real(r8), intent(in) :: &
         KOM(NHC,NOX,NP), & !// Equilibrium constants
         A0(NHC,NOX,NP), &  !// Aerosol concentration
         GC0(NHC,NOX,NP)    !// Gas phase concentration SOA
    real(r8), intent(in) :: &
         POA, &             !// Total POA mass concentration
         SOA0               !// Total SOA mass concentration before iteration

    !// Input/Output
    real(r8), intent(inout) :: M0

    !// Locals for iteration
    integer :: IT,N1,N2,N3
    real(r8) :: MM0, MX, MX1, MX2, RFM0
    real(r8) :: XP(NHC,NOX,NP)
    logical :: LITERATE, LITERATE2

    !// Parameters for iteration
    real(r8), parameter :: &
         ddm = 0.1_r8, &
         xlim = 1.e-3_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'FINDM0'
    !// --------------------------------------------------------------------

    !// New total aerosol concentration (M0), Iterative procedure to solve
    !// M0 from eq.5

    !// Calculate first products in the nominator of the first term of the
    !// expression on the left hand side of eq. 5
    !// There is no dhc in the line below, because its already included in
    !// GC0. This was done in OSLO_CHEM
    !do N3 = 1, NP
    !   do N2 = 1, NOX
    !      do N1 = 1, NHC
    !         XP(N1,N2,N3) = KOM(N1,N2,N3) * (A0(N1,N2,N3) + GC0(N1,N2,N3))
    !      end do
    !   end do
    !end do
    XP(:,:,:) = KOM(:,:,:) * (A0(:,:,:) + GC0(:,:,:))

    !// Aerosol can only form if there is existing POA or SOA, 
    !// i.e. no heterogeneous nucleation.
    !// Note that SUM(XP) must be greater than one if POA=0.
    !// Note also that after spinup, theres nowhere where POA is
    !// likely to be zero.
    LITERATE = ((sum(XP).ge.1._r8) .and. (SOA0.gt.0._r8)) .or. (POA.gt.0._r8)
    if (.not.LITERATE) then
       M0 = 0._r8
       return
    end if

    !// Initialize additional iteration flag used below.
    !// Instead of setting this to .false., I think we could rather
    !// test if M0==0.
    LITERATE2 = .true.


    !// Will do iteration
    !// --------------------------------------------------------------------

    !// RFM0= Left-hand side of eq. 5 (C.&S., 2002)
    CALL FM0(NHC,NOX,NP,KOM,XP,POA,M0,RFM0)

    !// Iteration counter
    IT = 0

    if (RFM0 .gt. 1._r8) then
       do while (RFM0 .gt. 1._r8)
          IT = IT + 1
          M0 = M0 * (1._r8 + ddm)

          call FM0(NHC,NOX,NP,KOM,XP,POA,M0,RFM0)

          !// Exit loop if M0=0
          if (M0 .eq. 0._r8) then
             !// Original CTM2 code seem to do the second iteration below,
             !// but actually at this point its RFM0 was either Infinity
             !// or NaN, since routine FM0 divides by M0.
             !// Thus the second iteration was not carried out.
             !// CTM3 does not divide by zero M0, but we skip iteration 2.
             write(6,'(a,i7,a)') f90file//':'//subr// &
                  ': RFM0>1: M0=0 after ',IT,' iterations; returning'
             LITERATE2 = .false.
             exit !// exit do while loop
          end if

          if (IT .gt. 10000) then
             write(6,'(a,i7,a)') f90file//':'//subr// &
                  ': Too many iterations in FINDM0, giving up (1)'
             write(6,'(a,4es16.6)') f90file//':'//subr// &
                  ': RFM0,M0,SUM(XP),POA=',RFM0,M0,sum(XP),POA
             stop
          end if
       end do !// do while (RFM0 .gt. 1._r8)

       MX = M0
       MM0 = M0 / (1._r8 + ddm)

    else if (RFM0 .lt. 1._r8) then
       do while (RFM0 .lt. 1._r8)
          IT = IT + 1
          M0 = M0 * (1._r8 - ddm)
          if (M0 .lt. 1.e-30_r8) then
             write(6,'(a,i7,a)') f90file//':'//subr//': setting M0 to 0'
             M0 = 0._r8
          end if

          call FM0(NHC,NOX,NP,KOM,XP,POA,M0,RFM0)

          if (M0 .eq. 0._r8) then
             !// Exit this iteration.
             write(6,'(a,i7,a)') f90file//':'//subr// &
                  ': RFM0<1: M0=0 after ',IT,' iterations; returning'
             !// Original CTM2 code leaves M0 as zero and skips the second
             !// iteration below. In CTM2, RFM0 was now either Infinity or
             !// NaN, but since we exit the loop and iteration 2 is not done,
             !// RFM0 is no longer required.
             !// FM0 in CTM3 will return a high value of RFM0 when
             !// POA>0 and M0==0. We still skip iteration 2.
             LITERATE2 = .false.
             exit !// exit do while loop
          end if

          if (IT .gt. 10000) then
             write(6,'(a,i7,a)') f90file//':'//subr// &
                  ': Too many iterations in FINDM0, giving up (2)'
             write(6,'(a,4es16.6)') f90file//':'//subr// &
                  ': RFM0,M0,SUM(XP),POA=',RFM0,M0,sum(XP),POA
             stop
          end if
       end do !// do while (RFM0 .lt. 1._r8)
    
       MX = M0
       MM0 = M0
       M0 = M0 / (1._r8 - ddm)
    else if (RFM0 .eq. 1._r8) then
       !// When RFM0=1, iteration 2 (do while loop) will not happen.
       LITERATE2 = .false.
       !// Do nothing to M0, it will remain unchanged
       MM0 = 0._r8 !// Not necessary because second iteration is not done
       MX = M0    !// Not necessary because second iteration is not done
    end if

    !// If first iteration gave M0=0, LITERATE2 was set to false,
    !// to stop further iteration.
    if (LITERATE2) then
       !// M0 is not zero, iterate again
       MX1 = MM0
       MX2 = M0

       !// Reset iteration counter
       IT = 0
       !Print*,'RFM0,M0,SUM(XP),POA before 3=',RFM0,M0,SUM(XP),POA
       do while (abs(1._r8 - RFM0) .gt. XLIM) 
          !//MX = MX1 + (MX2 - MX1) * 0.5_r8
          MX = (MX1 + MX2) * 0.5_r8
          IT = IT + 1

          call FM0(NHC,NOX,NP,KOM,XP,POA,MX,RFM0)

          if (RFM0 .GT. 1._r8) then
             MX1 = MX
          else
             MX2 = MX
          end if
          if (IT .gt. 10000) then
             write(6,'(a,i7,a)') f90file//':'//subr// &
                  ': Too many iterations in FINDM0, giving up (3)'
             write(6,'(a,4es16.6)') f90file//':'//subr// &
                  ': RFM0,M0,SUM(XP),POA=',RFM0,M0,sum(XP),POA
             stop
          end if
       end do

       M0 = MX
    end if !// if (LITERATE2) then

    !// --------------------------------------------------------------------
  end subroutine FINDM0
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine FM0(NHC,NOX,NP, KOM,XP,POA,M0,RFM0)
    !// --------------------------------------------------------------------
    !// Calculates lefthand side of Eq. 5 in Chung and Seinfeld (JGR, 2002)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NHC, NOX, NP
    real(r8), intent(in) :: KOM(NHC,NOX,NP), XP(NHC,NOX,NP)
    real(r8), intent(in) :: POA, M0
    !// Output
    real(r8), intent(out) :: RFM0

    !// Locals
    integer :: N1, N2, N3
    real(r8)  :: rsum
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'FM0'
    !// --------------------------------------------------------------------

    if (M0 .lt. 0._r8) then
       write(6,'(a,i7,a)') f90file//':'//subr// &
            ': M0 less than 0 ----STOP----'
       stop
    end if
    rsum = 0._r8
    do N3 = 1, NP
       do N2 = 1, NOX
          do N1 = 1, NHC
             rsum = rsum + XP(N1,N2,N3) / (1._r8 + KOM(N1,N2,N3) * M0)
          end do
       end do
    end do
    !// Return value
    RFM0 = RSUM
    !// Note that M0 may be zero at this point; need to check
    !// due to division by M0.
    if (M0 .gt. 0._r8) then
       RFM0 = RFM0 + POA/M0
    else if (M0 .eq. 0._r8) then
       !// Should return high value if POA>0 and M0=0
       if (POA .gt. 0._r8) RFM0 = RFM0 + 1000._r8
    end if

    !// --------------------------------------------------------------------
  end subroutine FM0
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine soa_diag_drydep(DRYDEP,DV1,I,J)
    !// --------------------------------------------------------------------
    !// Accumulate loss of condensed SOA through dry deposition.
    !// Wet deposition is diagnosed separately.
    !//
    !// DRYDEP is positive!
    !//
    !// Amund Sovde, July 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: TRACER_ID_MAX
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I,J
    real(r8), intent(in) :: DRYDEP(TRACER_ID_MAX), DV1
    !// Locals
    integer :: N, TRNR
    !// --------------------------------------------------------------------

    do N = 1, ndep_soa
       TRNR = trsp_idx(soa_deps(N))
       !// Diagnose is in units of grams [g]
       TOTAL_SOA_DDEP(I,J,N) = TOTAL_SOA_DDEP(I,J,N) &
            - DRYDEP(soa_deps(N)) & !// molecules/cm3
            * DV1 * 1.e6_r8 / AVOGNR * TMASS(TRNR) !// g
    end do

    !// --------------------------------------------------------------------
  end subroutine soa_diag_drydep
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine soa_diag_separate(PROD,I,J,LSTART,LEND)
    !// --------------------------------------------------------------------
    !// Accumulate condensation or evaporation of SOA.
    !// Positive PROD, i.e. condensation, is diagnosed as production,
    !// while negative PROD is diagnosed as evaporation.
    !//
    !// PROD is the difference (after-before) of the separation routine,
    !// so for each grid box it contains only the change in steady state.
    !// Positive production must then be regarded as condensation, and
    !// negative production as evaporation.
    !//
    !// Called from both tropchem and stratchem.
    !// Unit of PROD is mass [kg].
    !//
    !// Amund Sovde, July 2014, July 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: TRACER_ID_MAX
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I,J,LSTART,LEND
    real(r8), intent(in) :: PROD(TRACER_ID_MAX,LPAR) !// [kg/gridbox]
    !// Locals
    integer :: N,L, NSOA
    !// --------------------------------------------------------------------

    do N = 1, ndep_soa
       NSOA = soa_deps(N)
       do L = LSTART, LEND
          !// Diagnose is in units of grams [g]
          if (PROD(NSOA,L) .gt. 0._r8) then
             !// Diagnose positive production as production/condensation
             TOTAL_SOA_PROD(I,J,L,N) = TOTAL_SOA_PROD(I,J,L,N) &
                  + PROD(NSOA,L) * 1000._r8  !// kg -> g
          else
             !// Diagnose negative production as evaporation
             TOTAL_SOA_EVAP(I,J,L,N) = TOTAL_SOA_EVAP(I,J,L,N) &
                  + PROD(NSOA,L) * 1000._r8  !// kg -> g
          end if
       end do
    end do

    do N = 1, ngas_soa
       NSOA = soa_gas(N)
       do L = LSTART, LEND
          !// Diagnose is in units of grams [g]
          if (PROD(NSOA,L) .gt. 0._r8) then
             !// Production = evaproation of aerosols
             TOTAL_SOAG_SEPPROD(I,J,L,N) = TOTAL_SOAG_SEPPROD(I,J,L,N) &
                  + PROD(NSOA,L) * 1000._r8  !// kg -> g
          else
             !// Loss = condensation to aerosols
             TOTAL_SOAG_SEPCOND(I,J,L,N) = TOTAL_SOAG_SEPCOND(I,J,L,N) &
                  + PROD(NSOA,L) * 1000._r8  !// kg -> g
          end if
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine soa_diag_separate
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine soa_diag_soagas(PROD,I,J,LSTART,LEND)
    !// --------------------------------------------------------------------
    !// Accumulate chemically produced gasphase SOA species (SOAGASes).
    !//
    !// PROD is the difference (after-before) of chemistry.
    !//
    !// Called from both tropchem and stratchem.
    !// Unit of PROD is mass [kg].
    !//
    !// Amund Sovde, May 2016
    !// --------------------------------------------------------------------
    use cmn_size, only: TRACER_ID_MAX
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I,J,LSTART,LEND
    real(r8), intent(in) :: PROD(TRACER_ID_MAX,LPAR) !// [kg/gridbox]
    !// Locals
    integer :: N,L, NSOA
    !// --------------------------------------------------------------------

    do N = 1, ngas_soa
       NSOA = soa_gas(N)
       do L = LSTART, LEND
          !// Diagnose is in units of grams [g]
          if (PROD(NSOA,L) .gt. 0._r8) then
             !// Diagnose positive production as production/condensation
             TOTAL_SOAG_PROD(I,J,L,N) = TOTAL_SOAG_PROD(I,J,L,N) &
                  + PROD(NSOA,L) * 1000._r8  !// kg -> g
          end if
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine soa_diag_soagas
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine soa_diag2file(NDAY,NDAYI)
    !// --------------------------------------------------------------------
    !// Put out SOA diagnoses to file.
    !// Called from p-main.
    !//
    !// IMPORTANT:
    !//   Controlled by BUDGET calendar!
    !//   This is due to usage of STTTND and STTTN0.
    !//
    !//   Written to file is (in unit [g]):
    !//   - Accumulated production from NDAY0 to NDAY
    !//   - Accumulated evaporation from NDAY0 to NDAY
    !//   - Accumulated dry deposition from NDAY0 to NDAY
    !//   - Accumulated large scale wet scavenging from NDAY0 to NDAY
    !//   - Accumulated convective wet scavenging from NDAY0 to NDAY
    !//   - instantaneous burden
    !//
    !// Amund Sovde, June 2014, July 2012
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_ctm, only: STT, MPBLK, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, &
         JDATE, JMON, JYEAR
    use cmn_chem, only: TNAME
    use cmn_diag, only: NTND, NTND_CNSCAV, STTTND, STTTN0, TLDIAG, NDAY0
    use cmn_oslo, only: CONVWASHOUT, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NDAYI
    !// Locals
    integer :: TRNR,I,J,II,JJ,L,N,MP
    real(r8) :: tot_prod, tot_evap, tot_ddep, tot_lscav, tot_cscav, tot_burden
    real(r8) :: tot_SepProdGas, tot_SepCondGas
    real(r8) :: SCAV3D(IPAR,JPAR,LPAR,NDEP_SOA)
    !// For file stuff
    integer :: ifnr
    logical :: fnr_ok
    character(len=80) :: filename
    character(len=4)  :: cyear
    !// --------------------------------------------------------------------

    !// Open result file
    !// Find non-used file number for input file
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do

    !// Open file (it was initialised in subroutine soa_init)
    open(ifnr,file=filesoabudget,form='unformatted',position='append')
    write(ifnr) NDAY,NDAY0           !// This day and start of diagnose
    write(ifnr) real(TOTAL_SOA_PROD, r4) !// Production from SOA_separate
    write(ifnr) real(TOTAL_SOA_EVAP, r4) !// Evaporation from SOA_separate
    write(ifnr) real(TOTAL_SOA_DDEP, r4) !// Loss due to dry deposition

    tot_prod = sum(TOTAL_SOA_PROD)
    tot_evap = sum(TOTAL_SOA_EVAP)
    tot_ddep = sum(TOTAL_SOA_DDEP)

    !// Even though wet loss is diagnosed using the CTM3 core diagnostic
    !// arrays, namely the real*4 STTTND, we do a separate diagnotic here
    !// real*8. For convective scavenging, the CTM3 diagnose CONVWASHOUT is
    !// used. This means that SOA 3D diagnose is tied to BUDGET calendar
    !// in LxxTyy-file.

    !// Add large scale wet scav loss to 3D loss
    tot_lscav = sum(TOTAL_SOA_LSSCAV)
    write(ifnr) real(TOTAL_SOA_LSSCAV, r4) !// Large scale scavenging


    !// Add convective wet scav loss to 3D loss
    tot_cscav = 0._r8
    SCAV3D(:,:,:,:) = 0._r8
    if (NTND_CNSCAV .gt. 0) then
       do MP = 1, MPBLK
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do N = 1, ndep_soa
                   TRNR = trsp_idx(SOA_DEPS(N))
                   do L = 1, LPAR
                      SCAV3D(I,J,L,N) = &
                           + CONVWASHOUT(L,TRNR,II,JJ,MP) * 1000._r8 !// kg -> g
                   end do
                end do
             end do
          end do
       end do
       do N = 1, ndep_soa
          TRNR = trsp_idx(SOA_DEPS(N))
          tot_cscav = tot_cscav + STTTN0(TRNR,NTND_CNSCAV) * 1000._r8 !// kg -> g
       end do
    end if
    write(ifnr) real(SCAV3D, r4) !// Convective scavenging

    !// Burden [g] (instantaneous values)
    SCAV3D(:,:,:,:) = 0._r8
    do N = 1, ndep_soa
       TRNR = trsp_idx(SOA_DEPS(N))
       SCAV3D(:,:,:,N) = STT(:,:,:,TRNR) * 1000._r8
    end do
    write(ifnr) real(SCAV3D, r4) !// Burden (instantaneous)
    tot_burden = sum(SCAV3D)
    close(ifnr)


    !// Print out total production and loss
    write(6,'(a,5es11.3,f5.1)') '* SOA AVG P/E/D/L/C/life:', &
         tot_prod, tot_evap, tot_ddep, tot_lscav, tot_cscav, &
         -tot_burden / (tot_ddep + tot_lscav + tot_cscav) &
           * real(NDAY-NDAY0+1, r8)

    !// Diagnoses
    TOTAL_SOA_DDEP(:,:,:) = 0._r8
    TOTAL_SOA_PROD(:,:,:,:) = 0._r8
    TOTAL_SOA_EVAP(:,:,:,:) = 0._r8
    TOTAL_SOA_LSSCAV(:,:,:,:) = 0._r8

    !// Reset old values for nopsdiag
    oldProd(:) = 0._r8
    oldEvap(:) = 0._r8
    oldDDep(:) = 0._r8
    !// Must reset oldLscv and oldCscv also!
    !// After this routine is completed, p-main will call TBGT_G to
    !// reset STTTND and STTTN0.
    oldLscv(:) = 0._r8
    oldCscv(:) = 0._r8

    !// Write daily accumulated data to file
    !// This routine is called after JDAY/JMON/JYEAR have been updated
    if (JDATE.eq.1 .and. JMON.eq.1) then
       write(cyear,'(i4.4)') JYEAR-1 !// Has been updated to next year
    else
       write(cyear,'(i4.4)') JYEAR   !// Current year
    end if
    filename = 'soa_dailybudgets_'//cyear//'.dta'
    open(ifnr,file=filename,form='unformatted')
    write(ifnr) NDEP_SOA, 366
    write(ifnr) dailyBrdn !// Daily average
    write(ifnr) dailyProd !// Accumulated g/day
    write(ifnr) dailyEvap !// Accumulated g/day
    write(ifnr) dailyDDep !// Accumulated g/day
    write(ifnr) dailyLscv !// Accumulated g/day
    write(ifnr) dailyCscv !// Accumulated g/day
    close(ifnr)
    !// Reset if new year
    !// Note that this will fail if diagnostics are not done 1Jan.
    if (JDATE.eq.1 .and. JMON.eq.1) then
       dailyProd(:,:) = 0._r8
       dailyEvap(:,:) = 0._r8
       dailyDDep(:,:) = 0._r8
       dailyLscv(:,:) = 0._r8
       dailyCscv(:,:) = 0._r8
       dailyBrdn(:,:) = 0._r8
    end if


    !// SOA GASES
    !// --------------------------------------------------------------------
    !// Open file (it was initialised in subroutine soa_init)
    open(ifnr,file=filesoagbudget,form='unformatted',position='append')
    write(ifnr) NDAY,NDAY0           !// This day and start of diagnose
    write(ifnr) real(TOTAL_SOAG_PROD, r4) !// Produced in chemistry
    write(ifnr) real(TOTAL_SOAG_SEPPROD, r4) !// Produced in SOA_separate
    write(ifnr) real(TOTAL_SOAG_SEPCOND, r4) !// Condensed in SOA_separate

    tot_prod = sum(TOTAL_SOAG_PROD)
    tot_SepProdGas = sum(TOTAL_SOAG_SEPPROD)
    tot_SepCondGas = sum(TOTAL_SOAG_SEPCOND)

    !// Even though wet loss is diagnosed using the CTM3 core diagnostic
    !// arrays, namely the real*4 STTTND, we do a separate diagnotic here
    !// real*8. For convective scavenging, the CTM3 diagnose CONVWASHOUT is
    !// used. This means that SOA 3D diagnose is tied to BUDGET calendar
    !// in LxxTyy-file.

    !// Add large scale wet scav loss to 3D loss
    tot_lscav = sum(TOTAL_SOAG_LSSCAV)
    write(ifnr) real(TOTAL_SOAG_LSSCAV, r4) !// Large scale scavenging


    !// Add convective wet scav loss to 3D loss
    tot_cscav = 0._r8
    SCAV3D(:,:,:,:) = 0._r8
    if (NTND_CNSCAV .gt. 0) then
       do MP = 1, MPBLK
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do N = 1, ngas_soa
                   TRNR = trsp_idx(SOA_GAS(N))
                   do L = 1, LPAR
                      SCAV3D(I,J,L,N) = &
                           + CONVWASHOUT(L,TRNR,II,JJ,MP) * 1000._r8 !// kg -> g
                   end do
                end do
             end do
          end do
       end do
       do N = 1, ngas_soa
          TRNR = trsp_idx(SOA_GAS(N))
          tot_cscav = tot_cscav + STTTN0(TRNR,NTND_CNSCAV) * 1000._r8 !// kg -> g
       end do
    end if
    write(ifnr) real(SCAV3D, r4) !// Convective scavenging

    !// Burden [g] (instantaneous values)
    SCAV3D(:,:,:,:) = 0._r8
    do N = 1, ngas_soa
       TRNR = trsp_idx(SOA_GAS(N))
       SCAV3D(:,:,:,N) = STT(:,:,:,TRNR) * 1000._r8
    end do
    write(ifnr) real(SCAV3D, r4) !// Burden (instantaneous)
    tot_burden = sum(SCAV3D)
    close(ifnr)


    !// Print out total production and loss
    write(6,'(a,5es11.3,f5.1)') '* SOAG AVG P/SP/SC/L/C/life:', &
         tot_prod, tot_sepprodgas, tot_lscav, tot_cscav, &
         -tot_burden / (tot_sepcondGas + tot_lscav + tot_cscav) &
           * real(NDAY-NDAY0+1, r8)

    !// Diagnoses
    TOTAL_SOAG_PROD(:,:,:,:) = 0._r8
    TOTAL_SOAG_SEPPROD(:,:,:,:) = 0._r8
    TOTAL_SOAG_SEPCOND(:,:,:,:) = 0._r8
    TOTAL_SOA_LSSCAV(:,:,:,:) = 0._r8

    !// Reset old values for nopsdiag
    oldProdGas(:) = 0._r8
    oldSepProdGas(:) = 0._r8
    oldSepCondGas(:) = 0._r8
    !// Must reset oldLscv and oldCscv also!
    !// After this routine is completed, p-main will call TBGT_G to
    !// reset STTTND and STTTN0.
    oldLscvGas(:) = 0._r8
    oldCscvGas(:) = 0._r8

    !// Write daily accumulated data to file
    !// This routine is called after JDAY/JMON/JYEAR have been updated
    if (JDATE.eq.1 .and. JMON.eq.1) then
       write(cyear,'(i4.4)') JYEAR-1 !// Has been updated to next year
    else
       write(cyear,'(i4.4)') JYEAR   !// Current year
    end if
    filename = 'soag_dailybudgets_'//cyear//'.dta'
    open(ifnr,file=filename,form='unformatted')
    write(ifnr) NGAS_SOA, 366
    write(ifnr) dailyBrdnGas !// Daily average
    write(ifnr) dailyProdGas !// Accumulated g/day
    write(ifnr) dailySepProdGas !// Accumulated g/day
    write(ifnr) dailySepCondGas !// Accumulated g/day
    write(ifnr) dailyLscvGas !// Accumulated g/day
    write(ifnr) dailyCscvGas !// Accumulated g/day
    close(ifnr)
    !// Reset if new year
    !// Note that this will fail if diagnostics are not done 1Jan.
    if (JDATE.eq.1 .and. JMON.eq.1) then
       dailyProdGas(:,:) = 0._r8
       dailySepProdGas(:,:) = 0._r8
       dailySepCondGas(:,:) = 0._r8
       dailyLscvGas(:,:) = 0._r8
       dailyCscvGas(:,:) = 0._r8
       dailyBrdnGas(:,:) = 0._r8
    end if


    !// --------------------------------------------------------------------
  end subroutine soa_diag2file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine soa_nopsdiag(NDAY,NMET,NOPS,JDAY,DTOPS)
    !// --------------------------------------------------------------------
    !// Print diags of SOA each NOPS.
    !// Note that this diagnostic uses the same diagnostic arrays containing
    !// accumulated values, i.e. the same arrays as in soa_diag2file.
    !//
    !// IMPORTANT:
    !//   If you do changes here, make sure you know what you are doing!
    !//   This routine must not be called after soa_diag2file.
    !//
    !// Amund Sovde, June-July 2014
    !// --------------------------------------------------------------------
    use cmn_ctm, only: STT
    use cmn_diag, only: NTND, STTTN0, TLDIAG,NTND_CNSCAV,NTND_LSSCAV
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NMET, NOPS, JDAY
    real(r8), intent(in)  :: DTOPS
    !// Locals
    integer :: NSCAV, TRNR,I,J,II,JJ,L,N,MP
    real(r8) :: tot_lscav, tot_cscav, tot_prod, tot_evap, tot_ddep, &
              tot_brd, ZDT, R1
    real(r8) :: tot_lscavGas, tot_cscavgas, tot_prodGas, tot_sepCondGas, &
         tot_sepProdGas, tot_brdGas
    real(r8), dimension(NDEP_SOA) :: &
         newProd,newEvap,newDDep, newLscv,newCscv, burden
    real(r8), dimension(NGAS_SOA) :: &
         newProdGas,newLscvGas,newCscvGas, newSepProdGas, newSepCondGas, &
         burdenGas
    real(r8) :: tmplifetime(8)
    logical, parameter :: LPrntBud = .false.
    logical, parameter :: LPrntLif = .false.
    !// --------------------------------------------------------------------

    !// Time step is NOPS step
    ZDT = 1._r8/DTOPS

    !// Initialise
    newProd(:) = 0._r8
    newEvap(:) = 0._r8
    newDDep(:) = 0._r8
    newLscv(:) = 0._r8
    newCscv(:) = 0._r8

    newProdGas(:) = 0._r8
    newLscvGas(:) = 0._r8
    newCscvGas(:) = 0._r8
    newSepProdGas(:) = 0._r8
    newSepCondGas(:) = 0._r8

    !// All global values used here (also STTTN0) are real8, as they
    !// should be.

    !// Find burden of LS scav since last NOPS
    if (NTND_LSSCAV .gt. 0) then
       !// SOA AEROSOLS
       do N = 1, ndep_soa
          TRNR = trsp_idx(SOA_DEPS(N))
          !// LS scav accumulated up to end of this time step (convert kg -> g)
          newLscv(N) = STTTN0(TRNR,NTND_LSSCAV) * 1000._r8
       end do
       !// What was scavenged by LS since during the last NOPS
       tot_lscav = sum(newLscv(:) - oldLscv(:))

       !// SOA GASES
       do N = 1, ngas_soa
          TRNR = trsp_idx(SOA_GAS(N))
          !// LS scav accumulated up to end of this time step (convert kg -> g)
          newLscvGas(N) = STTTN0(TRNR,NTND_LSSCAV) * 1000._r8
       end do
       !// What was scavenged by LS since during the last NOPS
       tot_lscavGas = sum(newLscvGas(:) - oldLscvGas(:))

    end if


    !// Find burden of CNV scav since last NOPS
    tot_cscav = 0._r8
    if (NTND_CNSCAV .gt. 0) then
       !// SOA AEROSOLS
       do N = 1, ndep_soa
          TRNR = trsp_idx(SOA_DEPS(N))
          !// CNV scav accumulated up to end of this time step (convert kg -> g)
          newCscv(N) = STTTN0(TRNR,NTND_CNSCAV) * 1000._r8
       end do
       !// What was scavenged by CNV since during the last NOPS
       tot_cscav = sum(newCscv(:) - oldCscv(:))

       !// SOA GASES
       do N = 1, ngas_soa
          TRNR = trsp_idx(SOA_GAS(N))
          !// CNV scav accumulated up to end of this time step (convert kg -> g)
          newCscvGas(N) = STTTN0(TRNR,NTND_CNSCAV) * 1000._r8
       end do
       !// What was scavenged by CNV since during the last NOPS
       tot_cscavGas = sum(newCscvGas(:) - oldCscvGas(:))
    end if

    !// How much has been produced up to end of this time step [g]
    do N = 1, ndep_soa
       newProd(N) = sum(TOTAL_SOA_PROD(:,:,:,N))
    end do
    !// What was accumulated since last time
    tot_prod = sum(newProd(:) - oldProd(:))

    !// How much has been evaporated up to end of this time step [g]
    do N = 1, ndep_soa
       newEvap(N) = sum(TOTAL_SOA_EVAP(:,:,:,N))
    end do
    !// What was accumulated since last time
    tot_evap = sum(newEvap(:) - oldEvap(:))

    !// How much has been lost to drydep up to end of this time step [g]
    do N = 1, ndep_soa
       newDDep(N) = sum(TOTAL_SOA_DDEP(:,:,N))
    end do
    !// What was accumulated since last time
    tot_ddep = sum(newDDep(:) - oldDDep(:))

    !// Sum up burden
    do N = 1, ndep_soa
       TRNR = trsp_idx(SOA_DEPS(N))
       burden(N) = sum(STT(:,:,:,TRNR)) * 1000._r8 !// [g]
    end do
    tot_brd = sum(burden)


    !// How much has been produced chemically up to end of this time step [g]
    do N = 1, ngas_soa
       newProdGas(N) = sum(TOTAL_SOAG_PROD(:,:,:,N))
    end do
    !// What was accumulated since last time
    tot_prodGas = sum(newProdGas(:) - oldProdGas(:))

    !// How much soagas has come from evaproating SOA?
    do N = 1, ngas_soa
       newSepProdGas(N) = sum(TOTAL_SOAG_SEPPROD(:,:,:,N))
    end do
    !// What was accumulated since last time
    tot_sepProdGas = sum(newSepProdGas(:) - oldSepProdGas(:))

    !// How much soagas has been lost by condensating to SOA?
    do N = 1, ngas_soa
       newSepCondGas(N) = sum(TOTAL_SOAG_SEPCOND(:,:,:,N))
    end do
    !// What was accumulated since last time
    tot_sepCondGas = sum(newSepCondGas(:) - oldSepCondGas(:))

    !// Sum up burden SOAGAS
    do N = 1, ngas_soa
       TRNR = trsp_idx(SOA_GAS(N))
       burdenGas(N) = sum(STT(:,:,:,TRNR)) * 1000._r8 !// [g]
    end do
    tot_brdGas = sum(burdenGas)

    !// SOA AEROSOLS
    !// Print burden as grams.
    !// Print tendencies/diagnostics as gram per seconds.
    !// Annual mean production was typically 55-69Tg/yr in CTM2. To convert
    !// from g/s to Tg/yr, multiply by 3.1536d-5 (g/s*1d-12*3600*24*365).
    !//    P = 9.52d5 g/s = 30 Tg/yr
    !//    P = 1.75d6 g/s = 55 Tg/yr
    if (LPrntBud) &
         write(6,'(a,7es10.2)') 'SOA B/P/E/D/L/C/N:', &
         tot_brd, & !// g
         tot_prod * ZDT, tot_evap * ZDT, tot_ddep * ZDT, & !// g/s
         tot_lscav * ZDT, tot_cscav * ZDT, &               !// g/s
         (tot_prod + tot_evap + tot_ddep + tot_lscav + tot_cscav) * ZDT

    !// Print lifetimes [days]
    R1 = ZDT * 86400._r8 !// [1/s * 86400s/day]
    !// The lifetimes are calculated from the loss rates of each process,
    !// and the total lifetime from all loss rates.
    !// Assuming steady state, production rate = loss rate, so
    !// we can also find lifetime from burden and production rate.
    !//
    !// In addition a few other "lifetimes" are calculated
    tmplifetime(1) = -tot_brd/(tot_evap * R1)  !// evaporation
    tmplifetime(2) = -tot_brd/(tot_ddep * R1)  !// dry deposition
    tmplifetime(3) = -tot_brd/(tot_lscav * R1) !// large scale scavenging
    tmplifetime(4) = -tot_brd/(tot_cscav * R1) !// convective scavenging
    !// Total lifetime without evaporation (for curiosity only)
    tmplifetime(5) = -tot_brd/((tot_ddep + tot_lscav + tot_cscav) * R1)
    !// Total lifetime prod+evap (for curiosity only)
    tmplifetime(6) = tot_brd/((tot_prod + tot_evap) * R1)
    !// Total lifetime from all loss rates
    tmplifetime(7) = -tot_brd/((tot_ddep+tot_lscav+tot_cscav+tot_evap)*R1)
    !// Total lifetime from production rate
    tmplifetime(8) = tot_brd/(tot_prod * R1)
    !/ Check for values not fitting output format (f6.1,7f6.1)
    do n=1,8
       if (tmplifetime(n) .gt. 999._r8) tmplifetime(n) = 999._r8
       if (tmplifetime(n) .lt. -99._r8) tmplifetime(n) = -99._r8
    end do

    !// Print Burden[Gg] and lifetimes [days]
    if (LPrntLif) &
         write(6,'(a,8f6.1)') 'SOA E/D/L/C/TL-E/P+E/TL/P:', &
         tmplifetime(1:8)   !// The lifetimes calculated above


    !// SOA GASES
    if (LPrntBud) &
         write(6,'(a,7es10.2)') 'SOAG B/P/L/C/SP/SC/N:', &
         tot_brdGas, & !// g
         tot_prodGas * ZDT, tot_lscavGas * ZDT, tot_cscavGas * ZDT, & !// g/s
         tot_sepProdGas * ZDT, tot_sepCondGas * ZDT, &                !// g/s
         (tot_prodGas + tot_lscavGas + tot_cscavGas &
          + tot_sepProdGas + tot_sepCondGas) * ZDT                    !// g/s

    !// In addition a few other "lifetimes" are calculated
    tmplifetime(1) = tot_brdGas/(tot_prodGas * R1) !// Chem. prod.
    tmplifetime(2) = -tot_brdGas/(tot_lscavGas * R1) !// large scale scavenging
    tmplifetime(3) = -tot_brdGas/(tot_cscavGas * R1) !// convective scavenging
    tmplifetime(4) = tot_brdGas/(tot_sepProdGas * R1) !// separation
    tmplifetime(5) = -tot_brdGas/(tot_sepCondGas * R1) !// separation
    !// Total lifetime from all loss rates
    tmplifetime(6) = -tot_brdGas / &
         ((tot_lscavGas + tot_cscavGas + tot_sepCondGas)*R1)
    !// Total lifetime from production rate
    tmplifetime(7) = tot_brdGas/((tot_prodGas + tot_sepProdGas) * R1)
    tmplifetime(8) = 0._r8
    !/ Check for values not fitting output format (f7.1,7f6.1)
    do n=1,8
       if (tmplifetime(n) .gt. 999._r8) tmplifetime(n) = 999._r8
       if (tmplifetime(n) .lt. -99._r8) tmplifetime(n) = -99._r8
    end do

    !// Print lifetimes [days]
    if (LPrntLif) &
         write(6,'(a,8f6.1)') 'SOAG P/L/C/SP/SC/TL/TP:', &
         tmplifetime(1:7)   !// The lifetimes calculated above


    !// Accumulate daily diagnostics
    !// These will be written to file in soa_diag2file
    if (NMET.eq.1 .and. NOPS.eq.1) then
       !// Need to set to zero at the beginning of day, otherwise we
       !// there may be numbers present when starting on a new year.
       dailyProd(:,JDAY) = 0._r8
       dailyEvap(:,JDAY) = 0._r8
       dailyDDep(:,JDAY) = 0._r8
       dailyLscv(:,JDAY) = 0._r8
       dailyCscv(:,JDAY) = 0._r8
       dailyBrdn(:,JDAY) = 0._r8

       dailyProdGas(:,JDAY) = 0._r8
       dailyLscvGas(:,JDAY) = 0._r8
       dailyCscvGas(:,JDAY) = 0._r8
       dailySepProdGas(:,JDAY) = 0._r8
       dailySepCondGas(:,JDAY) = 0._r8
       dailyBrdnGas(:,JDAY) = 0._r8

    end if
    !// Then accumulate
    dailyProd(:,JDAY) = dailyProd(:,JDAY) + (newProd(:) - oldProd(:))
    dailyEvap(:,JDAY) = dailyEvap(:,JDAY) + (newEvap(:) - oldEvap(:))
    dailyDDep(:,JDAY) = dailyDDep(:,JDAY) + (newDDep(:) - oldDDep(:))
    dailyLscv(:,JDAY) = dailyLscv(:,JDAY) + (newLscv(:) - oldLscv(:))
    dailyCscv(:,JDAY) = dailyCscv(:,JDAY) + (newCscv(:) - oldCscv(:))
    dailyBrdn(:,JDAY) = dailyBrdn(:,JDAY) + burden(:) / R1 !// Make average

    dailyProdGas(:,JDAY) = dailyProdGas(:,JDAY) &
         + (newProdGas(:) - oldProdGas(:))
    dailyLscvGas(:,JDAY) = dailyLscvGas(:,JDAY) &
         + (newLscvGas(:) - oldLscvGas(:))
    dailyCscvGas(:,JDAY) = dailyCscvGas(:,JDAY) &
         + (newCscvGas(:) - oldCscvGas(:))
    dailySepProdGas(:,JDAY) = dailySepProdGas(:,JDAY) &
         + (newSepProdGas(:) - oldSepProdGas(:))
    dailySepCondGas(:,JDAY) = dailySepCondGas(:,JDAY) &
         + (newSepCondGas(:) - oldSepCondGas(:))
    dailyBrdnGas(:,JDAY) = dailyBrdnGas(:,JDAY) &
         + burdenGas(:) / R1 !// Make average

    !// Update old values for diagnose next time
    oldProd(:) = newProd(:)
    oldEvap(:) = newEvap(:)
    oldDDep(:) = newDDep(:)
    oldLscv(:) = newLscv(:)
    oldCscv(:) = newCscv(:)

    oldProdGas(:) = newProdGas(:)
    oldLscvGas(:) = newLscvGas(:)
    oldCscvGas(:) = newCscvGas(:)
    oldSepProdGas(:) = newSepProdGas(:)
    oldSepCondGas(:) = newSepCondGas(:)
    !// --------------------------------------------------------------------
  end subroutine soa_nopsdiag
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine soa_diag_lsscav(BTT,BTTBCK,MP)
    !// --------------------------------------------------------------------
    !// Accumulate SOA removed by large scale scavenging.
    !//
    !// Called from subroutine wdep_diag_ls in diagnostics_scavenging.f90.
    !// Unit of PROD is mass [kg].
    !//
    !// Amund Sovde, August 2014
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    !// In
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in):: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in):: BTTBCK

    !// Locals
    integer :: N,TRNR, I,J,II,JJ, L
    double precision :: RTMP
    !// --------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        do N = 1, ndep_soa
          TRNR = trsp_idx(soa_deps(N))
          do L = 1, LPAR
             RTMP = BTT(L,TRNR,II,JJ) - BTTBCK(L,TRNR,II,JJ)
             !// Diagnose is in units of grams [g]
             TOTAL_SOA_LSSCAV(I,J,L,N) = TOTAL_SOA_LSSCAV(I,J,L,N) &
                  + RTMP * 1.e3_r8  !// kg -> g
          end do
        end do

        do N = 1, ngas_soa
          TRNR = trsp_idx(soa_gas(N))
          do L = 1, LPAR
             RTMP = BTT(L,TRNR,II,JJ) - BTTBCK(L,TRNR,II,JJ)
             !// Diagnose is in units of grams [g]
             TOTAL_SOAG_LSSCAV(I,J,L,N) = TOTAL_SOAG_LSSCAV(I,J,L,N) &
                  + RTMP * 1.e3_r8  !// kg -> g
          end do
        end do

      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine soa_diag_lsscav
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module soa_oslo
!//=========================================================================
