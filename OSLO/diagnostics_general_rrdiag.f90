!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Amund Sovde Haslerud, January 2017
!//=========================================================================
!// CTM diagnostics
!//=========================================================================
module diagnostics_general_rrdiag
  !//-----------------------------------------------------------------------
  !// MODULE: diagnostics_general_rrdiag
  !// DESCRIPTION: Routines for diagnostics on aerosol uptake.
  !// ----------------------------------------------------------------------
  !// Diagnostics for the Oslo CTM3.
  !// Contains:
  !//   subroutine diag_ohch4n2o_init
  !//   subroutine init_lifetime
  !//   subroutine du_columns
  !//   subroutine write_ducolumns
  !//   subroutine init_daily_diag
  !//   subroutine daily_diag_output
  !//   subroutine nops_diag
  !//   subroutine mp_diag
  !//   subroutine TBGT_2FILE
  !//   subroutine REPORTS_CHEMISTRY
  !//   subroutine sumup_burden_and_lifetime
  !//   subroutine report_burden_and_lifetime
  !//   subroutine report_ch4n2o
  !//   subroutine diag_burden_snapshot
  !//   subroutine write_snapshot
  !//   subroutine write_snapshot_the
  !//   subroutine ch4n2o_burden
  !//   subroutine ch4_loss3
  !//   subroutine n2o_loss3
  !//   subroutine diags_tpset
  !//   subroutine report_negO3
  !//   subroutine tnd_emis_daily
  !//   subroutine tnd_emis2file
  !//   subroutine save_chemPL
  !//   subroutine chembud_output
  !//
  !//
  !// RR diagnostics was meant for evaluating the aerosol surface area /
  !// uptake coefficients. To produce output, locate the
  !//   if (LTBGT .or. LEND) then
  !// in pmain, just before end of NDAY-loop and insert
  !//   call output_rr(JYEAR,JMON,JDATE,NDAY)
  !//
  !// Amund Sovde Haslerud, January 2017
  !// Ole Amund Sovde, February 2015, August 2013, February 2012,
  !//              October 2008
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, MPBLK, IDBLK, JDBLK
  use cmn_parameters, only: M_AIR
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Diagnose aerosol uptake rates
  real(r8), dimension(LPAR,IPAR,JPAR) :: RRN2O5DIAG, RRHO2DIAG, RRQAERDIAG
  real(r8) :: RRcount
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file='diagnostics_general_rrdiag.f90'
  !// ----------------------------------------------------------------------
  save !// All variables are to be saved.
  private
  public diagnose_rr_init, diagnose_rr, output_rr
  !// ----------------------------------------------------------------------

contains


  !// ----------------------------------------------------------------------
  subroutine diagnose_rr_init()
    !// --------------------------------------------------------------------
    !// Initialize aerosol uptake reaction diagnostics.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------

    RRN2O5DIAG(:,:,:) = 0._r8
    RRHO2DIAG(:,:,:)  = 0._r8
    RRQAERDIAG(:,:,:) = 0._r8
    !RRNO2DIAG(:,:,:)  = 0._r8

    RRcount = 0._r8

    !// --------------------------------------------------------------------
  end subroutine diagnose_rr_init
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  subroutine diagnose_rr(I,J, MP, RR_N2O5_H2O_AER, RR_HO2_AER, &
       RR_QAER, RR_NO2_SOOT)
    !// --------------------------------------------------------------------
    !// Add up aerosol uptake rates.
    !// Only done in troposphere.
    !// If stratosphere to be included, RRcount must be updated with care.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: I, J, MP
    real(r8), dimension(LPAR), intent(in) :: &
         RR_N2O5_H2O_AER, RR_HO2_AER, RR_QAER, RR_NO2_SOOT 
    !// --------------------------------------------------------------------

    RRN2O5DIAG(:,I,J) = RRN2O5DIAG(:,I,J) + RR_N2O5_H2O_AER(:)
    RRHO2DIAG(:,I,J)  = RRHO2DIAG(:,I,J) + RR_HO2_AER(:)
    RRQAERDIAG(:,I,J) = RRQAERDIAG(:,I,J) + RR_QAER(:)
    !RRNO2DIAG(:,I,J)  = RRNO2DIAG(:,I,J) + RR_NO2_SOOT(:)

    if (MP .eq. 1 .and. J.eq.1 .and. I.eq.1) RRcount = RRcount + 1._r8
    !// --------------------------------------------------------------------
  end subroutine diagnose_rr
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine output_rr(YEAR, MONTH, DATE, NDAY)
    !// --------------------------------------------------------------------
    !// Put diagnosed (average) RR to file.
    !// Follows budget calendar.
    !// Called from pmain.
    !// --------------------------------------------------------------------
    use cmn_ctm, only: AREAXY
    use cmn_diag, only: NDAY0
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: YEAR, MONTH, DATE, NDAY
    !// Locals
    integer :: ifnr
    character(len=80) :: filename
    character(len=8) :: datestamp
    real(r8) :: zcount
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='output_rr'
    !// --------------------------------------------------------------------

    !// Make average
    zcount = 1._r8 / RRcount

    RRN2O5DIAG(:,:,:) = RRN2O5DIAG(:,:,:) * zcount
    RRHO2DIAG(:,:,:)  = RRHO2DIAG(:,:,:) * zcount
    RRQAERDIAG(:,:,:) = RRQAERDIAG(:,:,:) * zcount
    !RRNO2DIAG(:,:,:)  = RRNO2DIAG(:,:,:) * zcount

    !// Find file number
    ifnr = get_free_fileid()

    !// YYYYMMDD for filenames, taken from JYEAR,JMON,JDATE, are
    !// the arguments of this subroutine.
    !// They are already updated by calendar.
    write(datestamp(1:8),'(i4.4,2i2.2)') YEAR,MONTH,DATE
    filename = 'rr_avg'//datestamp//'.dta'

    !// Write file
    open(ifnr,file=filename,form='unformatted')
    write(ifnr) IPAR,JPAR,LPAR,NDAY0,NDAY
    write(ifnr) AREAXY

    write(ifnr) 3 !// Increase to 4 if RRNO2DIAG is included
    write(ifnr) RRN2O5DIAG
    write(ifnr) RRHO2DIAG
    write(ifnr) RRQAERDIAG
    !write(ifnr) RRNO2DIAG
    close(ifnr)

    write(6,'(a)') f90file//':'//subr// &
            ': Put RR-averages to file '//trim(filename)

    !// --------------------------------------------------------------------
  end subroutine output_rr
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module diagnostics_general_rrdiag
!//=========================================================================
