c-----------------------------------------------------------------------
c---(p-series.f)-----UCIrvine CTM  p-code 5.5 (9/07)
c---subs:   TSER_0, TSER_1, TSER_2
c
c      subroutine TSER_0 (KUNF)
c      subroutine TSER_1 (KSER)
c      subroutine TSER_2 (KUNF,KSER)
c
c-----------------------------------------------------------------------
      subroutine TSER_0 (KUNF)
c-----------------------------------------------------------------------
c---initialize the time series unformatted file with 
c---     a list of stations, tracers, number of points per day
c-----------------------------------------------------------------------
      use cmn_size, only: LPAR, NPAR
      use cmn_ctm, only: NTM, NROPSM, NRMETD
      use cmn_chem, only: TNAME
      use cmn_diag, only: NBOXS, ISTA, JSTA, TSTAX, LBGTS
      implicit none
c-----------------------------------------------------------------------
      integer, intent(in) :: KUNF   !fortran unit for unformatted series
      integer N, NTDO, KTDO, NTRA(NPAR)
      character(len=8)  TNTRA(NPAR)

      write(KUNF) NBOXS,(ISTA(N),JSTA(N),N=1,NBOXS),(TSTAX(N),N=1,NBOXS)

         NTDO = 0
       do N = 1,NTM
        if (LBGTS(N)) then
         NTDO = NTDO + 1
         NTRA(NTDO) = N
c---         TNTRA(NTDO) = TNAME(N)
          write(TNTRA(NTDO),'(a8)') TNAME(N)
        endif
       enddo

      write(KUNF) NTDO, (NTRA(N), N=1,NTDO),(TNTRA(N),N=1,NTDO)
        
        KTDO = NROPSM*NRMETD

      write(KUNF) KTDO, LPAR

      return
      end
      

      
c-----------------------------------------------------------------------
      subroutine TSER_1 (KSER)
c-----------------------------------------------------------------------
c---stores DAILY time series of mass fraction for ALL tracers & L-levels 
c---   at specific stations (NSBOX)
c-----------------------------------------------------------------------
      use cmn_precision, only: r4
      use cmn_size, only: LPAR
      use cmn_ctm, only: STT, AIR, NTM
      use cmn_diag, only: NBOXS, NSTPAR, STTS4, ISTA, JSTA
      implicit none
c-----------------------------------------------------------------------
      integer, intent(inout) :: KSER
      integer II,JJ,NS,N,L

      KSER = KSER+1
        KSER = min(KSER, NSTPAR)
      do NS = 1,NBOXS
        II = ISTA(NS)
        JJ = JSTA(NS)
        do N = 1,NTM
         do L = 1,LPAR
           STTS4(L,KSER,NS,N) = real(STT(II,JJ,L,N)/AIR(II,JJ,L), r4)
         enddo
        enddo
      enddo

      return
      end



c-----------------------------------------------------------------------
      subroutine TSER_2 (KUNF,KSER)
c-----------------------------------------------------------------------
c---write out unformatted times series:
c---   only for specified days (JDO_S calendar) and tracers (LBGTS)
c---   KUNF = fortran unit for unformatted series
c---   KSER = no. of daily terms = NROPSM*NRMETD (e.g., 8, 24, 48=dim)
c
c-----------------------------------------------------------------------
      use cmn_size, only: LPAR
      use cmn_ctm, only: NTM, IYEAR, IDAY, JYEAR, JDAY
      use cmn_diag, only: NBOXS, LBGTS, STTS4
      implicit none
c-----------------------------------------------------------------------
      integer, intent(in) :: KUNF
      integer, intent(in) :: KSER
      integer N,K,L,NS

c---Note that for each tracer (1:NTM) dumped (from LBGTS(1:NTM)) there is
c---     a header line with year, day, and dimensions of the unf array
c---       IYEAR & IDAY are from start of run (IDAY can be > 366)
c---     a real*4 array dimension (KSER(eg,1:24 hr), LM, NBOXS(stations))

      do N=1,NTM
       if (LBGTS(N)) then
        write(KUNF) IYEAR,IDAY,JYEAR,JDAY,N, LPAR,KSER,NBOXS
        write(KUNF) (((STTS4(L,K,NS,N),K=1,KSER),L=1,LPAR),NS=1,NBOXS)
       endif
      enddo

      return
      end
