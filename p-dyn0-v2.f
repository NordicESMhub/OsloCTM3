c-----------------------------------------------------------------------
c---(pdyn0.f)      p-code 6.0a  Alternate 2 (Prather 5/2011)
c---Advective transport, sets up met fields, pressure filters, CFL limits
c---subroutines:  DYN0, PFILTR, EPZ_TQ, EPZ_UV, EPZ_P, CFLADV
c
c-----------------------------------------------------------------------
      subroutine DYN0(DTWIND)
C-----------------------------------------------------------------------
C sets up the advective fields, called from MAIN
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, LPAR
      use cmn_ctm, only: ALFA, BETA, GAMA, AIR, AIRX,
     &     XYZA, XYZB, XYA, XYB, IMEPZ, AREAXY
      use cmn_met, only: U, V, Q, T, P
      use cmn_parameters, only: G0
      implicit none
C-----------------------------------------------------------------------

      real(r8), intent(in) ::  DTWIND


      real(r8), dimension(IPAR,JPAR)      ::  PCTM, MERR,SUMAD, SUMAQ
      real(r8), dimension(IPAR+1,JPAR)    ::  AX
      real(r8), dimension(IPAR,JPAR+1)    ::  BX
      real(r8), dimension(IPAR,JPAR,LPAR) ::  AIRNEW, AIRQKG
      real(r8)  RMSERR1,RMSERR2, ZDTW,G100,AIRQAV, DEL0,DELI,DELB
     &       ,DALFAI(IPAR)

      integer :: IIX,JJX,NITR,IMH,I,J,L,II,III
      logical :: LSP,LNP,LEW
      integer :: IM, JM, LM

C-----------------------------------------------------------------------
C----------------- UNITS OF AIR MASS AND TRACER = (kg)  ----------------
C----Air mass (kg) is given by area (m^2) * pressure thickness (Pa) / g0
C---- AREAXY(I,J)= area of [I,J] (m^2)
C---- P(I,J)     = surf pressure (Pa) averaged in extended zone.
C---- Q(I,J,L)   = specific humidity of grid box (kg H2O / kg wet air)
C                   averaged in extended zone.
C---- AIR (I,J,L)= dry-air mass (kg) in each box as calculated in CTM
C----              at the beginning of each time step.
C---- PCTM(I,J)  = inferred wet-air (total) surf press (Pa) calc. in CTM
C----              (using SUMAQ & AIR-X-NEW)
C---- DTWIND     = time step (sec) that applies to the averaged wind fields
C----              i.e., the time between successive pressures.

C----LOCAL:
C---- AIRNEW(I,J,L)= new dry-air mass in each CTM box after horizontal
C----                divergence (ALFA+BETA) over time step DTWIND (sec)
C---- AIRX(I,J,L)= expected dry-air mass in each CTM box after calculating the
C----              vertical divergence (GAMA)   (also used for GCM dry mass)
C----            = XYZA(I,J,L) + XYZB(I,J,L)*PCTM(I,J) - AIRQKG(I,J,L)
C----Local:
C---- AIRQKG(I,J,L)= kg of H2O in each grid box from GCM P's & q's
C---- SUMAD(I,J) = column of dry air (kg)
C---- SUMAQ(I,J) = column of water (kg)
C----
C----Assume that we have "wet-air" mass fluxes across each boundary
C----   U(I,J,L)   ==> [I,J,L] ==>  U(I+1,J,L)   (kg/s)
C----   V(I,J,L)   ==> [I,J,L] ==>  V(I,J+1,L)   (kg/s)
C----
C----Convert to "dry-air" mass flux in/out of box using average Q at boundary
C----   ALFA(I,J,L)   ==> [I,J,L] ==>  ALFA(I+1,J,L)   (kg/s)
C----   BETA(I,J,L)   ==> [I,J,L] ==>  BETA(I,J+1,L)   (kg/s)
C----
C----Calculate convergence in each layer of dry air, compare with expected
C----   dry air mass (AIRX) and then calculate vertical dry-mass fluxes
C----   GAMA(I,J,L)   ==> [I,J,L] ==>  GAMA(I,J,L+1)    (kg/s)
C----
C----Horizontal pressure filter adjusts U & V to reduce error in [PCTM - PS]
C----     U + pressure filter ==> U#,  V + filter ==> V# (temporary)
C----   The pressure filter does nearest neighbor flux (adjusting ALFA/BETA)
C-----------------------------------------------------------------------
C----
C----Note that L->L+1 is upward (decreasing pressure) and that boundaries:
C----GAMA(I,J,1) = GAMA(I,J,LM+1) = 0  no flux across lower/upper boundaries
C----   BETA(I,1,L) = BETA(I,JM+1,L) = 0   no flux at S & N poles
c
c----alternate 1a:
C----   ALFA(1,J,L) = ALFA(IM+1,J,L) is NOT ZERO, but now uses the met-field
c----      values and adds additional flow to achieve uniform mass.
C----Dimensions for ALFA, BETA, GAMA are extended by +1 beyond grid to
C----   allow simple formulation of fluxes in/out of final grid box.
C
C-----GCM input U,V,PSG is ALWAYS of GLOBAL dimensions (IM x JM x LM)
C-----Indices of ALFA, BETA, GAMA, Q & PS are always GLOBAL also
C-----
C-----to do a subset 'WINDOW' calculation define an offset(I0) and size(IM<IPAR)
C
C-----then the tracer arrays(STT), and diagnostics are local (w.r.t. WINDOW)
C-----the ALFA,BETA,GAMA,P,...  need to be remapped to the WINDOW
c----->>>>>>>>this has not been done yet<<<<<<<<<<<<<<<<<
C-----------------------------------------------------------------------
C
      G100   = 100._r8 /G0
      ZDTW   = 1._r8 /DTWIND
      LSP    = .true.
      LNP    = .true.
      LEW    = .true.

c Set for local use (IM, JM, LM are removed from CTM3)
      IM = IPAR
      JM = JPAR
      LM = LPAR

c>>>>>>>>>>>temp fix:  do all the EPZ filtering on V, T, Q, P here
c>>>>>>>>>>>    do the U fix for EPZs using std sub, but not the poles
c
c ------- average  T(I,J,L) & Q(I,J,L) over extended polar zones.
c ------- ALSO NEED TO average T, and fluxes, and BL-ht, etc.
c ------- assume that IMEPZ(1 & JM) = IM to ensure polar box averaging
        call EPZ_TQ(T,Q,XYZA,XYZB,P,IMEPZ,IPAR,JPAR,LPAR,IM,JM,LM)
c ------- average U's and V's in EPZ_s
        call EPZ_UV(U,V,            IMEPZ,IPAR,JPAR,LPAR,IM,JM,LM)
c ------- finally average P's in EPZs
        call EPZ_P(P,               IMEPZ,IPAR,JPAR, IM,JM)
c<<<<<<<<<<<<<<

c------Assumes EPZ filtering of the met fields already been done in pwind.f

c---------SUMAQ(I,J) & SUMAD(I,J): column integral of dry-ar & water (kg)------
c---------AIRQKG(I,J,L): amount of water in each grid box (kg)----------
      do J = 1,JM
      do I = 1,IM
        SUMAQ(I,J) = 0._R8
      enddo
      enddo

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIRQKG(I,J,L) = Q(I,J,L)*(XYZA(I,J,L)+P(I,J)*XYZB(I,J,L))
        SUMAQ(I,J) = SUMAQ(I,J) + AIRQKG(I,J,L)
      enddo
      enddo
      enddo

c===ALFA, BETA, and GAMA are dry-air mass fluxes, water transport not included.

c---define ALFA's from U's and water vapor(Q) - NOW include polar boxes
c---    wrap ALFA at dateline
      do L = 1,LM
       do J = 1,JM
        do I = 2,IM
            AIRQAV  = 0.5_r8 * (Q(I-1,J,L) + Q(I,J,L))
          ALFA(I,J,L) = U(I,J,L) * (1._r8 - AIRQAV) * G100
        enddo
            AIRQAV  = 0.5_r8 * (Q(IM,J,L) + Q(1,J,L))
          ALFA(1,J,L) = U(1,J,L) * (1._r8 - AIRQAV) * G100
          ALFA(IM+1,J,L) = ALFA(1,J,L)   !wrap ALFA at dateline
       enddo
      enddo

c---define BETA's from V's and water vapor(Q), BETA = 0 over-the-pole
      do L = 1,LM
       do J = 2,JM
        do I = 1,IM
            AIRQAV  = 0.5_r8 * (Q(I,J-1,L) + Q(I,J,L))
          BETA(I,J,L) = V(I,J,L) * (1._r8 - AIRQAV) * G100
        enddo
       enddo
        do I = 1,IM
         BETA(I,1,L)    = 0._r8
         BETA(I,JM+1,L) = 0._r8
        enddo
      enddo

c===NB, POLAR CAP is separate from the Extended Polar Zones for met fields
c---Alternate 2:  separate the over-the-pole pie-boxes (NQ = 2*JM)
c---   calculate mass residuals after met-field BETA & ALFA, apply ALFA corr.
c---   to ensure over-the-pole boxes the same air mass but do not connect

c----calculate ALFAs in polar boxes to redistribute divergence uniformly
c---  add this onto exisiting ALFA's
      do L = 1,LM
       do J = 1,JM,JM-1     ! just do the poles
        DALFAI(:) = 0._r8
        do I = 1,IM
c--- DELB = excess mass accumulation in polar=pie boxes (I,1,L) and (I,JM,L)
          DELB = BETA(I,J,L)-BETA(I,J+1,L) + ALFA(I,J,L)-ALFA(I+1,J,L)
          DELI = DELB / float(IM)
          DEL0 = 0.5_r8*(DELB - DELI)
         do II=1,IM
          III = mod (I+II-1, IM) + 1    ! this loops only 1:IM, not IM+1
          DALFAI(III) = DALFAI(III) + DEL0 - DELI*float(II-1)
         enddo
        enddo
        do I = 1,IM
          ALFA(I,J,L) = ALFA(I,J,L) + DALFAI(I)
        enddo
         ALFA(IM+1,J,L) = ALFA(1,J,L)
       enddo
      enddo

c----calculate new DRY-air mass in each box expected at end of time step
c----N.B. presume meridion pipe flow, w/correct BETA(,1,) & BETA(,JM+1,)
      do J = 1,JM
      do I = 1,IM
        SUMAD(I,J) = 0._R8
      enddo
      enddo

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIRNEW(I,J,L) = AIR(I,J,L) + DTWIND *
     &    (ALFA(I,J,L) - ALFA(I+1,J,L) + BETA(I,J,L) - BETA(I,J+1,L))
        SUMAD(I,J) = SUMAD(I,J) + AIRNEW(I,J,L)
      enddo
      enddo
      enddo

c----------------------P FIXER--------------------------------------------------
C----begin P-fixer, where errors in projected surf pressure are mostly
c----corrected by an (dALFA,dBETA) flow that does not change GAMA
c
C-- (1) Define error in surface pressure(kg) expected at end of time step
C-- (2) Filter by error in adjacent boxes, weight by areas, adjust ALFA & BETA
C---  PCTM(I,J)=new CTM wet-air column based on dry-air convergence (Pascals)
C---  MERR(I,J)=pressure-error (kg) between CTM-GCM at new time (before filter)

c<<<<diagnostic only:
        RMSERR1 = 0._r8

      do J = 1,JM
      do I = 1,IM
        PCTM(I,J) = (SUMAD(I,J) + SUMAQ(I,J) - XYA(I,J)) / XYB(I,J)
        MERR(I,J) = (PCTM(I,J) - P(I,J)) * AREAXY(I,J) * G100

c<<<<diagnostic only:
        RMSERR1 = RMSERR1 + MERR(I,J)**2

      enddo
      enddo

        NITR   = 5
C-----------------------------------------------------------------------
      call PFILTR (MERR,AX,BX,AREAXY,IPAR,JPAR,IM,JM)
C-----------------------------------------------------------------------
C---------Calculate corrections to ALFA & BETA from AX and BX:----------
      do L = 1,LM
       do J = 1,JM
        do I = 1,IM+1
             IIX = MIN(I,IM)
          ALFA(I,J,L) = ALFA(I,J,L) +
     &             AX(I,J) * XYZB(IIX,J,L) / (XYB(IIX,J) * DTWIND)
        enddo
       enddo
      enddo
c
      do L = 1,LM
       do J = 2,JM
            JJX = J
            if (J+J .gt. JM)  JJX = J-1
        do I = 1,IM
            BETA(I,J,L) = BETA(I,J,L)
     &           + BX(I,J) * XYZB(I,JJX,L) / (XYB(I,JJX) * DTWIND)
        enddo
       enddo
      enddo

c===POLAR CAP fixes again
c---Alternate 2:  separate the over-the-pole pie-boxes (NQ = 2*JM)
c         BETA(I,1,L)    = 0._r8
c         BETA(I,JM+1,L) = 0._r8
c----re-calculate ALFAs in polar boxes to redistribute divergence uniformly
c---  add this onto exisiting ALFA's
      do L = 1,LM
       do J = 1,JM,JM-1     ! just do the poles
        DALFAI(:) = 0._r8
        do I = 1,IM
c--- DELB = excess mass accumulation in polar=pie boxes (I,1,L) and (I,JM,L)
          DELB = BETA(I,J,L)-BETA(I,J+1,L) + ALFA(I,J,L)-ALFA(I+1,J,L)
          DELI = DELB / float(IM)
          DEL0 = 0.5_r8*(DELB - DELI)
         do II=1,IM
          III = mod (I+II-1, IM) + 1    ! this loops only 1:IM, not IM+1
          DALFAI(III) = DALFAI(III) + DEL0 - DELI*float(II-1)
         enddo
        enddo
        do I = 1,IM
          ALFA(I,J,L) = ALFA(I,J,L) + DALFAI(I)
        enddo
         ALFA(IM+1,J,L) = ALFA(1,J,L)
       enddo
      enddo
C--------------------------end of pressure fixer-----------------------

c----RE-calculate new airmass in each box expected at end of time step
c----N.B. presume meridion pipe flow, w/correct BETA(,1,) & BETA(,JM+1,)

      do J = 1,JM
      do I = 1,IM
        SUMAD(I,J) = 0._r8
      enddo
      enddo

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIRNEW(I,J,L) = AIR(I,J,L) + DTWIND *
     &    (ALFA(I,J,L) - ALFA(I+1,J,L) + BETA(I,J,L) - BETA(I,J+1,L))
        SUMAD(I,J) = SUMAD(I,J) + AIRNEW(I,J,L)
      enddo
      enddo
      enddo

        RMSERR2 = 0._r8

      do J = 1,JM
      do I = 1,IM
        PCTM(I,J) = (SUMAD(I,J) + SUMAQ(I,J) - XYA(I,J)) / XYB(I,J)
        MERR(I,J) = (PCTM(I,J) - P(I,J)) * AREAXY(I,J) * G100

c<<<<diagnostic only:
        RMSERR2 = RMSERR2 + MERR(I,J)**2

      enddo
      enddo

c<<<<diagnostic only:
        RMSERR1 = RMSERR1/float(IM*JM)
        RMSERR2 = RMSERR2/float(IM*JM)
       write(6,'(a,1p,e10.3,a,e10.3)') ' rms P-err (kg):',sqrt(RMSERR1)
     &                        ,'    after P-fixer:', sqrt(RMSERR2)



C---------GAMA:  redistribute the new dry-air mass consistent with the new CTM
C---------   surface pressure(PCTM), rigid upper b.c., no change in PCTM
C---------AIRX(I,J,L) = dry-air mass expected, based on PCTM and eta-levels
c---Polar Cap OK, because of DYN2V meridian pipe flow, have set BETA(J=1 & JM+1)

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIRX(I,J,L)= XYZA(I,J,L) + PCTM(I,J)*XYZB(I,J,L) - AIRQKG(I,J,L)
      enddo
      enddo
      enddo

c---fix GAMA = 0 at top-of-atmos, ensure GAMA = 0 at surface (numerical noise)
      do J = 1,JM
      do I = 1,IM
        GAMA(I,J,LM+1) = 0._r8
        GAMA(I,J,1)    = 0._r8
      enddo
      enddo

      do L = LM,2,-1
      do J = 1,JM
      do I = 1,IM
        GAMA(I,J,L) = GAMA(I,J,L+1) - (AIRNEW(I,J,L) - AIRX(I,J,L))
      enddo
      enddo
      enddo

      do L = 2,LM
      do J = 1,JM
      do I = 1,IM
        GAMA(I,J,L)  = GAMA(I,J,L) * ZDTW
      enddo
      enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine EPZ_UV(U,V,IMEPZ,ID,JD,LD,IM,JM,LM)
c-----------------------------------------------------------------------
c   redistributes the U,V fluxes to smooth over Extended Polar Zones
c        U and V are "wet-air" mass fluxes in/out of box.
c        U(I,J,K)   ==> [I,J,K] ==>  U(I+1,J,K)   (kg/s)
c        V(I,J,K)   ==> [I,J,K] ==>  V(I,J+1,K)   (kg/s)
c        IMEPZ(J)   = # of boxes in extended zones
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      integer,intent(in) ::  ID,JD,LD, IM,JM,LM
      integer,intent(in) ::  IMEPZ(JD)
      real(r8),dimension(ID,JD,LD),intent(inout) :: U, V

      real(r8) :: UVNET,SUMV,ZIMZ
      integer :: I,J,L,IEPZ,IMZ,IEND
c-----------------------------------------------------------------------

c---average V's along poleward side of an Extended Polar Zone
      do J = 2,JM
        if (2*J .LE. JM)  then
          IMZ = IMEPZ(J)      ! S.Hem.
         else
          IMZ = IMEPZ(J-1)    ! N.Hem.
        endif
        if (IMZ .gt. 1)  then
          ZIMZ     = 1._r8 / real(IMZ, r8)
          do L = 1,LM
          do IEPZ = 1,IM,IMZ      ! # EPZ_s = IM/IMZ
            SUMV   = 0._r8
            do I = IEPZ,IEPZ+IMZ-1
              SUMV = SUMV + V(I,J,L)
            enddo
            SUMV   = SUMV * ZIMZ
            do I = IEPZ,IEPZ+IMZ-1
              V(I,J,L) = SUMV
            enddo
          enddo
          enddo
        endif   ! check on IMZ
      enddo     ! end of J loop

      do J = 2,JM-1
        IMZ = IMEPZ(J)
        if (IMZ .GT. 1)  then
          ZIMZ      = 1._r8 /real(IMZ, r8)
          do L = 1,LM
C----IEPZ loops stride thru EPZs: 1, 1+IMZ, 1+2*IMZ,..,
c----calculate net UV convergence and then internal U's
          do IEPZ = 1,IM,IMZ
              IEND  = mod(IEPZ+IMZ,IM)
              UVNET = U(IEPZ,J,L) - U(IEND,J,L)
            do I = IEPZ,IEPZ+IMZ-1       !I all boxes in a single EPZ
              UVNET = UVNET + V(I,J,L) - V(I,J+1,L)
            enddo
              UVNET = UVNET * ZIMZ           ! net UV convergence per grid box
            do I = IEPZ+1,IEPZ+IMZ-1
              U(I,J,L) = U(I-1,J,L) + V(I-1,J,L) - V(I-1,J+1,L) - UVNET
            enddo
          enddo
          enddo
        endif
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine EPZ_TQ(T,Q,XYZA,XYZB,P,IMEPZ,ID,JD,LD,IM,JM,LM)
c-----------------------------------------------------------------------
c   average Q(I,J,L) and T(I,J,L) over extended polar zones.
c   must have UNAVERAGED P(I,J) from met fields for this to work accurately
c   ****must have IMEPZ(1 & JM) = IM
c         T(I,J,L)   = temperature of grid box (C ?)
c         Q(I,J,L)   = specific humidity of grid box (kg H2O / kg wet air)
c         IMEPZ(J)   = # of boxes in extended zones
c         wet air mass kg (I,J,L) = XYZA() + XYZB() * P()
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      integer,intent(in) ::  ID,JD,LD, IM,JM,LM
      integer,intent(in) ::  IMEPZ(JD)
      real(r8),dimension(ID,JD,LD),intent(in)    :: XYZA, XYZB
      real(r8),dimension(ID,JD),intent(in)       :: P
      real(r8),dimension(ID,JD,LD),intent(inout) :: Q, T

      real(r8) :: SUMA,SUMT,SUMQ,AWET
      integer :: I,J,L,IEPZ,IMZ
c-----------------------------------------------------------------------
      do L = 1,LM
      do J = 2,JM-1
          IMZ  = IMEPZ(J)
        if (IMZ .GT. 1)  then
c---------average t,Q over EPZ_s------------------------------------------
          do IEPZ = 1,IM,IMZ                        ! # EPZ_s = IM/IMZ
              SUMA = 0._r8
              SUMQ = 0._r8
              SUMT = 0._r8
            do I = IEPZ,IEPZ+IMZ-1
              AWET = XYZA(I,J,L) + P(I,J) * XYZB(I,J,L)
              SUMT = SUMT + T(I,J,L)*AWET
              SUMQ = SUMQ + Q(I,J,L)*AWET
              SUMA = SUMA + AWET
            enddo
              SUMQ = SUMQ / SUMA
              SUMT = SUMT / SUMA
            do I = IEPZ,IEPZ+IMZ-1
              Q(I,J,L) = SUMQ
              T(I,J,L) = SUMT
            enddo
          enddo
        endif
      enddo
c---------average t,Q over poles------------------------------------------
          do J = 1,JM,JM-1
              SUMA = 0._r8
              SUMQ = 0._r8
              SUMT = 0._r8
            do I = 1,IM
              AWET = XYZA(I,J,L) + P(I,J) * XYZB(I,J,L)
              SUMT = SUMT + T(I,J,L)*AWET
              SUMQ = SUMQ + Q(I,J,L)*AWET
              SUMA = SUMA + AWET
            enddo
              SUMQ = SUMQ / SUMA
              SUMT = SUMT / SUMA
            do I = 1,IM
              Q(I,J,L) = SUMQ
              T(I,J,L) = SUMT
            enddo
          enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine EPZ_P(P,IMEPZ,ID,JD, IM,JM)
c-----------------------------------------------------------------------
c   average P(I,J) over extended polar zones.
c   must be done last in case grid size needed to average other quantites
c   ****must have IMEPZ(1 & JM) = IM
c         P(I,J)   = pressure (hPa) of grid box
c         IMEPZ(J)   = # of boxes in extended zones
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      integer,intent(in) ::  ID,JD, IM,JM
      integer,intent(in) ::  IMEPZ(JD)
      real(r8),dimension(ID,JD),intent(inout) :: P
c
      real(r8) :: SUMP,ZIMZ
      integer :: I,J,IEPZ,IMZ
c-----------------------------------------------------------------------
      do J = 2,JM-1
          IMZ  = IMEPZ(J)
          ZIMZ = 1._r8 / real(IMZ, r8)
        if (IMZ .GT. 1)  then
c---------average P over EPZ_s------------------------------------------
          do IEPZ = 1,IM,IMZ
              SUMP = 0._r8
            do I = IEPZ,IEPZ+IMZ-1
              SUMP = SUMP + P(I,J)
            enddo
              SUMP = SUMP * ZIMZ
            do I = IEPZ,IEPZ+IMZ-1
              P(I,J) = SUMP
            enddo
          enddo
        endif
      enddo
          ZIMZ = 1._r8 / real(IM, r8)
          do J = 1,JM,JM-1
c---------average P over poles -------------------------------------------
              SUMP = 0._r8
            do I = 1,IM
              SUMP = SUMP + P(I,J)
            enddo
              SUMP = SUMP * ZIMZ
            do I = 1,IM
              P(I,J) = SUMP
            enddo
          enddo

      return
      end
c
c
c ----------------------------------------------------------------------
      subroutine PFILTR (MERR,ALFAX,BETAX,AXY,ID,JD,IM,JM)
c ----------------------------------------------------------------------
c---new version 5.3+, only global filtering, no local--assume all fields global!c---------FILTER is smoothing presssure error, the error in pressure
c---------between predicted Ps(CTM) and Ps(GCM).
c ------- MERR(ID,JD)        mass error (total kg in each (I,J)-box)
c ------- AXY(ID,JD)         area of grid box (I,J)
c ------- ALFAX(ID+1,JD)     corrections to ALFA based on MERR
c ------- BETAX(ID,JD+1)     corrections to BETA based on MERR
c
      use cmn_precision, only: r8
      implicit none
      integer,intent(in)  ::  ID,JD,IM,JM
      real(r8), intent(in)  ::   AXY(ID,JD)
      real(r8), intent(inout)::  MERR(ID,JD)
      real(r8), intent(out) ::   ALFAX(ID+1,JD),BETAX(ID,JD+1)

      real(r8) :: ERRAVG,SUMAXY,ALFAVG,  BETAXJ(321),ALFAXI(641)
      integer :: I,J

          ALFAX(:,:) = 0._r8
          BETAX(:,:) = 0._r8

c---redistribute mass error uniformly along I-th meridional stripe SPole->NPole
      do I = 1,IM
          ERRAVG = 0._r8
          SUMAXY = 0._r8
        do J = 1,JM
          ERRAVG = ERRAVG + MERR(I,J)
          SUMAXY = SUMAXY + AXY(I,J)
        enddo
          ERRAVG = ERRAVG / SUMAXY
c---calculate perturbation to BETA: calc air to move out of box [J-1]
c---    leave behind the average mass error distributed over meridion
          BETAXJ(1) = 0._r8
c         BETAXJ(JM+1) = 0._r8
        do J = 2,JM
          BETAXJ(J) = BETAXJ(J-1) + MERR(I,J-1) - AXY(I,J-1)*ERRAVG
        enddo
c---add onto any other corrections to BETAX and update MERR
        do J = 2,JM
          BETAX(I,J) = BETAX(I,J) + BETAXJ(J)
        enddo
        do J = 1,JM
          MERR(I,J) = MERR(I,J) + (BETAX(I,J)-BETAX(I,J+1))
        enddo
      enddo

c---redistribute mass error uniformly along J-th latitude belt (loop)
      do J = 1,JM
          ERRAVG = 0._r8
          SUMAXY = 0._r8
        do I = 1,IM
          ERRAVG = ERRAVG + MERR(I,J)
          SUMAXY = SUMAXY + AXY(I,J)
        enddo
          ERRAVG = ERRAVG / SUMAXY
c---calculate perturbation to ALFA: calc air to move out of box [I-1]
          ALFAXI(1) = 0._r8
          ALFAXI(IM+1) = 0._r8
        do I = 2,IM
          ALFAXI(I) = ALFAXI(I-1) + MERR(I-1,J) - AXY(I-1,J)*ERRAVG
        enddo
c---because this is a loop, remove <ALFAXI> to avoid net rotation
          ALFAVG = 0._r8
        do I = 2,IM
          ALFAVG = ALFAVG + ALFAXI(I)
        enddo
          ALFAVG = ALFAVG / float(IM)
        do I = 1,IM+1
          ALFAXI(I) = ALFAXI(I) - ALFAVG
        enddo
c---add onto any other corrections to ALFAX and update MERR
        do I = 1,IM+1
          ALFAX(I,J) = ALFAX(I,J) + ALFAXI(I)
        enddo
        do I = 1,IM
          MERR(I,J) = MERR(I,J) + (ALFAX(I,J)-ALFAX(I+1,J))
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine CFLADV (DTOPS,NADV)
c-----------------------------------------------------------------------
c---global Lifshitz (not CFL) time-step limiter for advection step:
c     NADV = number of steps per Op-Split time(DTOPS)
c     follows the sequence W-V-U (now the core Op-Split sequence)
c     does not allow box to fall below (1 - CFLLIM)*AIR at any step.
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, LPAR
      use cmn_ctm, only: ALFA, BETA, GAMA, AIR, AIRX, CFLLIM
      use cmn_diag, only: MNADV
      use cmn_met, only: U, V, Q, T, P
      implicit none
C-----------------------------------------------------------------------
      real(r8), intent(in)  ::  DTOPS
      integer, intent(out) :: NADV
c
      real(r8), dimension(IPAR,JPAR,LPAR) :: AIRMN, AIR2, AIR3, AIR4
      real(r8) :: AD1MIN, DTFIX, ARATIO, DTMAX, DTDYN, AIRMIN, CNST
      integer :: I,J,L
      integer :: IM, JM, LM
c-----------------------------------------------------------------------
      DTFIX  = 5._r8*60._r8     ! do test steps with 5 min time step

c Set for local use (IM, JM, LM are removed from CTM3)
      IM = IPAR
      JM = JPAR
      LM = LPAR

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIRMN(I,J,L)  = min(AIR(I,J,L),AIRX(I,J,L))
      enddo
      enddo
      enddo

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIR2(I,J,L) = AIRMN(I,J,L) + DTFIX*(GAMA(I,J,L)-GAMA(I,J,L+1))
      enddo
      enddo
      enddo

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIR3(I,J,L) = AIR2(I,J,L) + DTFIX*(ALFA(I,J,L)-ALFA(I+1,J,L))
      enddo
      enddo
      enddo

      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AIR4(I,J,L) = AIR3(I,J,L) + DTFIX*(BETA(I,J,L)-BETA(I,J+1,L))
      enddo
      enddo
      enddo

      ARATIO  = 1._r8
      CNST   = 1._r8 - 0.5_r8*DTFIX/DTOPS
c                      ^-- should be smaller than CFLLIM (0.95)
      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
        AD1MIN = CNST * AIRMN(I,J,L)
        AIRMIN = min(AIR2(I,J,L),AIR3(I,J,L),AIR4(I,J,L),AD1MIN)
        ARATIO = min(ARATIO,AIRMIN/AIR(I,J,L))
      enddo
      enddo
      enddo

c-------- adjust time step ---------------------------------------------
      DTMAX  = DTFIX * CFLLIM / (1._r8 - ARATIO)
      NADV   = int(DTOPS/DTMAX) + 1
      !// Possible override for the number of global transport steps
      !// per NOPS. In T42, NADV is usually 1 (1 hour step) while
      !// limiting by 2 (30 min) makes CTM3 more similar to CTM2
      !NADV = max(NADV, 2)
      DTDYN  = DTOPS / real(NADV, r8)

      !// Sum up to find average NADV
      MNADV   = MNADV + NADV

      return
      end
