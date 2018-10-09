
c-----------------------------------------------------------------------
      Program PMEANC
c-----------------------------------------------------------------------
      include 'cmn_h.f'
      include 'cmn_w.f'
      include 'cmn_s.f'

      real*8  ZNN,ZIJ, TL(LPAR),TSUM,ZOFL(LPAR+1),DELZ, SUMXYZ
c---day number, starting and ending day
      integer NDAY, NDAYI, NDAYE, NMET, NN,I,J,L
      character*80 INFILE1
      logical  LYEAR

C------- input and initialization, includes Chem input -----------------
      PMEAN(:,:) = 0.d0
      TL(:) = 0.d0
      NN = 0
      ZIJ = 1.d0 / dble(IPAR*JPAR)

      call INPUT (NDAYI, NDAYE)
      call WIND (0,1,1)

c     INFILE1 = 'P_mean'//MODEL(4:10)//'.dat'
      INFILE1 = 'P_mean'//MODEL(3:10)//'.dat'
      if (LGAUGRD .and. MODEL(4:6).eq.'1x1')  then
        write (INFILE1(7:9),'(1HN,I2)') JPAR/2
        write(6,*) 'EC1x1 in Gauss grid'
      endif

      do NDAY = NDAYI,NDAYE-1

        do NMET = 1,NRMETD
          call WIND (1,NDAY,NMET)     ! read in met fields
          NN = NN + 1
          do J = 1,JPAR
          do I = 1,IPAR
            PMEAN(I,J) = PMEAN(I,J) + P(I,J)
          enddo
          enddo
          do L = 1,LPAR
            TSUM = 0.d0
            do J = 1,JPAR
            do I = 1,IPAR
              TSUM = TSUM + T(I,J,L)
            enddo
            enddo
            TL(L) = TL(L) + ZIJ*TSUM
          enddo
        enddo
        IDAY = NDAY+1
        call CALENDR (IYEAR,IDAY, LLPYR,LFIXMET,MYEAR,
     &      JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)
      enddo

      ZNN = 1.d0 / dble(NN)
      do J = 1,JPAR
      do I = 1,IPAR
        PMEAN(I,J) = ZNN * PMEAN(I,J)
      enddo
      enddo
      do L = 1,LPAR
        TL(L) = ZNN * TL(L)
      enddo

      ZOFL(1)  = 0.d0
      do L = 2,LPAR+1
        DELZ  = -29.36d0 * TL(L-1) * log(ZEDG(L)/ZEDG(L-1))
        ZOFL(L) = ZOFL(L-1) + DELZ
      enddo

      open (10, file=INFILE1, form='formatted',status='new')
      write(10,'(2A,2(I3,A),2I5,A,2I4)') MODEL(3:10)
     &   ,': mean PSF ((P(I,J),I=1,'
     &   ,IPAR,'),J=1,',JPAR,'), YEAR',IYEAR,JYEAR,'  DAY',NDAYI,JDAY-1
      write(10,'(10f8.2)') ((PMEAN(I,J),I=1,IPAR),J=1,JPAR)
      write(10,'(2A)') MODEL(3:10),
     & ' Area in km^2 (AREAXY(1,J),J=1,JPAR)'
      write(10,'(10f8.4)') (1.d-9*AREAXY(1,J),J=1,JPAR)
      write(10,'(2A)') MODEL(3:10),
     & ' Longitude (edge_pnt deg.) (XDEDG(I),I=1,IPAR+1)'
      write(10,'(10f8.2)') (XDEDG(I),I=1,IPAR+1)
      write(10,'(2A)') MODEL(3:10),
     & ' Latitude (edge_pnt deg.) (YDEDG(J),J=1,JPAR+1)'
      write(10,'(10f8.2)') (YDEDG(J),J=1,JPAR+1)
      write(10,'(2A)') MODEL(3:10),
     & ' Longitude (mid_pnt deg.) (XDGRD(I),I=1,IPAR)'
      write(10,'(10f8.2)') (XDGRD(I),I=1,IPAR)
      write(10,'(2A)') MODEL(3:10),
     & ' Latitude (mid_pnt deg.) (YDGRD(J),J=1,JPAR)'
      write(10,'(10f8.2)') (YDGRD(J),J=1,JPAR)
      write(10,'(2A)') MODEL(3:10),
     & ' Altitude (edge km)'
      write(10,'(10f8.3)') (1.d-3*ZOFL(L),L=1,LPAR+1)
      write(10,'(2A)') MODEL(3:10),
     & ' mean TEM (T(L),L=1,LPAR)'
      write(10,'(10f8.2)') (TL(L),L=1,LPAR)
      write(10,'(2A)') MODEL(3:10),
     & ' ETAA: L=1,LPAR+1'
      write(10,'(10f8.3)') (ETAA(L),L=1,LPAR+1)
      write(10,'(2A)') MODEL(3:10),
     & ' ETAB: L=1,LPAR+1'
      write(10,'(10f8.5)') (ETAB(L),L=1,LPAR+1)
      close(10)

      open (10, file='P_meanT319L60.unf', form='unformatted'
     &        ,status='new')
      write(10) PMEAN,AREAXY,XDEDG,YDEDG,XDGRD,YDGRD,ZOFL,TL,ETAA,ETAB
      stop
      end


c-----------------------------------------------------------------------
      subroutine INPUT(NDAYI,NDAYE)
c-----------------------------------------------------------------------
c--- read std input(5) file that contains basic data for CTM run
c--- read  file names for additional data sets and call the readers:
c---     tracer specifications and supporting data
c---     met field data and file/directory info
c-----------------------------------------------------------------------
      include 'cmn_h.f'
      include 'cmn_w.f'
      integer, intent(out) ::  NDAYI, NDAYE

      character*80 INFILE1,INFILE2,INFILE3
      real*8   XLNG(2), YLAT(2), GM00Z
      integer  I,J,L,M,N, IMX,JMX,LMX, NOPSTL, K
      integer  POLAVG(25)
      integer  II, JMPOLR
      logical  LYEAR, LISLSCP2

c-------------------------------PRIMARY CTM DATA (unit=5)----------------------
c---READ input on unit=5
      read(5,'(a80)') RTITLE
        write(6,'(a10,a80)') 'CTM run>>>', RTITLE
        write(6,'(2a)') 'Met mdl>>>',MODEL
      read(5,*)
      read(5,'(i5)') IYEAR
      read(5,'(2i5)') NDAYI
      read(5,'(2i5)') NDAYE
        write(6,'(a,3i8)') ' base yr, day begin/end:',IYEAR,NDAYI,NDAYE
        GMTAU  = 0.d0
      read(5,*)
      read(5,'(l5)') LLPYR
      read(5,'(l5)') LFIXMET
      read(5,'(l5)') LCLDQMD
      read(5,'(l5)') LCLDQMN
      read(5,'(l5)') LCLDRANA
      read(5,'(l5)') LCLDRANQ
      read(5,'(i5)') RANSEED
      read(5,'(3i5)') NRMETD,NROPSM,NRCHEM
      read(5,'(i5)') LMTSOM
        LMTSOM = min(3,max(1, LMTSOM))
      read(5,'(f5.3)') CFLLIM
      read(5,*)
      read(5,'(a)') MPATH1
      read(5,'(a)') MPATH2
      read(5,'(a)') MFILE3

c---update the calendar to the new, upcoming day = 1st day of CTM run
      call CALENDR (IYEAR,NDAYI, LLPYR,LFIXMET,MYEAR,
     &      JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)

      read(5,'(3i5)') IMX,JMX,LMX
      read(5,'(3i5)') IM, JM, LM
      read(5,'(i5)') JMPOLR
      read(5,'(l5)') LGAUGRD
      read(5,'(f5.1)')  GM00Z

c---read in eta levels of met-field grid (collapse later)
c---if LM < LPARW (met field layers) use LMMAP to collapse met-layers on CTM L
        read(5,*)
      do L = 1,LPARW+1
        read(5,'(5x,i5,f20.10,f17.10)') LMMAP(L),ETAAW(L),ETABW(L)
      enddo

C---SET UP GRID (I,J,L)
      call SET_GRID (GM00Z,JMPOLR)

      return
      end


c-----------------------------------------------------------------------
      subroutine SET_GRID (GM0000,JMPOLR)
c-----------------------------------------------------------------------
c---calls subs: LABELG and GAUSST
c        plus (EC spectral subs) LEGGEN, FFT_DDSS, FFT_SET
c
c---NEW general grid set up, global only, handles EC and GISS grids
c---Latitude (J=1:JM) grid:
c      LGAUGRD = .true. - use Gaussian grid for latitude, else uniform grid
c      JMPOLR (regualr grid only)  =0 =>polar box same delta-lat as others,
c                                  =1 =>polar box is halfsize (eg GISS JM=46)
c---Longitude (I=1:IM) grid:   always regular, uniform spacing
c      GM0000 = I-coord of the Greenwich Merid: 
c          1.5 => mid-box[1] on GM
c          IPARW/2+1 => left-box[1] at Dateline
c---currently for EC, GM0000=1.5;   for GISS4x5, GM0000=37.0 & JMPOLR=1
c
c---Pressure/Altitude (L=1:LM) grid:  based on eta coords
c
c-----------------------------------------------------------------------
      include 'cmn_h.f'
      include 'cmn_s.f'
      include 'cmn_w.f'
      real*8,  intent(in):: GM0000
      integer, intent(in):: JMPOLR

      integer I,J,L,LW, IMB2,JMB2
      real*8  DEL_JD,DEL_I,DEL_ID


      real*8  SUMXYZ, AIRGLB

      real*8  WGAULAT(JPARW),WGAUWT(JPARW), WYEDG1P1(JPARW+1)
      real*8  YEDG1P1(JPAR+1)
      real*8, dimension(JPARW) :: YGRDW, YDGRDW
      real*8, dimension(JPARW+1) :: YEDGW, YDEDGW
      real*8, dimension(IPARW) :: XGRDW, XDGRDW
      real*8, dimension(IPARW+1) :: XEDGW, XDEDGW

c-----------------------------------------------------------------------
          JMB2 = JPARW/2
          IMB2 = IPARW/2
          if (mod(IPARW,2).ne.0) call EXITC(' >>>>IMPARW must be even')

c-----------------------------------------------------------------------
c---LONGITUDE grid is always regular and set by IPARW and GM000
        DEL_ID = 360.d0/dble(IPARW)
      do I = 1,IPARW+1
        XDEDGW(I) = (dble(I) - GM0000) * DEL_ID
      enddo
      do I = 1,IPARW
        XDGRDW(I) = 0.5d0*(XDEDGW(I)+XDEDGW(I+1)) 
      enddo
      do I = 1,IPARW+1
        XEDGW(I) = XDEDGW(I) * CPI180
      enddo
      do I = 1,IPARW
        XGRDW(I) = XDGRDW(I) * CPI180
      enddo
c---for computation XEDGW must be montonic, for labeling shift to neg. longit.
      do I = 1,IPARW+1
        if (XDEDGW(I) .gt. 180.d0) XDEDGW(I) = XDEDGW(I) - 360.d0
      enddo
      do I = 1,IPARW
        if (XDGRDW(I) .gt. 180.d0) XDGRDW(I) = XDGRDW(I) - 360.d0
      enddo
c---set dateline box
      IDTLN = mod(int(GM0000)+IPAR/2-1, IPAR) + 1

c-----------------------------------------------------------------------
c---LATITUDE grid can be Gauss-pt, or regular, incl. half-size boxes at pole
c---   WYEDG1P1 = sin(latitude) of box edges (-1 to +1) [temporary]

      if (mod(JPARW,2) .ne. 0) call EXITC(' >>>>JPARW must be even')

      if (MODEL(1:2) .eq. 'EC')  then

        if (.not. LGAUGRD) then
          DEL_JD = 180.d0/dble(JPARW - JMPOLR)
          YDEDGW(1) = -90.d0
          YDEDGW(2) = -90.d0 + DEL_JD/dble(1 + JMPOLR)
          do J = 3,JPARW
            YDEDGW(J) = YDEDGW(2) + dble(J-2) * DEL_JD
          enddo
          YDEDGW(JPARW+1) = 90.d0
          do J = 1,JPARW+1
            YEDGW(J) = YDEDGW(J)*CPI180
            WYEDG1P1(J) = sin(YDEDGW(J)*CPI180)
          enddo
          do J = 1,JPARW
            YGRDW(J) = asin(0.5d0*(WYEDG1P1(J)+WYEDG1P1(J+1)))
            YDGRDW(J) = ZPI180 * YGRDW(J)
          enddo
          do J = 1,JPARW+1
            WGLYE(J) = Asin(WYEDG1P1(J))
          enddo
          do J = 1,JPARW
            WGLYG(J) = YGRDW(J)
          enddo

        else

          call GAUSST (JPARW,WGAULAT,WGAUWT)

          WYEDG1P1(1) = -1.d0
          do J = 2,JMB2
            WYEDG1P1(J) = WYEDG1P1(J-1) + WGAUWT(J-1)
          enddo
          WYEDG1P1(JMB2+1) = 0.d0
          do J = JMB2+2,JPARW+1
            WYEDG1P1(J) = -WYEDG1P1(JPARW+2-J)
          enddo

          do J = 1,JPARW+1
            WGLYE(J) = Asin(WYEDG1P1(J))
          enddo
          do J = 1,JPARW
            WGLYG(J) = Asin(WGAULAT(J))
          enddo

          do J = 1,JPARW+1
            YEDGW(J) = asin(WYEDG1P1(J))
            YDEDGW(J) = ZPI180 * YEDGW(J)
          enddo
          do J = 1,JPARW
            YGRDW(J) = asin(WGAULAT(J))
            YDGRDW(J) = ZPI180 * YGRDW(J)
          enddo

        endif

c-----------------------------------------------------------------------
c---EC: set up Spectral to Grid transform functions: Assoc-Leg Polys
        do J = 1,JMB2
          call LEGGEN(ALP(1,J),YGRDW(J),NTRUNW,NMMAX)
        enddo
        do J = 1,JMB2
          call LEGGEN(ALPV(1,J),YEDGW(J+1),NTRUNW,NMMAX)
        enddo
c---transforms from Vorticiy and Divergence to U and V
          call FFT_DDSS (NTRUNW,A0,NMMAX,DD,SS)
c---FFT transforms, use grid
          call FFT_SET (TRIG,IFAX,IPARW)
c---trig functions for U (shift to boundaries)
          DEL_I = C2PI/dble(IPARW)
        do I = 1,IMB2
          TRIGU(2*I-1) = cos(DEL_I*(dble(I)-0.5d0))
          TRIGU(2*I  ) = sin(DEL_I*(dble(I)-0.5d0))
        enddo

      else
c---not EC grids
        DEL_JD = 180.d0/dble(JPARW - JMPOLR)
        YDEDGW(1) = -90.d0
        YDEDGW(2) = -90.d0 + DEL_JD/dble(1 + JMPOLR)
        do J = 3,JPARW
          YDEDGW(J) = YDEDGW(2) + dble(J-2) * DEL_JD
        enddo
        YDEDGW(JPARW+1) = 90.d0
        do J = 1,JPARW+1
          YEDGW(J) = YDEDGW(J)*CPI180
          WYEDG1P1(J) = sin(YDEDGW(J)*CPI180)
        enddo
        do J = 1,JPARW
          YGRDW(J) = asin(0.5d0*(WYEDG1P1(J)+WYEDG1P1(J+1)))
          YDGRDW(J) = ZPI180 * YGRDW(J)
        enddo
      endif
c-----------------------------------------------------------------------

      do J = 1,JPARW
        do I = 1,IPARW
          AREAXYW(I,J)= A0*A0*(XEDGW(I+1)-XEDGW(I))*
     &                 (WYEDG1P1(J+1)-WYEDG1P1(J))
        enddo
      enddo

c--- Define lat-lon of ctm grid and weightings if horizontal resolution degraded

      call DBLDBL(YGRD,YDGRD,XGRD,XDGRD,YEDG,YDEDG,XEDG,XDEDG,
     &            YGRDW,YDGRDW,XGRDW,XDGRDW,YEDGW,YDEDGW,XEDGW,XDEDGW,
     &            ZDEGI,ZDEGJ,IMAP,JMAP,
     &            IPARW,JPARW,IPAR,JPAR,IDGRD,JDGRD,LDEG)

        YEDG1P1(1) = -1.d0
      do J = 2,JPAR
        YEDG1P1(J) = sin(YEDG(J))
      enddo
        YEDG1P1(JPAR+1) = 1.d0

c-----------------------------------------------------------------------
c---PRESSURE/ALTITUDE grid, remap the global wind grid (W)
c---if LM < LPARW (met field layers) use LMMAP to collapse met-layers on CTM L

c---calculate new ETAA & ETAB for CTM(1:LM) with remapped vertical grid(1:LPARW)
        do LW = 1,LPARW+1
          ETAAW(LW) = 0.01d0*ETAAW(LW)
        enddo
          ETAA(1) = ETAAW(1)
          ETAB(1) = ETABW(1)
        do LW = 1,LPARW
          L = LMMAP(LW)
          ETAA(L+1) = ETAAW(LW+1)
          ETAB(L+1) = ETABW(LW+1)
        enddo
          LM = LMMAP(LPARW)

c---calculate weighting for LW=1:LPARW to merged CTM layer L=1:LM
        do LW = 1,LPARW
          L = LMMAP(LW)
          XLMMAP(LW) = ((ETAAW(LW)-ETAAW(LW+1)) + 
     &        1.d3*(ETABW(LW)-ETABW(LW+1)) ) / ((ETAA(L)-ETAA(L+1)) +
     &                  1.d3*(ETAB(L)-ETAB(L+1)) )
        enddo

c---DIAGNOSTICS AND WEIGHTING FACTORS

c>>>>>this new gridding does not yet deal with nested, doubled grids!
c---AREAS:  use the IM,JM rather than IPARW,JPARW in case of re-grid
      do J = 1,JM
        DISTY(J) = A0*(YEDG(J+1)-YEDG(J))
      enddo
      do J = 1,JM
        DISTX(J) = 
     &   A0 * (XEDG(2)-XEDG(1)) * 0.5d0 * (cos(YEDG(J))+cos(YEDG(J+1)))
      enddo
      do J = 1,JM
       do I = 1,IM
        AREAXY(I,J)= A0*A0*(XEDG(I+1)-XEDG(I))*(YEDG1P1(J+1)-YEDG1P1(J))
       enddo
      enddo
      SUMXYZ  = 0.d0
      do J = 1,JM
       do I = 1,IM
         SUMXYZ = SUMXYZ + AREAXY(I,J)
       enddo
      enddo
      write(6,'(A,1PE12.5)') 'total surface area of T319 ',SUMXYZ
      SUMXYZ  = 0.d0
      do J = 1,JM
        SUMXYZ = SUMXYZ + AREAXY(1,J)
      enddo
      write(6,'(A,1PE12.5)') 'total surface area of T319 ',640d0*SUMXYZ

c---air mass factors for each grid box:
c---    air mass(I,J,L) in kg = XYZA() + XYZB()*surface pressure (mbar!!)
            SUMXYZ  = 0.d0
      do L = 1,LM
        do J = 1,JM
          do I = 1,IM
            XYZA(I,J,L) = (ETAA(L)-ETAA(L+1))*AREAXY(I,J)*1.d2/G0
            XYZB(I,J,L) = (ETAB(L)-ETAB(L+1))*AREAXY(I,J)*1.d2/G0
            SUMXYZ  = SUMXYZ + XYZB(I,J,L)
          enddo
        enddo
      enddo
        do J = 1,JM
          do I = 1,IM
            XYZA(I,J,LM+1) = ETAA(LM+1)*AREAXY(I,J)*1.d2/G0
            XYZB(I,J,LM+1) = ETAB(LM+1)*AREAXY(I,J)*1.d2/G0
          enddo
        enddo
c---column air mass factors XYA, XYB do not include > model top
      do J = 1,JM
        do I = 1,IM
          XYA(I,J)   = 0.d0
          XYB(I,J)   = 0.d0
          do L = 1,LM
            XYA(I,J) = XYA(I,J) +XYZA(I,J,L)
            XYB(I,J) = XYB(I,J) +XYZB(I,J,L) 
          enddo
        enddo
      enddo
      
      Call TRUNG8(PMEANW,PMEAN,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD,
     &            JDGRD,IPARW,JPARW,IM,JM,1,1)

C-----areas and air mass
        AREAG = 0.d0
        PMEANG = 0.d0
      do J=1,JM
      do I=1,IM
        AREAG = AREAG + AREAXY(I,J)
        PMEANG = PMEANG + PMEAN(I,J)*AREAXY(I,J)
      enddo
      enddo
        PMEANG = PMEANG/AREAG
        AIRGLB    = 0.d0
      do L=1,LM
      do J=1,JM
      do I=1,IM
        AIRGLB  = AIRGLB + XYZA(I,J,L) + XYZB(I,J,L)*PMEAN(I,J)
      enddo
      enddo
      enddo

c---create new 3-D P-grid for interpolating sources, etc. using PMEAN
c---    PIJL(I,J,L=1:LM+1) is pressure at edges
      do L=1,LM+1
      do J=1,JM
      do I=1,IM
        PIJL(I,J,L) =  ETAA(L) + ETAB(L)*PMEAN(I,J)
      enddo
      enddo
      enddo

c---------------labels:-------------------------------------------------
      do J=1,JM
        call LABELG(YDGRD(J),TLAT(J),1)
      enddo
      do J=1,JM+1
        call LABELG(YDEDG(J),TLATE(J),1)
      enddo
      do I=1,IM
        call LABELG(XDGRD(I),TLNG(I),2)
      enddo
      do I=1,IM+1
        call LABELG(XDEDG(I),TLNGE(I),2)
      enddo     
      do L=1,LM+1
          ZEDG(L) = ETAA(L) + ETAB(L)*1000.d0 
          call LABELG(ZEDG(L),TALTE(L),3)
      enddo
      do L=1,LM
          ZGRD(L) = 0.5d0*(ZEDG(L)+ZEDG(L+1))
          call LABELG(ZGRD(L),TALT(L),3)
      enddo

c-----------------write-out grid data-----------------------------------
c---will need to shift absolute global grid to actual window calc later!
      write(6,*) '  ------------------grid----------------------'
      write(6,'(a,1p,e15.8)') ' Surface presse(hPa): ', PMEANG
      write(6,'(a,1p,e15.8)') ' AIR mass wet (kg)  : ', AIRGLB
      write(6,'(a,1p,e15.8)') ' AREA globe(m^2)    : ', AREAG
      write(6,'(A)') 'Area of grid boxes (m^2) 1:JM '
      write(6,'(1P10E10.3)') (AREAXY(1,J),J=1,JM)

      write(6,'(a)') ' MID-POINT of GRID BOXES'
      write(6,'(a,i5)') '     J=1:',JM
      write(6,'(5(i5,f9.5,f6.1,a5))') 
     &      (J, YGRD(J), YDGRD(J), TLAT(J), J=1,JM)
      write(6,'(a,i5)') '     I=1:',IM
      write(6,'(5(i5,f9.5,f6.1,a5))')
     &      (I, XGRD(I), XDGRD(I), TLNG(I), I=1,IM)
      write(6,'(a,i8)') '  IDATELINE=', IDTLN
      write(6,'(a)')  '   LW   wting  ==> L (lbl)'
      write(6,'(i5,f10.5,i5,a5)')
     &     (LW, XLMMAP(LW),LMMAP(LW),TALT(LMMAP(LW)), LW=1,LPARW)
      write(6,*) ' EDGES of GRID'
      write(6,'(a,i5)') '     J=1:',JM+1
      write(6,'(5(i5,f9.5,f6.1,a5))') 
     &      (J, YEDG(J), YDEDG(J), TLATE(J), J=1,JM+1)
      write(6,'(a,i5)') '     I=1:',IM+1
      write(6,'(5(i5,f9.5,f6.1,a5))') 
     &      (I, XEDG(I), XDEDG(I), TLNGE(I), I=1,IM+1)
      write(6,'(a,i5)') '     L=1:',LM+1
      write(6,'(1x,a1,i2,a2,a5,f9.3,a3)') 
     &      ('L',L-1,'.5', TALTE(L), ZEDG(L),' mb', L=1,LM+1)

      return
      end


c-----------------------------------------------------------------------
      subroutine LABELG(X,TITLX,NDEX)
c-----------------------------------------------------------------------
c---generate char*4 titles for lat(NDEX=1), long(NDEX=2) or vert(NDEX=3)
c---
      implicit none
      real*8, intent(in) ::  X
      character*4, intent(out) :: TITLX
      integer, intent(in) :: NDEX

      character*1, parameter, dimension(2)  :: TITLNS =['N','S']
      character*1, parameter, dimension(2)  :: TITLEW =['E','W']
      character*1, parameter, dimension(12) :: TITLNN =[' ','1','2',
     &            '3','4','5','6','7','8','9','0','-']
      real*8    XABS
      integer   I1000,I0100,I0010,I0001,I100,I010,I001,I10,I01,IXABS

      XABS = abs(X)
c---latitude (NDEX=1)
      if (NDEX.EQ.1) then
       IXABS = int(XABS + 0.5d0)
       I10 = IXABS/10
       I01 = IXABS - 10*I10
       if (I01.eq.0) I01 = 10
       I10 = min(10,I10)
       TITLX(1:1) = TITLNN(1)
       TITLX(2:2) = TITLNN(I10+1)
       TITLX(3:3) = TITLNN(I01+1)
       if (X.ge.0.d0) then
        TITLX(4:4) = TITLNS(1)
       else
        TITLX(4:4) = TITLNS(2)
       endif
      endif
c---longitude (NDEX=2)
      if (NDEX.eq.2) then
       IXABS = int(XABS + 0.5d0)
       I100 = IXABS/100
       I010 = (IXABS - 100*I100)/10
       I001 = IXABS - 100*I100 - 10*I010
       if (I100.gt.0 .AND. I010.eq.0) I010 = 10
       if (I001.eq.0) I001=10
       I100 = min(10,I100)
       TITLX(1:1) = TITLNN(I100+1)
       TITLX(2:2) = TITLNN(I010+1)
       TITLX(3:3) = TITLNN(I001+1)
       if (X.ge.0.d0) then
        TITLX(4:4) = TITLEW(1)
       else
        TITLX(4:4) = TITLEW(2)
       endif
      endif
C---pressure (NDEX=3) - nearest mbar only for now
      if (NDEX.eq.3) then
       if (XABS.gt.1.d0) then
        IXABS = int(XABS + 0.5d0)
        I1000 = IXABS/1000
        I0100 = (IXABS - 1000*I1000)/100
        I0010 = (IXABS - 1000*I1000 - 100*I0100)/10
        I0001 =  IXABS - 1000*I1000 - 100*I0100 - 10*I0010
        I1000 = min(I1000,10)
        if (I1000.gt.0 .AND. I0100.eq.0) I0100 = 10
        if (I0100.gt.0 .AND. I0010.eq.0) I0010 = 10
        if (I0001.eq.0) I0001 = 10
        TITLX(1:1) = TITLNN(I1000+1)
        TITLX(2:2) = TITLNN(I0100+1)
        TITLX(3:3) = TITLNN(I0010+1)
        TITLX(4:4) = TITLNN(I0001+1)
       else
        TITLX(1:1) = TITLNN(1)
        TITLX(2:2) = TITLNN(1)
        TITLX(3:3) = TITLNN(1)
        XABS = 10.d0*XABS
        TITLX(3:3) = TITLNN(12)
        if(XABS.LT.1.d0) then
          XABS = 10.d0*XABS
          TITLX(2:2) = TITLNN(12)
          if(XABS.LT.1.d0) then
            XABS = 10.d0*XABS
            TITLX(1:1) = TITLNN(12)
          endif
        endif
        I0001 = int(XABS + 0.50d0)
        I0001 = min(I0001,9)
        TITLX(4:4) = TITLNN(I0001+1)
       endif
      endif

      return
      end


C---------------------------------------------------------------------
      subroutine DBLDBL (YGRD,YDGRD,XGRD,XDGRD,YEDG,YDEDG,XEDG,XDEDG,
     &            YGRDW,YDGRDW,XGRDW,XDGRDW,YEDGW,YDEDGW,XEDGW,XDEDGW,
     &            ZDEGI,ZDEGJ,IMAP,JMAP,
     &            IPARW,JPARW,IPAR,JPAR,IDGRD,JDGRD,LDEG)
C---------------------------------------------------------------------
c---Define lat-lon of ctm grid and weightings if horizontal resolution degraded
c---   for example, take 1x1 met data and run at 2x2 or 2x3
C---------------------------------------------------------------------
C     IDGRD    Number of high-res longitude grid boxes in each low-res box
c     JDGRD    Number of high-res latitude ...
C     IMAP     Index of high-res longitude grid box for each IMAPN
C     ZDEGI    Fraction of low-res box accounted for by high-res box
C---------------------------------------------------------------------
      implicit NONE

      integer, Intent(in) ::  IPARW,JPARW,IPAR,JPAR,IDGRD,JDGRD
      real*8, Intent(in) ::   YGRDW(JPARW), YDGRDW(JPARW)
     &                       ,YEDGW(JPARW+1), YDEDGW(JPARW+1)
     &                       ,XGRDW(IPARW), XDGRDW(IPARW)
     &                       ,XEDGW(IPARW+1), XDEDGW(IPARW+1)
      real*8, Intent(out) ::  YGRD(JPAR), YDGRD(JPAR)
     &                       ,YEDG(JPAR+1), YDEDG(JPAR+1)
     &                       ,XGRD(IPAR), XDGRD(IPAR)
     &                       ,XEDG(IPAR+1), XDEDG(IPAR+1)
      real*8, Intent(out) ::  ZDEGI(IDGRD,IPAR),ZDEGJ(JDGRD,JPAR)
      integer, Intent(out) :: IMAP(IDGRD,IPAR),JMAP(JDGRD,JPAR)
      logical, Intent(out) :: LDEG

      real*8, parameter ::  CPI    = 3.141592653589793d0
      real*8, parameter ::  ZPI180 = 180.d0/CPI

      integer I, J, K, JW
      real*8  Y1P1W(361), Y1P1(181), ZDY

      if (IDGRD.eq.1 .and. JDGRD.eq.1 )  then
c--- For native resolution no degradation required
        LDEG = .false.
        ZDEGI(:,:) = 1.d0
        ZDEGJ(:,:) = 1.d0

        do I = 1,IPAR
          IMAP(1,I) = I
          XGRD(I)   = XGRDW(I)
          XDGRD(I)  = XDGRDW(I)
          XEDG(I)   = XEDGW(I)
          XDEDG(I)  = XDEDGW(I)
        enddo
          XEDG(IPAR+1)  = XEDGW(IPAR+1)
          XDEDG(IPAR+1) = XDEDGW(IPAR+1)

        do J = 1,JPAR
          JMAP(1,J) = J
          YGRD(J)   = YGRDW(J)
          YDGRD(J)  = YDGRDW(J)
          YEDG(J)   = YEDGW(J)
          YDEDG(J)  = YDEDGW(J)
        enddo
          YEDG(JPAR+1)  = YEDGW(JPAR+1)
          YDEDG(JPAR+1) = YDEDGW(JPAR+1)

      else
c  Otherwise, degrade
        LDEG = .true.
          write(6,'(a,2i5)') 'Degrade Longitude resolution:', IPARW,IPAR
          write(6,'(a,2i5)') 'Degrade Latitude  resolution:', JPARW,JPAR

        do I = 1,IPAR
          do K = 1,IDGRD
            IMAP(K,I)  = K + (I-1)*IDGRD
            ZDEGI(K,I) = 1.d0 / dble(IDGRD)
          enddo
        enddo
        do J = 1,JPAR
          do K = 1,JDGRD
            JMAP(K,J)  = K + (J-1)*JDGRD
          enddo
        enddo

        do I = 1,IPAR
          XEDG(I)  = XEDGW(IMAP(1,I))
          XDEDG(I) = XDEDGW(IMAP(1,I))
        enddo
          XEDG(IPAR+1)  = XEDGW(IPARW+1)
          XDEDG(IPAR+1) = XDEDGW(IPARW+1)
        do I = 1,IPAR
          XGRD(I)  = 0.5d0*(XEDG(I)+XEDG(I+1))
          XDGRD(I) = ZPI180*XGRD(I)
        enddo

        do J = 1,JPAR 
          YEDG(J)  = YEDGW(JMAP(1,J))
          YDEDG(J) = YDEDGW(JMAP(1,J))
        enddo
          YEDG(JPAR+1)  = YEDGW(JPARW+1)
          YDEDG(JPAR+1) = YDEDGW(JPARW+1)
        do J = 1,JPAR
          YGRD(J)  = 0.5d0*(YEDG(J)+YEDG(J+1))
          YDGRD(J) = ZPI180*YGRD(J)
        enddo

          Y1P1(1) = -1.d0
        do J = 2,JPAR
          Y1P1(J) = sin(YEDG(J))
        enddo
          Y1P1(JPAR+1) = 1.d0
          Y1P1W(1) = -1.d0
        do J = 2,JPARW 
          Y1P1W(J) = sin(YEDGW(J))
        enddo
          Y1P1W(JPARW+1) = 1.d0

        do J = 1,JPAR
            ZDY = 1.d0 / (Y1P1(J+1)-Y1P1(J))
          do K = 1,JDGRD
            JW = JMAP(K,J)
            ZDEGJ(K,J)  = (Y1P1W(JW+1)-Y1P1W(JW)) * ZDY
          enddo
        enddo

c---Output mapping
          write(6,'(a)') 'Longitude re-mapping: '
        do I = 1,3
          write(6,1000) I,IDGRD,(IMAP(K,I),ZDEGI(K,I),K=1,IDGRD)
 1000     format(2i3,6(i4,' (',f6.4,')'))
        enddo
        do I = IPAR-1,IPAR
          write(6,1000) I,IDGRD,(IMAP(K,I),ZDEGI(K,I),K=1,IDGRD)
        enddo
        write(6,'(a)') 'Latitude  re-mapping: '
        do J = 1,3
          write(6,1000) J,JDGRD,(JMAP(K,J),ZDEGJ(K,J),K=1,JDGRD)
        enddo
        do J = JPAR-1,JPAR
          write(6,1000) J,JDGRD,(JMAP(K,J),ZDEGJ(K,J),K=1,JDGRD)
        enddo
      endif

      return
      end


c---------------------------------------------------------------------
      subroutine GAUSST(NG,XP,WT)
c---------------------------------------------------------------------
c---calculate NG Gaussian quadrature points (XP) & weights (WT)
c---     for interval (X1,X2)  from GISS-Lacis code.
c---tested against EC version GAUAW, simpler, same to 1.d-13 or better)

      implicit none
      integer, intent(in) ::  NG
      real*8, dimension(NG), intent(out) :: XP, WT

      real*8, parameter :: PI = 3.141592653589793d0
      real*8, parameter :: PS = 1.013211836423378d-01
      real*8, parameter :: DXL = 1.d-16
      real*8, parameter :: X1 = -1.d0,  X2 = 1.d0

      real*8 XMID,XHAF,DNG,DN,DM,DI,C,Z,ZZ,PN,PNI,PNJ,PNK,DX,X
      integer NN,N2,N,I,J

        XMID = (X2+X1)/2.D0
        XHAF = (X2-X1)/2.D0
        DNG = NG
        NN = NG/2
        N2 = NN*2
      if (N2.eq.NG) goto 110

        XP(NN+1) = XMID
        WT(NN+1) = 1.d0
      if (NG.lt.2) return

        PN = 1.d0
        N = 0
  100 N = N+2
        DN = N
        DM = DN-1.d0
        PN = PN*(DM/DN)
      if (N.lt.N2) goto 100

        WT(NN+1) = 2.d0*XHAF/(DNG*PN)**2
  110 I = 0
        C = PI/sqrt(DNG*(DNG+1.d0)+0.5d0-PS)/105.d0
  120 I = I+1
        DI = I
        Z = PS/(4.d0*DI-1.d0)**2
        ZZ = (105.d0+Z*(210.d0-Z*(2170.d0-Z*(105812.d0-12554474.d0*Z))))
        X = DCOS(ZZ*C*(DI-0.25d0))
  130 N = 1
        DM = 1.d0
        PNI = 1.d0
        PNJ = X
  140 N = N+1
        DN = N
        PNK = ((DM+DN)*X*PNJ-DM*PNI)/DN
        PNI = PNJ
        PNJ = PNK
        DM = DN
      if (N.lt.NG) goto 140

        DX = PNJ*(1.d0-X*X)/DNG/(PNI-X*PNJ)
        X = X-DX
      if (abs(DX).gt.DXL) goto 130

        J = NG+1-I
        XP(I) = XMID-XHAF*X
        XP(J) = XMID+XHAF*X
        WT(I) = 2.d0*XHAF*(1.d0-X*X)/(DNG*PNI)**2
        WT(J) = WT(I)
      if (I.lt.NN) goto 120

      return
      end



c----------------------------------------------------------------------
      subroutine LEGGEN (PLEG,PLAT,NT,NMMAX)
c----------------------------------------------------------------------
c   generate Legendre functions, rewritten from SPLEG1 (ECMWF model) 
c       NT = CTM truncation
c       NMMAX = max number of spectral coeffs = (T+1)*(T+4)/2
c       PLAT  = Latitude in radians
c
c   Adapted to CTM:    J. K. Sundet   July 1994
c----------------------------------------------------------------------
      implicit none
      integer,intent(in)      ::  NT,NMMAX
      real*8, intent(in)      ::  PLAT
      real*8, dimension(NMMAX), intent(out) :: PLEG

c  for T159, NMMAX=(159+1)*(159+4)/2 = 13040, for T319, = 51680
      real*8, dimension(51680) :: ZHLP1, ZHLP2, ZHLP3
      real*8 ZSIN,ZCOS,ZF1M,ZRE1,ZF2M,ZM,ZN,ZE1,ZE2
      integer NT1,J1M,JLM,JM,JCN,NPC
     
      if (NMMAX.gt.51680) then
        call EXITC('>>>>LEGGEN not enough scratch space T>319')
      endif

        NT1 = NT + 1
        ZSIN   = sin(PLAT)
        ZCOS   = sqrt(1.d0 - ZSIN*ZSIN)
        JLM    = 2
        PLEG(1)= 1.d0
        ZF1M   = sqrt(3.d0)
        PLEG(2)= ZF1M*ZSIN
        NPC    = 2
      do JM = 2,NT1
        ZM       = JM - 1.d0
        ZHLP1(JM)= sqrt(2.d0*ZM + 3.d0)
        ZHLP2(JM)= 1.d0/sqrt(2.d0*ZM)
        NPC      = NPC + 1
      enddo
        ZHLP1(1) = sqrt(3.d0)

      do JM = 1,NT1
        J1M = JM - 1
        ZM  = J1M
        ZRE1= ZHLP1(JM)
        ZE1 = 1.d0/ZRE1
       if (J1M.ne.0) then
          ZF2M = ZF1M*ZCOS*ZHLP2(JM)
          ZF1M = ZF2M*ZRE1
          JLM  = JLM + 1
          PLEG(JLM)= ZF2M
          JLM  = JLM + 1
          PLEG(JLM)= ZF1M*ZSIN
          NPC  = NPC + 1
       endif
       do JCN = J1M+2,NT1
         ZN   = JCN
         ZHLP3(JCN)= sqrt((4.d0*ZN*ZN - 1.d0)/(ZN*ZN - ZM*ZM))
         NPC  = NPC + 1
       enddo
       do JCN = J1M+2,NT1
         ZE2  = ZHLP3(JCN)
         JLM  = JLM + 1
         PLEG(JLM)= ZE2*( ZSIN*PLEG(JLM-1) - ZE1*PLEG(JLM-2) )
         ZE2  = 1.d0/ZE2
         ZE1  = ZE2
         NPC  = NPC + 1
       enddo
      enddo

      return
      end


      
c------------------------------------------------------------------------
      subroutine FFT_SET (TRIGS,IFAX,N)
c------------------------------------------------------------------------
c---computes factors of n & trigonometric functions required by FFT991
c
c    Adapted to CTM, March 1996: J.K. Sundet
C------------------------------------------------------------------------
      implicit none
      integer,intent(in)                  :: N
      integer,dimension(10),intent(out)   :: IFAX
      real*8, dimension(N), intent(out)   :: TRIGS

      real*8   ANGLE,DEL,FN
      integer  JFAX(10),  I,J,K,NHL,NU,  IFAC,NFAX,II1,II2
c---LFAX are the factorization for FFTs, notes say the order is 8,6,5,4,3,2?
      integer, parameter, dimension(7) :: LFAX = [6,8,5,4,3,2,1]

      real*8, parameter ::  CPI = 3.141592653589793d0

C-----------------------------------------------------------------------
        FN  = dble(N)
        NHL = int(FN*0.5d0 - 1.d0)
        DEL = 2.d0*CPI/FN
      do K=0,NHL
        ANGLE        = dble(K)*DEL
        TRIGS(2*K+1) = cos(ANGLE)
        TRIGS(2*K+2) = sin(ANGLE)
      enddo

c---this section of factorization is a mess, only updated exit calls
C Find factors of N (8,6,5,4,3,2; only one 8 allowed)
         NU = N
         K  = 0
         I  = 1
      do 100 II1=1,999
        if (NU.GT.1) then
          IFAC = LFAX(I)
          if (mod(NU,IFAC).eq.0) then
            K = K + 1
            JFAX(K) = IFAC
            if (IFAC.ne.8) then
C Apart from 8, there might be other factors that appears more that once
              do II2=1,999
                NU = NU/IFAC
                if (mod(NU,IFAC).eq.0) then
                  K = K + 1
                  JFAX(K) = IFAC
                else
                  goto 100
                endif
              enddo
            else
C Eight is encountered
              if(K.GT.1) then
                JFAX(1) = 8
                JFAX(K) = 6
              endif
              NU = NU/IFAC
            endif
          endif
C There are only six factors available
          I = I + 1
          if (I.gt.7) then
          call EXITIJL('>>>>FFT_SET: N has illegal factors',N,N,N)
          endif
          goto 100
        elseif (NU.eq.1) then
C Reverse order of factors
          NFAX    = K
          IFAX(1) = NFAX
          do J=1,NFAX
            IFAX(NFAX+2-J) = JFAX(J)
          enddo
          IFAX(10) = N
          goto 99
C NFAC.LT.1
        else
          call EXITIJL('>>>>FFT_SET: N has illegal factors',N,N,N)
        endif

  100 continue
   99 continue

      return
      end


c-----------------------------------------------------------------------
      subroutine FFT_DDSS (NT,A0,NMMAX,DD,SS)
c-----------------------------------------------------------------------
c Calculate the constant arrays needed for evaluating spectral coeffisients 
c     for u and v when the voritcity and divergence is given.
c
c   DD(N,M) = -A0 * sqrt( (N*N-M*M)/(4*N*N-1) ) / N  "divergence" factors  
c   SS(N,M) = -A0 * M / ( N*(N+1) )                  "vorticity" factors
c
c  Revised:  J. K. Sundet, January 1994 (Method from ECMWF manual)
c----------------------------------------------------------------------
      implicit none
      integer, intent(in)                  :: NT, NMMAX
      real*8, intent(in)                   :: A0
      real*8, dimension(NMMAX),intent(out) :: DD, SS

      real*8  ZM,ZN
      integer NTP1, IMD,IMN, JM,JN
c-----------------------------------------------------------------------
       NTP1 = NT+1
       IMD  = 1
       IMN  = 1
      do JM=1,NTP1
          ZM = JM - 1.d0
        do JN=JM,NTP1
           ZN = JN - 1.d0
          if (JN.GT.1) then
            DD(IMD) = -sqrt((ZN*ZN - ZM*ZM)/(4.d0*ZN*ZN - 1.d0))/ZN*A0
            SS(IMN) = -ZM*A0/(ZN*(ZN + 1.d0))
          else
            DD(IMD) = 0.d0
            SS(IMN) = 0.d0
          endif
            IMN = IMN + 1
            IMD = IMD + 1
        enddo
          ZN  = ZN + 1.d0
          DD(IMD) = -sqrt((ZN*ZN - ZM*ZM)/(4.d0*ZN*ZN - 1.d0))/ZN*A0
          IMD = IMD + 1
      enddo
      return
      end


c-----------------------------------------------------------------------
      subroutine CALENDR (IYEAR,IDAY, LLPYR,LFIXMET,MYEAR,
     &         JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)
c-----------------------------------------------------------------------
c--- compute the year/month/day counting from day# IDAY of year IYEAR
c--- note that IYEAR is the reference year and IDAY >>365 for long runs
c---     JYEAR = current year
c---     JDAY = day of year (1:365 or 366)
c---     LLPYR = allow for leap year & LFIXMET = recycle met fields(no leap yr!)
c---     JMON = no. of month (1:12)
c---     TMON = 3-char label of month
c---     JDATE = day of the month (1:31)
c---     LYEAR = .true. = current year is being treated as a leap year
c---     TMET = 3-char label for day of year (Feb 29 = '901')
c---     SOLDEC = solar declination (degrees)
c---          SINDEC = sin(solar declination)
c---     SOLDIS = distance to sun (A.U.)
c-----------------------------------------------------------------------
      implicit none
      integer, intent(in)      :: IYEAR,IDAY
      logical, intent(in)      :: LLPYR,LFIXMET
      integer, intent(out)     :: JYEAR,JDAY,JMON,JDATE,MYEAR
      logical, intent(out)     :: LYEAR
      real*8,  intent(out)     :: SOLDEC,SOLDIS
      character*3, intent(out) :: TMON,TMET

      integer KDAY,KYEAR, JMET,I
      real*8 DFYR, G,SOLLNG,SINDEC

      integer, parameter, dimension(365) ::  JMOFD = [(1,I=1,31), 
     &  (2,I=32,59),     (3,I=60,90),    (4,I=91,120), (5,I=121,151), 
     &  (6,I=152,181),  (7,I=182,212),  (8,I=213,243), (9,I=244,273), 
     & (10,I=274,304), (11,I=305,334), (12,I=335,365) ]

      integer, parameter, dimension(366) ::  JMOFDL= [(1,I=1,31), 
     &  (2,I=32,60),     (3,I=61,91),    (4,I=92,121), (5,I=122,152), 
     &  (6,I=153,182),  (7,I=183,213),  (8,I=214,244), (9,I=245,274), 
     & (10,I=275,305), (11,I=306,335), (12,I=336,366) ]

      integer, parameter, dimension(13)  :: JDOFM =
     & [0,31,59,90,120,151,181,212,243,273,304,334,365]
      integer, parameter, dimension(13)  :: JDOFML =
     & [0,31,60,91,121,152,182,213,244,274,305,335,366]

      character*3, parameter, dimension(12) :: AMON = ['JAN','FEB',
     & 'MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

c--------------sun-earth data from US Naval Obs-------------------------
c---small shifts from year-to-year, just use year 2000
c---  e(deg) = 23.439 - 0.00013 * (Y - 2000)  (obliquity)
c---Earth-sun distance astronomical units (AU)
c---  g(deg) = 356.543 + 0.98560028*JDAY  (deg, JDAY=day of year=1=Jan 1.5)
c---  R = 1.00014 - 0.01671 cos g - 0.00014 cos 2g
c---  q(deg) = 279.473 + 0.98564736*JDAY
c---  L = q + 1.915 sin g + 0.020 sin 2g  (L= Sun's apparent eclipt long)
c--- sin d = sin e * sin L (solar declination)
      real*8, parameter :: CPI180 = 0.01745329252d0
      real*8, parameter :: COBLIQ = 0.3977725d0
      real*8, parameter :: ZPI180 = 1.d0/CPI180

c---leap years are simply assumed to be divisible by 4 (ie, 1904 - 2096)
        KDAY = IDAY
        KYEAR = IYEAR
   10 continue
        LYEAR = (mod(KYEAR,4).eq.0) .and. LLPYR .and. .not.LFIXMET
      if (LYEAR) then
        if (KDAY.le.366) goto 12
        KDAY = KDAY - 366
      else
        if (KDAY.le.365) goto 12
        KDAY = KDAY - 365
      endif
        KYEAR = KYEAR + 1
      goto 10
   12 continue
c---have reached current year
         JYEAR = KYEAR
         JDAY = KDAY
      if (LFIXMET)  then
         MYEAR = IYEAR
      else
         MYEAR = JYEAR
      endif
      if (LYEAR) then
         JMON = JMOFDL(JDAY)
         JDATE = JDAY - JDOFML(JMON)
         DFYR = dble(JDAY)/366.
        if (JDAY.gt.60) then
         JMET = JDAY - 1
        elseif (JDAY.eq.60) then
         JMET = 901
        else
         JMET = JDAY
        endif
      else
         JMON = JMOFD(JDAY)
         JDATE = JDAY - JDOFM(JMON)
         DFYR = dble(JDAY)/365.
         JMET = JDAY
      endif
      TMON = AMON(JMON)

      write (TMET(1:3),'(I3.3)')  JMET

c---solar declination and distance to sun: from US Naval Obs
        G = CPI180*(356.543d0 + 0.98560028d0*dble(JDAY))
      SOLDIS = 1.00014d0 - 0.01671d0*cos(G) - 0.00014d0*cos(G+G)
        SOLLNG = 279.473d0 + 0.98564736d0*dble(JDAY) 
     &            + 1.915d0*sin(G) + 0.02d0*sin(G+G)
        SINDEC = COBLIQ * sin(SOLLNG*CPI180)
      SOLDEC = asin(SINDEC)*ZPI180
c
      return
      end
      

c----------------------------------------------------------------------
      subroutine EXITC (MESSAG)
c-----------------------------------------------------------------------
      implicit none
      character*(*),intent(in)::   MESSAG
        write(6,'(a)')  MESSAG
      stop
      end

      
c----------------------------------------------------------------------
      subroutine EXITL (MESSAG,LABEL)
c-----------------------------------------------------------------------
      implicit none
      character*(*),intent(in)::   MESSAG, LABEL
        write(6,'(2a)')  MESSAG,LABEL
      stop
      end


c----------------------------------------------------------------------
      subroutine EXITIJL (MESSAG,I,J,L)
c-----------------------------------------------------------------------
      implicit none
      character*(*),intent(in)::   MESSAG
      integer, intent(in):: I,J,L
        write(6,'(a)')  MESSAG
        write(6,'(3i10)') I,J,L
      stop
      end
