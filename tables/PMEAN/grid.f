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
      include 'cmn.f'
      real*8,  intent(in):: GM0000
      integer, intent(in):: JMPOLR

      integer I,J,L,LW, IMB2,JMB2, IDTLN
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
      
C-----areas and air mass
        AREAG = 0.d0
      do J=1,JM
      do I=1,IM
        AREAG = AREAG + AREAXY(I,J)
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
