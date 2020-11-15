C****************************************************************************
C                Quenching Weights for Heavy Quarks 
C                (Multiple Soft Scattering Approximation)
C        	             January 26, 2005
C
C       Refs:
C    N.Armesto, A.Dainese, C.A.Salgado and U.A.Wiedemann, hep-ph/0501225
C    N.Armesto, C.A.Salgado and U.A.Wiedemann, hep-ph/0312106
C 
C
C   This package contains quenching weights for gluon radiation off 
C   massive quarks in the multiple soft scattering approximation.
C
C   qwmassmult(mmmm,rrrr,xxxx,continuous,discrete) 
C   returns the quenching weight for a quark with M/E = mmmm 
C   traversing a medium with transport coeficient q and length L. 
C   The input values are rrrr=0.5*q*L^3 and xxxx=DeltaE/wc, where
C   wc=0.5*q*L^2 and DeltaE is the energy radiated. The output values are
C   the continuous and discrete (prefactor of the delta function) parts
C   of the quenching weights.
C
C   In order to use this routine, the files cont_mass_mult.all and 
C   disc_mass_mult.all need to be in the working directory. 
C
C   An initialization of the tables is done by calling initmassmult when
C   qwmassmult() is called for the first time.
C
C   Please, send us any comment:
C
C       nestor.armesto@cern.ch
C       andrea.dainese@pd.infn.it
C       carlos.salgado@cern.ch 
C       urs.wiedemann@cern.ch 
C 
C-----------------------------------------------------------------------------

      subroutine qwmassmult(mmmm,rrrr,xxxx,continuous,discrete)
C
C     Compute weights for massive quark (mmmm=M/E), rrrr=0.5*q*L^3
C     and xxxx=w/wc, wc=0.5*q*L^2
C
      real*8           xxx(416), mmm(8), rrr(10), dw(8,10), cw(8,10,416)
      common /data/    xxx, mmm, rrr, dw, cw
*
      real*8           rrrr, xxxx, mmmm, continuous, discrete
      real*8           rrin, xxin, mmin
      integer          nrlow, nrhigh, nxlow, nxhigh, nmlow, nmhigh
      real*8           rrhigh, rrlow, rfraclow, rfrachigh
      real*8           xxhigh, xxlow, xfraclow, xfrachigh
      real*8           mmhigh, mmlow, mfraclow, mfrachigh
      real*8           crlowmlow, crlowmhigh, crhighmlow, crhighmhigh
      real*8           cmlow, cmhigh
      real*8           dmlow, dmhigh
*
      mmin = mmmm
      rrin = rrrr
      xxin = xxxx

*
C     find the two values of M/E to be used for the interpolation
*
      if (mmin.le.mmm(8)) mmin = 1.05d0*mmm(8);
      if (mmin.ge.mmm(1)) mmin = 0.95d0*mmm(1);
*
c      print*, 'mmm(8) = ', mmm(8)
c      print*, 'mmm(1) = ', mmm(1)
c      print*, 'mmin = ', mmin

      do 666,nm=1,8
         if (mmin.lt.mmm(nm)) then
            mmhigh = mmm(nm)
         else
            mmhigh = mmm(nm-1)
            mmlow  = mmm(nm)
            nmlow  = nm
            nmhigh = nm-1
            goto 665
         endif
 666  enddo
 665  continue
      mmin = mmmm

      mfraclow  = (mmhigh-mmin)/(mmhigh-mmlow)
      mfrachigh = (mmin-mmlow)/(mmhigh-mmlow)
*
C     find the two values of R to be used for the interpolation
*
      if (rrin.le.rrr(10)) rrin = 1.05d0*rrr(10)
      if (rrin.ge.rrr(1))  rrin = 0.95d0*rrr(1)


      do 777,nr=1,10
         if (rrin.lt.rrr(nr)) then
            rrhigh = rrr(nr)
         else
            rrhigh = rrr(nr-1)
            rrlow  = rrr(nr)
            nrlow  = nr
            nrhigh = nr-1
            goto 776
         endif
 777  enddo
 776  continue
      rrin = rrrr

      rfraclow  = (rrhigh-rrin)/(rrhigh-rrlow)
      rfrachigh = (rrin-rrlow)/(rrhigh-rrlow)
*
C     find the two values of x to be used for the interpolation
*
      if (xxin.le.xxx(1))   xxin = 1.05d0*xxx(1)
      if (xxin.ge.xxx(416)) xxin = 0.95d0*xxx(416)

      do 888,nx=1,416
         if (xxin.gt.xxx(nx)) then
            xxlow = xxx(nx)
         else
            xxlow = xxx(nx-1)
            xxhigh  = xxx(nx)
            nxhigh  = nx
            nxlow = nx-1
            goto 887
         endif
 888  enddo
 887  continue
      xxin = xxxx

      xfraclow  = (xxhigh-xxin)/(xxhigh-xxlow)
      xfrachigh = (xxin-xxlow)/(xxhigh-xxlow)

*
C     calculate interpolation for continuous and discrete weight
*
      crlowmlow   = xfraclow*cw(nmlow,nrlow,nxlow)
     &             +xfrachigh*cw(nmlow,nrlow,nxhigh)
      crlowmhigh  = xfraclow*cw(nmhigh,nrlow,nxlow)
     &             +xfrachigh*cw(nmhigh,nrlow,nxhigh)
      crhighmlow  = xfraclow*cw(nmlow,nrhigh,nxlow)
     &             +xfrachigh*cw(nmlow,nrhigh,nxhigh)
      crhighmhigh = xfraclow*cw(nmhigh,nrhigh,nxlow)
     &             +xfrachigh*cw(nmhigh,nrhigh,nxhigh)

      cmlow = rfraclow*crlowmlow + rfrachigh*crhighmlow
      cmhigh = rfraclow*crlowmhigh + rfrachigh*crhighmhigh

      continuous = mfraclow*cmlow + mfrachigh*cmhigh

      dmlow  = rfraclow*dw(nmlow,nrlow) 
     &        +rfrachigh*dw(nmlow,nrhigh)
      dmhigh = rfraclow*dw(nmhigh,nrlow) 
     &        +rfrachigh*dw(nmhigh,nrhigh)

      discrete = mfraclow*dmlow + mfrachigh*dmhigh
      if (discrete.lt.0.d0) discrete = 0.d0

c      print*, 'cont = ', continious
c      print*, 'disc = ', discrete


      end
C-------------------------------------------------------------------------
      subroutine initmassmult
*     
C     Read-in data
C     mmm(8)   = 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001
C     rrr(10)  = 1M, 400k, 100k, 40k, 10k, 4k, 1k, 200, 50, 10 
C     xxx(416) = 0 -> 20
*
      real*8           xxx(416), mmm(8), rrr(10), dw(8,10), cw(8,10,416)
      common /data/    xxx, mmm, rrr, dw, cw
      
      Character*255 ENVDIR
      Common /ENVCOM/ ENVDIR
      Character*255 Filename1
      Character*255 Filename2
*
      Filename1=Trim(ENVDIR)//'/PyQM/qweight/cont_mass_mult.all'
      open(unit=20,file=Filename1,status='OLD',err=90)
      do 109,nm=1,8 
         do 110,nx=1,416
            read (20,*) xxx(nx),
     &           cw(nm,1,nx), cw(nm,2,nx), cw(nm,3,nx), 
     &           cw(nm,4,nx), cw(nm,5,nx), cw(nm,6,nx), 
     &           cw(nm,7,nx), cw(nm,8,nx), cw(nm,9,nx),
     &           cw(nm,10,nx)
110     continue
 109  continue
      close(20)
*
      Filename2=Trim(ENVDIR)//'/PyQM/qweight/disc_mass_mult.all'
      open(unit=22,file=Filename2,status='OLD',err=91)
      do 111,nm=1,8
         do 112,nr=1,10
            read (22,*) mmm(nm), rrr(nr), dw(nm,nr)
 112     continue
 111  continue
      close(22)
*
      goto 888
 90   print*, 'input - output error' 
 91   print*, 'input - output error #2' 
 888  continue

      end

