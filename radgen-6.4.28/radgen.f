
C **********************************************************************
C
C     RADATA
C
C **********************************************************************

      BLOCK DATA RADATA
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "ppicom.inc"

      data
     +     amm/2.7928456d0/,
     +     amn/-1.913148d0/,
     +     chbar/.197328d0/,
     +     barn/.389379d6/,
     +     aml/.511000d-3/,
     +     aml2/.261112d-6/,
     +     al2/.522240d-6/,
     +     isf20/8/,
     +     pi/3.1415926d0/,
     +     pi2/9.869604d0/,
     +     alfa/.729735d-2/,
     +     amc2/1.151857d0/,
     +     amp/.938272d0/,
     +     amh/.938272d0/,
     +     i2/1,1,1,2,3,3,1,2/,
     +     i1/3,3,4,4,3,3,3,3/

      end

C **********************************************************************
C
C     RADGEN
C
C **********************************************************************

      SUBROUTINE RADGEN(e1,q2l,nul,ys,xs,phil,PhRAD,q2tr,anutr,WEIGHT)
      IMPLICIT NONE
***********************************************************
*  INPUT :
*        e1    : beam energy
*        q2,nu : uncorrected kinematics of scattered lepton
*        phil  : azimuthal angle of scattered lepton
*  OUTPUT :
*          PhRAD(4) : real rad. photon 4 vector
*          q2tr     : true Q2
*          Utr      : true U
***********************************************************

      include "cmpcom.inc"
      include "radgen.inc"
      include "phiout.inc"
      include "tailcom.inc"
      include "radgenkeys.inc"

      REAL PhRAD(*),q2tr,anutr
      REAL e1,q2l,nul,phil,xs,ys,weight
      INTEGER itagen
*

* reset proton mass (abr)
* necessary to get correct kinematics after first elastic event !
      amp=0.938272d0

C... Don't recalculate x and y, but use the ones already calculated
C      ys=geny
C      xs=genx
C      ys=nul/e1
C      xs=q2l/(2.*amp*nul)

C      write(*,*)"radgen: ",ys, q2l, xs, nul, phil, e1
C Do we want to use a lookup table?
      if(ixytest.ge.0)then
        if(kill_elas_res.ne.2)then
          call itafun(ys,xs,plrun,pnrun,ixytest,itagen)
        else
          itagen=0
          sigcor=1.
        endif
      else
	itagen=-1
      endif

      if(itagen.ne.0.and.abs(ixytest).ne.2)then
	iphi=0
	call mpolrad(e1,ys,xs,1.,plrun,pnrun,itagen)
	if(ixytest.eq.1)ixytest=-2
      endif
      call radgam_pol(e1,ys,xs,phil,ixytest,itagen,q2tr,anutr)
      
* set kinematics of possibly radiated real gamma
      phrad(1)=sngl(dplabg(1))
      phrad(2)=sngl(dplabg(2))
      phrad(3)=sngl(dplabg(3))
C to avoid rounding errors !
CCC      phrad(4)=ys*e1-anutr
      phrad(4)=sqrt(phrad(1)*phrad(1)+phrad(2)*phrad(2)+
     +              phrad(3)*phrad(3))

* set the weight
      if(kill_elas_res.eq.2)then
        weight=1.
      else
        weight = sigcor
      endif

      return
      end

C **********************************************************************
C
C     PHIDST_POL
C
C **********************************************************************

      subroutine phidst_pol
*     -----------------
*
*     calculate phi-distribution of bremsstrahlung photons
*     needed for correction of hadronic distributions
*
      include "double.inc"
      include "phiout.inc"
      include "cmpcom.inc"
      include "radgen.inc"
      include "ppicom.inc"
C
      dimension db(3),dx(31),dy(31),dz(31)
      external dfufi_pol
      data ndim /31/
C
C     boundaries for 2 integrations
C
      db(1)=0.d0
      db(2)=pi/10.d0
      db(3)=pi
      k=0
      dsum=0.d0
      do i=1,2
        dbnd1=db(i)
        dbnd2=db(i+1)
        call radquad_pol(dfufi_pol,dbnd1,dbnd2,dx,dy,dz,dh,ndim)
        do j=1,ndim
             if(i.eq.1.or.j.ne.1) k=k+1
             dphi(k)=dx(j)
             ddephi(k)=dh
             dsumph(k)=dsum+dz(j)
        enddo
        dsum=dsum+dz(ndim)
      enddo
      kmp=k
      return
      end

C **********************************************************************
C
C     RADGAM_POL
C
C **********************************************************************

      subroutine radgam_pol(erad,yrad,xrad,genphi,
     +     ixytest,itagen,q2tr,anutr)

      include "double.inc"
      include "phiout.inc"
      include "radgen.inc"
      include "cmpcom.inc"
      include "sxycom.inc"
      include "tailcom.inc"
      include "amf2com.inc"
      include "density.inc"

      real*8 podinl
      real temp             ! Temp variable for a denom

C     selection of mass dmj and mass index iadmj
      r1=rlu(0)
      if(ixytest.lt.0)then
        rr1=r1*(tbor+tine+tpro+tnuc)
C      print *,tine,tpro,tnuc
        if (rr1.gt.(tine+tpro+tnuc)) then
          ita=0
        elseif(rr1.gt.(tpro+tnuc)) then
          ita=1
          sicurr=rr1-tpro-tnuc
          call conk2(dble(erad),dble(xrad),dble(yrad),1)
        elseif(rr1.gt.tnuc)then
          ita=3
          sicurr=rr1-tnuc
          call conk2(dble(erad),dble(xrad),dble(yrad),1)
        else
          ita=2
          sicurr=rr1
          call conk2(dble(erad),dble(xrad),dble(yrad),2)
        endif
      else
        xs=dble(xrad)
        ys=dble(yrad)
        ita=itagen
        if(ita.eq.1)sicurr=r1*tine
        if(ita.eq.2)sicurr=r1*tnuc
        if(ita.eq.3)sicurr=r1*tpro
        if(ita.eq.2)then
          call conk2(dble(erad),dble(xrad),dble(yrad),2)
        else
          call conk2(dble(erad),dble(xrad),dble(yrad),1)
C           print *,an,xs,ys,ita,s,x
        endif
        isf1=1
        isf2=isf20
        isf3=1

C      print *,an,xs,ys,erad,xrad,yrad
C      sicurr=scgen
      endif

      if(ita.ne.0)then
*
*---     selection of theta angle
*
* maximum value of dsitkm can be smaller than 1 !
* ---> use modified random number (as in old generator)
* abr (03.01.02)
*       print *,dsitkm(ndxtkm(ita),ita)
        r1=r1*dsitkm(ndxtkm(ita),ita)
        do ktk=1,ndxtkm(ita)
C   print *,dsitkm(ktk,ita),r1,ndxtkm(ita)
          if (r1.le.dsitkm(ktk,ita)) then
            if (ktk.eq.1) then
              dtk1=r1/dsitkm(ktk,ita)
            else
              dtk1=(r1-dsitkm(ktk-1,ita))/
     +             (dsitkm(ktk,ita)-dsitkm(ktk-1,ita))
            endif
            dtk=dcvtkm(ktk,ita)+(dtk1-0.5d0)*ddetkm(ktk,ita)
            goto 30
          endif
        enddo
        dtk=dcvtkm(ndxtkm(ita),ita)
 30     continue
* up to this point identical to old generator
* exception: calculation of integrand
*
        taout=(sx-sqly*cos(dtk))/ap2
C       print *,taout,dtk,ita,ktk,r1
C       print *,sx,sqly,ap2
        if(ita.eq.1)then
          taa=taout
          iphi=0
          call tails(taout,tm)
          rcal=ap*demin
          rmax=(w2-amc2)/(1.+taout)
          do krr=1,nrr
            ddeler(krr,ktk)=(rmax-rcal)/nrr
            drcurr(krr,ktk)=rcal+ddeler(krr,ktk)*(krr-0.5)
C             print *,krr,ktk,drcurr(krr,ktk)
            if(krr.eq.1)then
             dsigmr(krr,ktk)=podinl(drcurr(krr,ktk))
            else
             dsigmr(krr,ktk)=dsigmr(krr-1,ktk)+podinl(drcurr(krr,ktk))
            endif
C           print *,dsigmr(krr,ktk),drcurr(krr,ktk),taout
          enddo
          r2=rlu(0)
          dsircu=r2*dsigmr(nrr,ktk)
          do krr=1,nrr
            if (dsircu.le.dsigmr(krr,ktk)) then
              if(krr.eq.1)then
                drr1=dsircu/dsigmr(krr,ktk)
              else
                drr1=(dsircu-dsigmr(krr-1,ktk))
     +               /(dsigmr(krr,ktk)-dsigmr(krr-1,ktk))
              endif
              drr=drcurr(krr,ktk)+(drr1-0.5d0)*ddeler(krr,ktk)
              goto 20
            endif
          enddo
          drr=drcurr(nrr,ktk)
 20       continue
          dom=drr/ap
        else
          dom=(sx-y)/ap/(1.+taout)
        endif
        rrout=ap*dom

* the following is again the same as in the old generator (abr 03.01.02)
*
*---selection of phi angle
*
*
C  print *,rrout,taout
C  pause
        iphi=1
        call phidst_pol
*
        r3=rlu(0)
        dphpoi=r3*dsumph(kmp)
        do kph=2,kmp
          if (dphpoi.le.dsumph(kph)) then
* avoid division by zero (abr)
*           dphk1=(dphpoi-dsumph(kph-1))/
*    +           (dsumph(kph)-dsumph(kph-1))
            temp = dsumph(kph) - dsumph(kph-1)
            if (temp .eq. 0) temp = 1e-20
            dphk1 = (dphpoi - dsumph(kph-1)) / temp
            dphk=dphi(kph)+(dphk1-1.0d0)*ddephi(kph)
            goto 40
          endif
        enddo
        dphk=dphi(kmp)
 40     continue
        r4=rlu(0)
        if (r4.gt.0.5) dphk=-dphk
        
*
*     radiative photon
*
        deg=dom
        dthg=dtk
        dphig=dphk
        sigma_total=sngl(tbor+tine+tpro+tnuc)
        
        q2tr=y+rrout*taout
        anutr=sx/ap-dom
C       write(6,*) 'ita,deg,dthg,dphig,q2tr,anutr,xtr',
C     &              ita,deg,dthg,dphig,q2tr,anutr,
C     &              q2tr/(2.*amp*anutr)
        
        dgpz=deg*dcos(dthg)
        dgpxy=deg*dsin(dthg)
        dgpx=dgpxy*dcos(dphig)
        dgpy=dgpxy*dsin(dphig)
*     
*---momentum components in the LAB-system:
*---two rotations needed - first within the scattering plane
*---around the y-axis and second around the new z-axis (beam
*---direction) by Phi of the scattering plane
*
        dgplx=-dgpz*dsts+dgpx*dcts
        dgply=dgpy
        dgplz=dgpz*dcts+dgpx*dsts
        dcphi=dcos(dble(genphi))
        dsphi=dsin(dble(genphi))
        dplabg(1)=dgplx*dcphi-dgply*dsphi
        dplabg(2)=dgplx*dsphi+dgply*dcphi
        dplabg(3)=dgplz
        
      else

        q2tr=2d0*amh*erad*xrad*yrad
        anutr=yrad*erad
C         write(*,'(i5,5f8.3)')ita,xs,ys,q2tr,anutr,sigcor

        dplabg(1)=0.0d0
        dplabg(2)=0.0d0
        dplabg(3)=0.0d0

      endif

* reset kinematics (abr)
*     call conk2(dble(erad),dble(xrad),dble(yrad),1)
* now done at beginning of sr radgen

      end

C **********************************************************************
C
C     DFUFI_POL
C
C **********************************************************************

      double precision function dfufi_pol(dx)
*     -----------------------------------
*     needed for correction of hadronic distributions
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "sxycom.inc"
      include "tailcom.inc"
      include "amf2com.inc"
      include "ppicom.inc"
      include "radgen.inc"

      phipoi=dx
      taa=taout
      call tails(taout,tm)
      dfufi_pol=an*alfa/pi*podinl(rrout)
      
      return
      end
      
C **********************************************************************
C
C     RADQUAD_POL
C
C **********************************************************************

      subroutine radquad_pol(dfunct,dlow,dup,dx,dy,dz,dh,ndim)
*     -------------------------------------------------
*
      include "double.inc"
      dimension dx(31),dy(31),dz(31)
      external dfunct
*
      dsum2=0.d0
      if (ndim.gt.1) then
        dh=(dup-dlow)/float(ndim-1)
        do i=1,ndim
          daux=dlow+dh*float(i-1)
          dx(i)=daux
          dy(i)=dfunct(daux)
        enddo
        do i=2,ndim
          dsum1=dsum2
          dsum2=dsum2+.5d0*(dx(i)-dx(i-1))*(dy(i)+dy(i-1))
          dz(i-1)=dsum1
        enddo
        dz(ndim)=dsum2
      elseif (ndim.eq.1) then
        dz(ndim)=dsum2
      endif
*     
      return
      end

C **********************************************************************
C
C     MPOLRAD
C
C **********************************************************************

      subroutine mpolrad(e1curr,yscurr,xscurr
     +     ,uncurr,plcurr,pncurr,itagen)
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "sxycom.inc"
      include "ppicom.inc"
      include "tailcom.inc"
      include "radgen.inc"
      include "radgenkeys.inc"
      include "deltacom.inc"
      include "density.inc"
      include "pypars.inc"

      real e1curr,yscurr,xscurr,uncurr,plcurr,pncurr

      dimension tai(5),si(2,3),si0(2,3),tls(2,3,4)

      e1=e1curr
      xs=xscurr
      ys=yscurr
      pl=plcurr
      pn=-pncurr
      un=uncurr

      if (abs(pncurr).eq.2.) then
c        negative polarization is twice as large as positive one
         qn=pncurr*1./2.
         if(qn.lt.0) qn=2.*qn
      else
         qn=0.
      endif

      call conk2(e1,xs,ys,1)
c
c
c delta is factorizing part of virtual and real leptonic bremsstrahlung
c
      call deltas(delta,delinf,tr)

      il=1
      in=1

      si(il,in)=0d0
      si0(il,in)=0d0
      tls(il,in,1)=0d0
      tls(il,in,2)=0d0
      tls(il,in,3)=0d0
      tls(il,in,4)=0d0
      
      isf1=1
      isf2=isf20
      isf3=1

      do ii=1,4
        tai(ii)=0d0
      enddo
      sib=0d0
      sia=0d0

      if(itagen.eq.-1)then
       ita1=1
       ita4=4
      else
       ita1=itagen
       ita4=itagen
      endif

      do 30 ita=ita1,ita4
c
c     sib is born cross section with polarized initial
c     lepton and proton
c     sia is contribution of anomalous magnetic moment.
c
* for photoproduction set elastic and quasielastic tails to zero
* abr (19.02.04)
      if (ntx.eq.ntpho .and. (ita.eq.2.or.ita.eq.3)) then
        tai(2)=0.0d0
        tai(3)=0.0d0
        goto 30
      endif

C ECA to use Antjes trick to kill the elastic and quasielastic tails
C always if we run pythia with radiative corrections
c      if (MSTP(199).eq.1) then
c        tai(2)=0.0d0
c        tai(3)=0.0d0
c        goto 30
c      endif

        if(ita.eq.2.and.kill_elas_res.eq.1)goto30
        if(ita.eq.3.and.ire.le.1)then
          tai(3)=0d0
C           write(9,'('' tai = 0.0  '')')
          goto 30
        end if
        if(ita.eq.1)then
          call bornin(sib,sia)
        endif
c
c     tai(1),tai(2),tai(3) are contributions of radiative tails:
c     1 - inelastic
c     2 - elastic
c     3 - quasielastic
c
        if(ita.eq.2) call conk2(e1,xs,ys,ita)

        if(kill_elas_res.ne.2) call qqt(sib,tai(ita))

* what is this ??? (abr 17.10.01)
* shouldnot this if statement be not equal here ???
* abr   if(ita.eq.2) call conk2(e1,xs,ys,1)
*       if(ita.ne.2) call conk2(e1,xs,ys,1)
* no, purpose is to reset kinematics, I think (abr 03.01.02)
* thus the original code is correct
        if(ita.eq.2) call conk2(e1,xs,ys,1)

 30   continue
C       extai1=exp(alfa/pi*delinf)
C       extai2=((sx-y/rtara)**2/s/(s-y/rtara))**tr
C       extai3=((sx-y)**2/s/(s-y))**tr

      extai1=1.
      extai2=1.
      extai3=1.
      delinf=0.

C       sig=sib*extai1*(1.+alfa/pi*(delta-delinf))+sia
C      +       +tai(4)+tai(5)
cilyichev only one-loop without multy-soft photon emission
c      sig=sib*redfac*(1.+vertex+vac+small)+sia
      sig=sib*(1.+log(redfac)+vertex+vac+small)+sia
     +     +tai(4)
     +     +tai(1)+(tai(2)*extai2+tai(3)*extai3)/rtara

      si(il,in)=si(il,in)+sig
      si0(il,in)=si0(il,in)+sib
      tls(il,in,1)=tls(il,in,1)+tai(1)
      tls(il,in,2)=tls(il,in,2)+tai(2)*extai2/rtara
      tls(il,in,3)=tls(il,in,3)+tai(3)*extai3/rtara
C       tls(il,in,4)=sib*extai1*(1.+alfa/pi*(delta-delinf))+sia
C      +     +tai(4)+tai(5)
      tls(il,in,4)=sib*redfac*(1.+vertex+vac+small)+sia
     +     +tai(4)

C      write(*,'(1x,8g11.4)')xs,ys,si0(il,in),si(il,in)
C      +  ,tls(il,in,1),tls(il,in,2),tls(il,in,3)
C      +  ,tls(il,in,4)
      if(itagen.eq.-1.or.itagen.eq.1)then
C        if(abs(tine-tls(il,in,1))/tine.gt.5d-2)
C      +     write(*,*)' tine',tine,tls(il,in,1),itagen
        tine=tls(il,in,1)
      endif
      if(itagen.eq.-1.or.itagen.eq.2)then
C        if(abs(tnuc-tls(il,in,2))/tnuc.gt.5d-2)
C      +     write(*,*)' tnuc',tnuc,tls(il,in,2),itagen
        tnuc=tls(il,in,2)
      endif
      if(itagen.eq.-1.or.itagen.eq.3)then
C      if(abs(tpro-tls(il,in,3))/tpro.gt.5d-2)
C      +     write(*,*)' tpro',tnuc,tls(il,in,3),itagen
        tpro=tls(il,in,3)
      endif
      if(itagen.eq.-1)then
        tbor=tls(il,in,4)
        sigrad=si(il,in)
        sig1g=si0(il,in)
      endif
C #ifdef TEST
C       write(*,'(1x,6g11.4)')sig1g,sigrad,tine,tnuc,tpro,tbor
C #endif
      sigcor=sigrad/sig1g
      if(kill_elas_res.eq.2)sigcor=1.
      
C #ifdef TEST
C       print *,' kill_elas_res=',kill_elas_res,amc2
C       print *,' sig1g=',sig1g,sib
C       print *,' sigrad=',sigrad,sig
C       print *,' tine =',tine
C       print *,' tai(4) =',tai(4)
C       print *,' tnuc =',tnuc
C       print *,' tpro =',tpro
C       print *,' tbor =',tbor
C       print *,' sig1g=',sig1g
C       print *,' sigcor=',sigcor
C       print *,' delta=',alfa/pi*delta
C       print *,' vac  =',vac
C       print *,' vertex=',vertex
C       print *,' small=',small
C       print *,' tai(4)= ',sib*log(redfac)
C       sigt1=sib*alfa/pi*delta+tai(5)
C       sigt2=sib*(log(redfac)+vertex+vac+small)
C       sigt3=sib*alfa/pi*(delta+delta5)
C       print *,' sigt=',sigt1,sigt2,sigt3
C #endif
*     print *,' kill_elas_res=',kill_elas_res,amc2
*     print *,' sig1g=',sig1g,sib
*     print *,' sigrad=',sigrad,sig
*     print *,' tine =',tine
*     print *,' tai(4) =',tai(4)
*     print *,' tnuc =',tnuc
*     print *,' tpro =',tpro
*     print *,' tbor =',tbor
*     print *,' sig1g=',sig1g
*     print *,' sigcor=',sigcor
*     print *,' delta=',alfa/pi*delta
*     print *,' vac  =',vac
*     print *,' vertex=',vertex
*     print *,' small=',small
*     print *,' redfac= ',redfac

      end

C **********************************************************************
C
C     CONK2
C
C **********************************************************************

      subroutine conk2(e1,xs,ys,iittaa)
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"

      if(iittaa.eq.2)then
        amp=amt
      else
        amp=amh
      endif
C      print *,amp,amh,iittaa
      call conkin(e1,xs*amh/amp,ys)

      return
      end

C **********************************************************************
C
C     CONKIN
C
C **********************************************************************

      subroutine conkin(e1,xss,yss)

c set of kinematical constants
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "polcom.inc"
      include "sxycom.inc"
      include "ppicom.inc"
      include "radgen.inc"
cilyichev additional include
      include "mc_set.inc"

      ap=2.*amp
      amp2=amp**2
      ap2=2.*amp**2
      s=ap*e1
      x=s*(1-yss)
      sx=s-x
      sxp=s+x
      y=s*xss*yss
      ym=y+al2
      tpl=s**2+x**2
      tmi=s**2-x**2
      w2=amp2+s-y-x
      als=s*s-al2*ap2
      alx=x*x-al2*ap2
      alm=y*y+4.*aml2*y
      aly=sx**2+4.*amp2*y
      sqls=dsqrt(als)
      sqlx=dsqrt(alx)
      sqly=dsqrt(aly)
      sqlm=dsqrt(alm)
      allm=dlog((sqlm+y)/(sqlm-y))/sqlm

      coe=xs/e1/1d3

      axy=pi*(s-x) * coe
      an=2.*alfa**2/sqls*axy*barn*amh/amp
      tamin=(sx-sqly)/ap2
      tamax=(sx+sqly)/ap2

      dcts=(s*sx + ap2*y)/sqly/sqls
* avoid values above 1.0 (for x values below 1e-07, abr 18.2.04)
*     dsts=sin( acos(dcts) )
      dsts=sin( acos(MIN(MAX(dcts,-1.0),+1.0)) )

      as=s/2./aml/sqls
      bs=0.
      cs=-aml/sqls
cilyichev instead of
c        ae=amp/sqls
c        be=0.
c        ce=-s/ap/sqls
c we need
      if ((mcSet_PTarget(1:1).eq.'L').or.
     +    (mcSet_PTarget(1:2).eq.'DT')) then
         ae=amp/sqls
         be=0.
         ce=-s/ap/sqls
      elseif(mcSet_PTarget(1:1).eq.'T')then
         sqn=dsqrt(s*x*y-aly*aml2-amp2*y*y)
         ae=(-s*x+ap2*ym)/sqls/sqn/2.
         be=sqls/sqn/2.
         ce=-(s*y+al2*sx)/sqls/sqn/2.
      endif

      apq=-y*(ae-be)+ce*sx
      apn=(y+4.*aml2)*(ae+be)+ce*sxp
      dk2ks=as*ym+al2*bs+cs*x
      dksp1=as*s+bs*x+cs*ap2
      dapks=2.*(al2*(as*ae+bs*be)+ap2*cs*ce+ym*(as*be+bs*ae)
     +     +s*(as*ce+cs*ae)+x*(bs*ce+cs*be))

      return
      end

C **********************************************************************
C
C     BORNIN
C
C **********************************************************************

      subroutine bornin(sibor,siamm)
c
c     sibor is born cross section with polarized initial
c     lepton and polarized target
c     siamm is contribution of anomalous magnetic moment.
C     the cross section is calculated as dsigma/dlogQ2dnu to fit
C     to the HERMES-disng generation

c
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "polcom.inc"
      include "sxycom.inc"
      include "ppicom.inc"
      include "tailcom.inc"

      dimension sfm0(8),tm(8)

      call strf(0d0,0d0,sfm0)
      tm(1)=-(2.*aml2-y)
      tm(2)=(-(amp2*y-s*x))/(2.*amp2)
      tm(3)=(2.*(apq*dk2ks-dapks*y)*aml)/amp
      tm(4)=apq/amp*(-(dk2ks*sx-2.*dksp1*y)*aml)/amp2
      tm(7)=(-(4.*aml2+3.*apn**2-3.*apq**2+y))/2.
      tm(8)=apq/amp*(-3.*(apn*sxp-apq*sx))/(2.*ap)
      ek=(3.*apq**2-y)/amp2
      tm(5)=-ek*tm(1)
      tm(6)=-ek*tm(2)
      ssum=0.
      do 1 isf=isf1,isf2,isf3
         ppol=1.
         if(isf.eq.3.or.isf.eq.4)ppol=pl*pn
C        if(isf.eq.3.or.isf.eq.4)ppol=-pn
         if(isf.ge.5)ppol=qn/6
        ssum=ssum+tm(isf)*sfm0(isf)*ppol
    1 continue
      sibor=ssum*an/y/y*2.
c     
c     formula (4) of kukhto and shumeiko paper
c
cc    res1=amp*ww1*(y+4.*aml2)-ww2*(s+x)**2/4./amp
cc    siamm=alfa/pi*al2*allm*(sibor+an*res1/y**2)
      siamm=0.

      return
      end

C **********************************************************************
C
C     DELTAS
C
C **********************************************************************

      subroutine deltas(delta,delinf,tr)
c
c delta is factorizing part of virtual and real leptonic bremsstrahlung
c
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "sxycom.inc"
      include "ppicom.inc"
      include "radgen.inc"
      include "deltacom.inc"

c
c    am2 : squared masses of charge leptons
c
      dimension am2(3)
      data am2/.26110d-6,.111637d-1,3.18301d0/

      del1=-ym*(alm*allm**2/2.+2.*fspen(2d0*sqlm/(y+sqlm))-pi2/2.)/sqlm
      del2=(3.*y/2.+4.*aml2)*allm-2.

      suml=0.
      do 10 i=1,3
        a2=2.*am2(i)
        sqlmi=dsqrt(y*y+2.*a2*y)
        allmi=dlog((sqlmi+y)/(sqlmi-y))/sqlmi
        suml=suml+2.*(y+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./y
 10   continue
      if(y.lt.1.d0)then
        aaa = -1.345d-9
        bbb = -2.302d-3
        ccc = 4.091
      elseif(y.lt.64d0)then
        aaa = -1.512d-3
        bbb =  -2.822d-3
        ccc = 1.218
      else
        aaa = -1.1344d-3
        bbb = -3.0680d-3
        ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*y)) *2*pi/alfa
      sum=suml+sumh
C       print *,' vacl_my_old=',alfa/pi*suml
C       print *,' vach_my_old=',alfa/pi*sumh
C       print *,' vac_my_old=',alfa/pi*sum
      sum=vacpol(y)
C       print *,' vac_my_new=',alfa/pi*sum

      aj0=2.*(ym*allm-1.)
c...elke okay the dlog((w2-amc2) can become negative if w2 is smaller than
C   amc2/1.151857d0/
c     write(*,*)w2,aml,amc2,aj0
      deltai=aj0*dlog((w2-amc2)/aml/dsqrt(w2))
c     write(*,*)"deltai=",deltai

      ss=x+y
      xx=s-y
      alss=ss**2-2.*w2*al2
      alxx=xx**2-2.*w2*al2
      sqlss=dsqrt(alss)
      sqlxx=dsqrt(alxx)
      allss=dlog((sqlss+ss)/(-sqlss+ss))/sqlss
      allxx=dlog((sqlxx+xx)/(-sqlxx+xx))/sqlxx
      dlm=dlog(y/aml2)
      sfpr=dlm**2/2.-dlm*dlog(ss*xx/aml2/w2)
     +     -(dlog(ss/xx))**2/2.+fspen((s*x-y*amp2)/ss/xx)-pi2/3.
      delta0=(ss*allss+xx*allxx)/2.+sfpr
      delta=deltai+delta0+del1+del2+sum
      delinf=(dlm-1.)*dlog((w2-amc2)**2/ss/xx)
      tr=alfa/pi*(dlm-1.)

      vac=alfa/pi*sum
      vertex=alfa/pi*del2
      small_old=alfa/pi*(pi2/6.-fspen(1.-amp2*y/s/x)
     +     + fspen(1.-s/x) + fspen(1.-x/s))
      small=alfa/pi*(-pi2/6.-log(ss*xx/s/x)*log(sx/y)
     +     +log(ss/s)*log(xx/x)
     +     -1.5*log(s/x)**2+2.*log(xx/x)*log(s/y)+2.*log(ss/s)*log(x/y)
     +     -2.*fspen(-y/x)-2.*fspen(y/s)+2.*fspen(y/sx)+fspen(s*x/ss/xx)
     +     -fspen(s*y/ss/sx)-fspen(x*y/sx/xx))
C       print *,' small_old=',small_old
C       print *,' small_new=',small
        redfac=exp(-alfa/pi*(dlm-1.)*log(s*x/(4.*amp2*demin**2) ))
        delta5=-(dlm-1.)*log((w2-amc2)**2*s*x/(4.*amp2*demin**2*ss*xx))
     +       -log(xx/x)*log(amp2*y/s**2)-log(ss/s)*log(amp2*y/x**2)
     +       -2.*fspen(-y/x)-2.*fspen(y/s)+2.*fspen(-tamin)
     +       +2.*fspen(-tamax)
     +       -fspen((-y-s*tamin)/xx)
     +       -fspen((-y-s*tamax)/xx)
     +       -fspen(( y-x*tamin)/ss)
     +       -fspen(( y-x*tamax)/ss)

      return
      end

C **********************************************************************
C
C     VACPOL
C
C **********************************************************************

      double precision function vacpol(y)

      implicit real*8(a-h,o-z)

      include "ppicom.inc"

c
c    am2 : squared masses of charge leptons
c
      dimension am2(3),a(5),b(5),c(5)
      data am2/.26110d-6,.111637d-1,3.18301d0/
      data a/0d0,0d0,0d0,1.2227d-3,1.64178d-3/
      data b/2.2877d-3,2.51507d-3,2.79328d-3,2.96694d-3,2.92051d-3/
      data c/4.08041425d0,3.09624477d0,2.07463133d0,1d0,1d0/

      suml=0.
      do 10 i=1,3
        a2=2.*am2(i)
        sqlmi=dsqrt(y*y+2.*a2*y)
        allmi=dlog((sqlmi+y)/(sqlmi-y))/sqlmi
        suml=suml+2.*(y+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./y
 10   continue

      if(y .lt. 4d0)then
       k=1
      elseif(y .lt. 16d0)then
       k=2
      elseif(y .lt. 100d0)then
       k=3
      elseif(y .lt. 8317.44d0)then
       k=4
      elseif(y .ge. 8317.44d0)then
       k=5
      else
       stop ' Y<0 in VACPOL'
      endif

      sumh = (a(k)+b(k)*log(1.+c(k)*y)) *2.*pi/alfa
      vacpol=suml+sumh

      end

C **********************************************************************
C
C     FSPENS
C
C **********************************************************************

      double precision function fspens(x)
c
c    spence function
c
      implicit real*8(a-h,o-z)

      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
 1    an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
 2    fspens=f

      return
      end

C **********************************************************************
C
C     FSPEN
C
C **********************************************************************

      double precision function fspen(x)

      implicit real*8(a-h,o-z)

      data f1/1.644934d0/

      if(x)8,1,1
 1    if(x-.5d0)2,2,3
 2    fspen=fspens(x)
      return
 3    if(x-1d0)4,4,5
 4    fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
 5    if(x-2d0)6,6,7
 6    fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
 7    fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
 8    if(x+1d0)10,9,9
 9    fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
 10   fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return

      end

C **********************************************************************
C
C     QQT
C
C **********************************************************************

      subroutine qqt(bo,tai)

      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "sxycom.inc"
      include "ppicom.inc"
      include "tailcom.inc"
      include "radgen.inc"

      external rv2
      dimension dbtk(8)
      data ep/1d-8/

      dsumtk=0.d0
      derrtk=0.d0
      isumtk=0
      if(ita.eq.1 .or. ita.eq.5)then
        tade=(w2-amc2)/(ap*demin) -1.d0
        costkm=(sx-ap2*tade)/sqly
        if(costkm.ge.1d0)then
          tai=0.
          return
        endif
        dtkmax=acos( max(-1d0,costkm ))
      else
        dtkmax=pi
      endif
      call intbtk2(dbtk,nbtk,dtkmax)
c     integrate each bin by ntk subbins
c     ntk=number of bins within the big bin dbtk(i)...dbtk(i+1)
      do 10 itk=1,nbtk
        call inttk2(isumtk,dbtk(itk),dbtk(itk+1),dsumtk,derrtk)
C     write(*,*)itk,isumtk,dsumtk

 10   continue
      if(ita.le.3)ndxtkm(ita)=isumtk
      tai=dsumtk

      end

C **********************************************************************
C
C     TAILS
C
C **********************************************************************

      subroutine tails(ta,tm)
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "polcom.inc"
      include "sxycom.inc"
      include "ppicom.inc"
      include "radgen.inc"
      include "bseocom.inc"

      dimension tm(8,6),ajm2(2),ajm3(3),ii(8)
      data ii/1,2,3,4,1,2,5,6/

      if(iphi.eq.0)then
        b2=(-aly*ta+sxp*sx*ta+2.*sxp*y)/2.
        b1=(-aly*ta-sxp*sx*ta-2.*sxp*y)/2.
        c1=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(s*ta+y)**2)
        c2=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(ta*x-y)**2)
        bb=1./sqly
        sc1=dsqrt(c1)
        sc2=dsqrt(c2)
        bi12=(sxp*(sx*ta+2.*y))/(sc1*sc2*(sc1+sc2))
        bi1pi2=1./sc2+1./sc1
        bis=-b1/sc1/c1+b2/sc2/c2
        bir=b2/sc2/c2+b1/sc1/c1
        b1i=-b1/aly/sqly
        b11i=(3.*b1**2-aly*c1)/2./aly**2/sqly
      else
C        write(*,*) ta,tamin,tamax,s,x,y,aml2,aly
C...elke: to avoid problems with the sqrt from rounding at the very low Q2
C         if using pythia6 ----> y becomes very small
        sqrtmb=((ta-tamin)*(tamax-ta)*(s*x*y-(y**2*amp2)-(aml2*aly)))
        if (sqrtmb.ge.0) then
            sqrtmb=sqrt(sqrtmb)
        else
           write(*,*)'radgen: ta,tamin,tamax,s,x,y,aml2,aly',
     &                 ta,tamin,tamax,s,x,y,aml2,aly
           write(*,*) "radgen: sqrtmb=",sqrtmb
           sqrtmb=0
        endif

        z1=(y*sxp+ta*(s*sx+ap2*y)-ap*cos(phipoi)*sqrtmb)/aly
        z2=(y*sxp+ta*(x*sx-ap2*y)-ap*cos(phipoi)*sqrtmb)/aly
        bb=1./sqly/pi
        bi12=bb/(z1*z2)
        bi1pi2=bb/z2+bb/z1
        bis=bb/z2**2+bb/z1**2
        bir=bb/z2**2-bb/z1**2
        b1i=bb*z1
        b11i=bb*z1**2
      endif
      sps=as+bs
      spe=ae+be
      ccpe=(ae-be)*ta+2.*ce
      ccps=(as-bs)*ta+2.*cs
      sis=(2.*bi1pi2*sps+bir*sps*ta+bis*ccps)/2.
      sir=( (2.*bi12*sps*ta+bir*ccps+bis*sps*ta))/2.
      si12=(bi12*ccps+bi1pi2*sps)/2.
      eis=(2.*bi1pi2*spe+bir*spe*ta+bis*ccpe)/2.
      eir=( (2.*bi12*spe*ta+bir*ccpe+bis*spe*ta))/2.
      ei12=(bi12*ccpe+bi1pi2*spe)/2.
      ois=((2.*bi1pi2+bir*ta)*(ccpe*sps+ccps*spe)+(ccpe*ccps+
     +     spe*sps*ta**2)*bis+8.*bb*spe*sps+4.*bi12*spe*sps*ta**2)/4.
      oir=( ((2.*bi12+bis)*(ccpe*sps+ccps*spe)*ta+(ccpe*ccps+
     +     spe*sps*ta**2)*bir+4.*bi1pi2*spe*sps*ta))/4.
      oi12=((ccpe*ccps+spe*sps*ta**2)*bi12+(ccpe*sps+ccps*spe)*
     +     bi1pi2+4.*bb*spe*sps)/4.
      eeis=((ccpe**2+spe**2*ta**2)*bis+8.*bb*spe**2+4.*bi12*spe
     +     **2*ta**2+4.*bi1pi2*ccpe*spe+2.*bir*ccpe*spe*ta)/4.
      eeir=( ((ccpe**2+spe**2*ta**2)*bir+4.*bi12*ccpe*spe*ta+4.
     +     *bi1pi2*spe**2*ta+2.*bis*ccpe*spe*ta))/4.
      eei12=((ccpe**2+spe**2*ta**2)*bi12+4.*bb*spe**2+2.*bi1pi2
     +     *ccpe*spe)/4.
      ei1pi2=(4.*bb*spe+bi12*spe*ta**2+bi1pi2*ccpe)/2.
      eei1i2=((ccpe**2+spe**2*ta**2)*bi1pi2+4.*(2.*ccpe-spe*ta)
     +     *bb*spe+8.*b1i*spe**2+2.*bi12*ccpe*spe*ta**2)/4.
      eb=((ccpe-spe*ta)*bb+2.*b1i*spe)/2.
      eeb=((ccpe-spe*ta)**2*bb+4.*(ccpe-spe*ta)*b1i*spe+4.*b11i
     +     *spe**2)/4.
      call ffu(1,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     +     ,eis,eir,ei12,ei1pi2,ta)
      call ffu(2,eb,eis,eir,ei12,ei1pi2,oir,ois,oi12
     +     ,eeis,eeir,eei12,eei1i2,ta)
      call ffu(3,eeb,eeis,eeir,eei12,eei1i2,0d0,0d0,0d0
     +     ,0d0,0d0,0d0,0d0,ta)
      ajm2(1)=apq/amp
      ajm2(2)=-1./amp
      ajm3(1)=(y-3.*apq**2)/amp2
      ajm3(2)=6.*apq/amp2
      ajm3(3)=-3./amp2
      do 15 i=1,8
        do 13 l=1,6
          tm(i,l)=0
 13     enddo
        do 10 k=1,i2(i)
          ajk=1.
          if(i.eq.4.or.i.eq.8)ajk=ajm2(k)
          if(i.eq.5.or.i.eq.6)ajk=ajm3(k)
          do 11 j=k,i1(i)+k-1
            tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,k)*ajk
            if((i.eq.5.or.i.eq.6).and.k.eq.2)
     +           tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,1)*ta/amp2
 11       continue
 10     continue
 15   continue
      
      return
      end
      
C **********************************************************************
C
C     FFU
C
C **********************************************************************

      subroutine ffu(n,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     +     ,eis,eir,ei12,ei1pi2,ta)
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "polcom.inc"
      include "sxycom.inc"
      include "ppicom.inc"
      include "bseocom.inc"

      hi2=aml2*bis-ym*bi12
      shi2=aml2*sis-ym*si12
      ehi2=aml2*eis-ym*ei12
      ohi2=aml2*ois-ym*oi12
      goto(10,20,30)n
 10   continue
      tm3(3,1,n)=(8.*(apq*dk2ks-dapks*y)*aml*hi2)/amp
      tm3(3,2,n)=(-2.*((2.*(bi12*dk2ks*ta-2.*shi2)*apq+(2.*shi2-
     +     sir*y+sis*ym)*apn+4.*dapks*hi2*ta)-4.*((2.*ei12-eis)*
     +     dk2ks-(si12-sis)*apn)*aml2)*aml)/amp
      tm3(3,3,n)=(2.*(((2.*si12+sir-sis)*apn*ta-2.*dk2ks*ei12*ta
     +     -6.*ohi2-oir*y+ois*ym)-4.*aml2*oi12)*aml)/amp
      tm3(3,4,n)=(2.*(2.*oi12-oir+ois)*aml*ta)/amp
      tm3(5,1,n)=-2.*(4.*aml2+3.*apn**2-3.*apq**2+y)*hi2
      tm3(5,2,n)=-2.*(6.*aml2*apn*eir-3.*apn**2*bi12*ta+3.*apn*
     +     apq*bi1pi2+6.*apq*ehi2+hi2*ta)
      tm3(5,3,n)=-(24.*aml2*eei12-6.*apn*ei1pi2-6.*apq*ei12*ta-
     +     2.*bb-bi12*ta**2)
 20   continue
      tm3(4,1,n)=(-4.*(dk2ks*sx-2.*dksp1*y)*aml*hi2)/amp2
      tm3(4,2,n)=(((2.*(sxp-2.*sx)*shi2+2.*bi12*dk2ks*sx*ta+8.*
     +     dksp1*hi2*ta-sir*sxp*y+sis*sxp*ym)-4.*(2.*bi12*dk2ks-bis*
     +     dk2ks-si12*sxp+sis*sxp)*aml2)*aml)/amp2
      tm3(4,3,n)=((((sxp*ta-ym)*sis-(sxp*ta-y)*sir+2.*bi12*dk2ks
     +     *ta+6.*shi2-2.*si12*sxp*ta)+4.*aml2*si12)*aml)/amp2
      tm3(4,4,n)=(-(2.*si12-sir+sis)*aml*ta)/amp2
      tm3(6,1,n)=(-3.*(apn*sxp-apq*sx)*hi2)/amp
      tm3(6,2,n)=(-3.*(2.*(apn*bir+eir*sxp)*aml2-(2.*bi12*sxp*ta
     +     -bi1pi2*sx)*apn+(bi1pi2*sxp+2.*hi2)*apq+2.*ehi2*sx))/(2.*
     +     amp)
      tm3(6,3,n)=(-3.*(8.*aml2*ei12-apn*bi1pi2-apq*bi12*ta-ei12*
     +     sx*ta-ei1pi2*sxp))/(2.*amp)
 30   continue
      tm3(1,1,n)=-4.*(2.*aml2-y)*hi2
      tm3(1,2,n)=4.*hi2*ta
      tm3(1,3,n)=-2.*(2.*bb+bi12*ta**2)
      tm3(2,1,n)=(((sxp**2-sx**2)-4.*amp2*y)*hi2)/(2.*amp2)
      tm3(2,2,n)=(2.*aml2*bir*sxp-4.*amp2*hi2*ta-bi12*sxp**2*ta+
     +     bi1pi2*sxp*sx+2.*hi2*sx)/(2.*amp2)
      tm3(2,3,n)=(2.*(2.*bb+bi12*ta**2)*amp2+4.*aml2*bi12-bi12*
     +     sx*ta-bi1pi2*sxp)/(2.*amp2)

      return
      end

C **********************************************************************
C
C     PODINL
C
C **********************************************************************

      double precision function podinl(r)
c
c     integrand (over r )
c
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "sxycom.inc"
      include "tailcom.inc"
      include "ppicom.inc"
      include "amf2com.inc"

      dimension sfm(8),iel(8)
      data iel/0,2,1,2,2,4,2,3/

      if(ita.ne.5) call strf(taa ,r,sfm)
      podinl=0.
      do 11 isf=isf1,isf2,isf3
        ppol=1.
        if(isf.eq.3.or.isf.eq.4)ppol=pl*pn
C        if(isf.eq.3.or.isf.eq.4)ppol=-pn
        if(isf.ge.5)ppol=qn/6
        if(ita.eq.2)ppol=ppol*(amt/amp)**iel(isf)
        irrend=i1(isf)+i2(isf)-1
        if(ita.eq.5)irrend=1
        do 1 irr=1,irrend
          if(ita.eq.5)then
            pres=-sfm0(isf)/2./y**2/r
          else
            pp=sfm(isf)
            if(irr.eq.1.and.ita.eq.4)pp=pp-sfm0(isf)*(1.+r*taa/y)**2
            pres=pp*r**(irr-2)/(y+r*taa)**2/2.
          endif
          podinl=podinl-tm(isf,irr)*pres*ppol
C           print *,' isf=',isf,' irr=',irr,podinl
C           write(*,*)' isf=',isf,' irr=',irr,podinl,sfm(isf)
 1      continue
 11   continue
      podinl=podinl*(1.+alfa/pi*vacpol(Y+R*taa))

      return
      end

C **********************************************************************
C
C     RV2
C
C **********************************************************************

      double precision function rv2(ta)
c
c     integrand (over ta )
c
      implicit real*8(a-h,o-z)
      external podinl

      include "cmpcom.inc"
      include "sxycom.inc"
      include "tailcom.inc"
      include "amf2com.inc"
      include "ppicom.inc"
      include "radgen.inc"

      taa=ta
      call strf(0d0,0d0,sfm0)
      call tails(ta,tm)
      rmin=1d-8
      rcal=ap*demin
      rmax=(w2-amc2)/(1.+ta)
      if(ita.eq.1)then
C         call dqn32(rmin,rmax,podinl,res)
        dsumtk=0.d0
        derrtk=0.d0
     
        dbmj2=rmax
        dbmj1=rcal
        nmj=nrr
        dr=(rmax-rcal)/nrr
        rcur=rcal-dr*.5d0
        dsum1=0.d0
c     loop over all bins
        do irr=1,nrr
          rcur=rcur+dr
          drcurr(irr,itkcur) = rcur
          ddeler(irr,itkcur) = dr
          dsum1=dsum1+podinl(rcur)
          dsigmr(irr,itkcur) = dr*dsum1
        enddo

        res=dr*dsum1

      elseif(ita.eq.2 .or. ita.eq.3)then
        aa=amt/amp
        if(ita.eq.3)aa=1.
        res=podinl((sx-y/aa)/(1d0+ta/aa))/(1.+ta/aa) /aa**2
      elseif(ita.eq.4)then
        rend=min(rcal,rmax)
        call dqn32(rmin,rend,podinl,res)
C       call qunc8(podinl,rmin,rend,1d-4,1d-9,res,er,nn,fl,3500)
      elseif(ita.eq.5)then
        res=podinl(1.d0/log(rmax/rcal))
      endif
      rv2=res
C       write(9,*)res,dr,rcur

      return
      end

C **********************************************************************
C
C     STRF
C
C **********************************************************************

      subroutine strf(ta,rr,sfm)
c
c     the programm calculates deep inelastic (ita=1),
c     elastic (ita=2), quasielastic (ita=3) structure functions
c     in kinematical point (ta,rr).
c     rr=sx-tt,
c     ta=(t-y)/rr,
c     where tt=t+amf2-amp2, amf2 is invarint mass of final hadrons
c
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "sxycom.inc"
      include "tailcom.inc"
      include "radgenkeys.inc"
      include "pypars.inc"

      real forgetp,forintp
      external forgetp,forintp
      dimension sfm(8)

      t=y+rr*ta
      if(ita.eq.1 .or. ita.eq.4 .or. ita.eq.5)then

* DIS case
        tt=sx-rr
        amf2=tt-t+amp2
        aks=t/tt
        anu=tt/ap
        epsi=ap2/tt

* get unpolarized cross sections
        itara=int(rtara+0.1)
        itarz=int(rtarz+0.1)
C....elke the order to call first mkr and than mkf2 is important for pythia6
C         to have R in mkf2 availabale
        call mkr(t,aks,dr)
C        if(aks.gt.0.99D0) then
C          write(*,*)"sx = ",sx , "y =", y, "ta = ", ta, "w2 =", w2
C          write(*,*)"q2 =",t, "xbj =",aks, "rr =", rr, "ita =", ita
C        endif
        call mkf2(t,aks,itara,itarz,f2,ff1)
        if (MSTP(199).eq.0) then
          f1=amp*(1.+anu**2/t)*f2/anu/(1.+dr)
        else
          f1=ff1
        endif

C        call mkf2(t,aks,itara,itarz,f2,ff1)
C        call mkr(t,aks,dr)
C        f1=amp*(1.+anu**2/t)*f2/anu/(1.+dr)

* get polarized cross sections
      if ((plrun.ne.0).and.(pnrun.ne.0)) then
         call mkasym (t,aks,itara,itarz,asym1,asym2)
c ilyichev g1=f1*asym1*fdilut(t,aks,itara)

cilyichev definition g1 and g2 via asym1, asym2 and f1.
         ga=sqrt(t)/anu
         g1=(asym1+ga*asym2)/(1+ga**2)*f1*fdilut(t,aks,itara)
         g2=(asym2-ga*asym1)/ga/(1+ga**2)*f1*fdilut(t,aks,itara)
         b1=0.d0
         b2=0.d0
         b3=0.d0
         b4=0.d0
      elseif((plrun.eq.0).and.(abs(pnrun).eq.2)) then
         call mkasym(t,aks,itara,itarz,asym1,asym2)
         b1=-3./2.*f1*asym1
         b2=2.d0*aks*b1
         b3=0.d0
         b4=0.d0
C        do we need
C       g1=0.d0
C       g2=0.d0
      else
        b1=0.d0
        b2=0.d0
        b3=0.d0
        b4=0.d0
        g1=0.d0
        g2=0.d0
      endif
       goto 10
      endif
c
c   rtarn,rtarz,rtara are n,z,a of the nucleus
c
      epsi=ap2/t
      rtarn=rtara-rtarz
c
c     tf is t in fermi**(-2)
c
      tf=t/chbar**2
c
c   gep,gmp are electric and magnetic form factors of proton
c   s.i.bilenkaya et al pisma zhetf 19(1974)613
c
      call ffpro(t,gep,gmp)
* add neutron electric and magnetic form factors
C-----fit of g.hohler et al.,nucl.phys. b114(1976)505.
C         FIT 8.2
c
      f1s= 0.71/(0.783**2+t)-0.64/(1.02**2+t)-0.13/(1.80**2+t)
      f2s=-0.11/(0.783**2+t)+0.13/(1.02**2+t)-0.02/(1.80**2+t)
c
      f1v= 0.05/(1.210**2+t)-0.52/(2.45**2+t)+0.28/(2.95**2+t)
      f2v=-1.99/(1.210**2+t)+0.20/(2.45**2+t)+0.19/(2.95**2+t)
c
      f1rho=(0.955+0.090/(1+t/0.355)**2)/(1+t/0.536)/2.
      f2rho=(5.335+0.962/(1+t/0.268)   )/(1+t/0.603)/2.
c
      f1v=f1v+f1rho
      f2v=f2v+f2rho
c
      f1p=f1s+f1v
      f1n=f1s-f1v
c
      f2p=f2s+f2v
      f2n=f2s-f2v
c
* keep original proton form factors and use only neutron from this fit
*     gep=f1p-t*f2p/4./amp2
*     gmp=f1p+f2p
      gen=f1n-t*f2n/4./amp2
      gmn=f1n+f2n

      if(ita.eq.2)then
        tau=t/4./amt**2
        tau1=1.+tau
        if(ire.eq.0) then
          f1=2.*amp2*tau*gmn**2
          f2=4.*amp2*tau*(gen**2+tau*gmn**2)/tau1
          g1=amp2*tau*2.*gmn*(gen+tau*gmn)/tau1
          g2=2.*amp2*tau**2*gmn*(gen-gmn)/tau1
          b1=0d0
          b2=0d0
          b3=0d0
          b4=0d0
        elseif(ire.eq.1)then
          f1=2.*amp2*tau*gmp**2
          f2=4.*amp2*tau*(gep**2+tau*gmp**2)/tau1
          g1=amp2*tau*2.*gmp*(gep+tau*gmp)/tau1
          g2=2.*amp2*tau**2*gmp*(gep-gmp)/tau1
          b1=0d0
          b2=0d0
          b3=0d0
          b4=0d0
        elseif (ire.eq.2)then
          call ffdeu(t,fcdeu,fmdeu,fqdeu)
          fc=fcdeu*amt
          fm=fmdeu*amt
          fq=fqdeu*amt
          f1=4./3.*tau*tau1*fm**2
          f2=4./9.*tau*(8.*tau**2*fq**2+6.*tau*fm**2+9.*fc**2)
          g1=-1./3.*tau*fm*(2.*tau*fq+6.*fc+3.*tau*fm)
          g2=-1./3.*tau**2*fm*(2.*tau*fq+6*fc-3.*fm)
          b1=2.*tau**2*fm**2
          b2=4.*fm**2*tau**2+
     +         16./3.*tau**3/tau1*(tau*fq+3.*fc-3.*fm)*fq
          b3=-4./3.*(3.*tau+2.)*fm**2*tau**2+
     +         16./9.*tau**3/tau1*(tau*fq+3.*fc-3.*fm)*fq
          b4=4./3.*(6.*tau+1.)*fm**2*tau**2+
     +         16./9.*tau**3/tau1*(tau*fq+3.*fc+3.*(3.*tau+2.)*fm)*fq
        elseif (ire.eq.3)then
          call ffhe3(t,ge,gm)
          f1=2.*amt**2*tau*gm**2
          f2=4.*amt**2*tau*(ge**2+tau*gm**2)/tau1
          g1=amt**2*tau*2.*gm*(ge+tau*gm)/tau1
          g2=2.*amt**2*tau**2*gm*(ge-gm)/tau1
          b1=0d0
          b2=0d0
          b3=0d0
          b4=0d0
        elseif(ire.eq.4)then
          call ffhe4(t,ge,gm)
*Test     call ffhe3(t,ge,gm)
          f1=0d0
          f2=4.*amp2*tau*(rtarz*ge)**2
          g1=0d0
          g2=0d0
          b1=0d0
          b2=0d0
          b3=0d0
          b4=0d0
        elseif((ire.eq.14).or.(ire.eq.84).or.(ire.eq.20))then
          ff=forgetp(sngl(t))
c         if (t.lt.0.002) then 
c           ff=forgetp(sngl(t))
c           print *,t,forgetp(sngl(t))
c           stop
c         endif
* copy from Polrad - abr
          f1=0d0
          f2=4.*amp2*tau*(rtarz*ff)**2
          g1=0d0
          g2=0d0
          b1=0d0
          b2=0d0
          b3=0d0
          b4=0d0
        endif
      elseif (ita.eq.3) then
        tau=t/4./amp**2
        tau1=1.+tau
        call ffquas(t,geun,gmun,gepo,gmpo)
        f1=2.*amp2*tau*gmun**2
        f2=4.*amp2*tau*(geun**2+tau*gmun**2)/tau1
        g1=amp2*tau*2.*gmpo*(gepo+tau*gmpo)/tau1
        g2=2.*amp2*tau**2*gmpo*(gepo-gmpo)/tau1
        b1=0.
        b2=0.
        b3=0.
        b4=0
      endif
 10   continue
      
      sfm(1)=f1+qn/6.*b1
      sfm(2)=epsi*(f2+qn/6.*b2)
      sfm(3)=epsi*(g1+g2)
      sfm(4)=epsi**2*g2
      sfm(5)=epsi**2*b1
      sfm(6)=epsi**3*(b2/3.+b3+b4)
      sfm(7)=epsi*(b2/3.-b3)
      sfm(8)=epsi**2*(b2/3.-b4)

      return
      end

C **********************************************************************
C
C     FFPRO
C
C **********************************************************************

      subroutine ffpro(t,gep,gmp)
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"

      gep=1.2742/(1.+t/0.6394**2)-.2742/(1.+t/1.582**2)
      gmp=(1.3262/(1.+t/0.6397**2)-.3262/(1.+t/1.3137**2))*amm
c     gep=1./((1.+.61*t)*(1.+2.31*t)*(1.+.04*t))
c     gmp=amm*gep

      end

C **********************************************************************
C
C     FFDEU
C
C **********************************************************************

      subroutine ffdeu(t,gc,gm,gq)
      implicit real*8(a-h,o-z)

      parameter (c2i3 = 2/3)
      dimension a(4),b(4),c(4)
      dimension al2ar(4),be2ar(4),ga2ar(4)
      data amd/1.8756280D0/chbar/0.19732858d0/
      data amp/0.9382796D0/
      data dmu/0.857406d0/dqu/25.84d0/
      integer para

      para=2

C      gd2=1d0/(1.d0+t/4./0.8952**2)**4
      gd2=1d0/(1.d0+t/4./0.89852**2)**4
      eta=t/4.d0/amd**2
      gd2e=gd2/(2d0*eta+1d0)
      sq2e=sqrt(2d0*eta)

      if(para.eq.1) then                                          
        al2ar(1)=1.8591*chbar**2                                 
        al2ar(4)=2d0*amd*0.58327d0                               
        be2ar(1)=19.586*chbar**2                                 
        be2ar(4)=2d0*amd*0.1d0                                   
        ga2ar(1)=1.0203*chbar**2                                 
        ga2ar(4)=2d0*amd*0.17338d0                               
      else                                                       
        al2ar(1)=1.5250*chbar**2                                 
        al2ar(4)=2d0*amd*0.24086d0                               
        be2ar(1)=43.677*chbar**2                                 
        be2ar(4)=2d0*amd*0.029138d0                              
        ga2ar(1)=1.8706*chbar**2                                 
        ga2ar(4)=2d0*amd*0.42693d0                               
      endif

      do i=2,3
        al2ar(i)=al2ar(1)+(al2ar(4)-al2ar(1))/3d0*(i-1)
        be2ar(i)=be2ar(1)+(be2ar(4)-be2ar(1))/3d0*(i-1)
        ga2ar(i)=ga2ar(1)+(ga2ar(4)-ga2ar(1))/3d0*(i-1)
      enddo

      if(para.eq.1) then                                          
        a(1)=2.4818*chbar**2                                     
        a(2)=-10.850*chbar**2                                    
        a(3)=6.4416*chbar**2                                     
      else                                                       
        a(1)=1.5706*chbar**2                                     
        a(2)=12.238*chbar**2
        a(3)=-42.046*chbar**2                                    
      endif
      a(4)=al2ar(4)*(1d0-a(2)/al2ar(2)-a(3)/al2ar(3)-a(1)/al2ar(1))

      if(para.eq.1) then                                          
        b(1)=-1.7654*chbar                                       
        b(2)=6.7874*chbar                                        
      else                                                       
        b(1)=0.0704*chbar                                        
        b(2)=0.1444*chbar                                        
      endif 

      bzn=1d0/be2ar(4)-1d0/be2ar(3)
      bbb=(2d0-dmu*amd/amp)/2d0/sqrt(2d0)/amd
      b(3)=(b(1)/be2ar(1)+b(2)/be2ar(2)-b(1)/be2ar(4)-b(2)/be2ar(4)
     +     -bbb)/bzn
      b(4)=-(b(1)/be2ar(1)+b(2)/be2ar(2)-b(1)/be2ar(3)-b(2)/be2ar(3)
     +     -bbb)/bzn
      ccc=(1d0-dmu*amd/amp-dqu)/4./amd**2

      if(para.eq.1) then                                         
        c(1)=-0.053830d0                                         
      else                                                       
        c(1)=-0.165770d0                                         
      endif

      znc2=ga2ar(1)*(ga2ar(3)*ga2ar(4)-ga2ar(2)*ga2ar(3)
     +     +ga2ar(2)**2-ga2ar(2)*ga2ar(4))
      c(2)=-ga2ar(2)/znc2*(c(1)*(
     +     ga2ar(3)*ga2ar(4)-ga2ar(1)*ga2ar(3)+ga2ar(1)**2
     +     -ga2ar(1)*ga2ar(4))-ccc*ga2ar(3)*ga2ar(4)*ga2ar(1) )
      znc3=ga2ar(1)*(ga2ar(3)-ga2ar(2))*(ga2ar(4)-ga2ar(3))
      c(3)=ga2ar(3)/znc3*(c(1)*(
     +     ga2ar(2)*ga2ar(4)-ga2ar(1)*ga2ar(4)+ga2ar(1)**2
     +     -ga2ar(1)*ga2ar(2))-ccc*ga2ar(2)*ga2ar(4)*ga2ar(1) )
      znc4=ga2ar(1)*(ga2ar(4)-ga2ar(2))*(ga2ar(4)-ga2ar(3))
      c(4)=-ga2ar(4)/znc4*(c(1)*(
     +ga2ar(2)*ga2ar(3)-ga2ar(1)*ga2ar(3)+ga2ar(1)**2-ga2ar(1)*ga2ar(2)
     +     )-ccc*ga2ar(2)*ga2ar(3)*ga2ar(1) )

      g00=0d0
      gp0=0d0
      gpm=0d0
      sqt=sqrt(t)
      do i=1,4
       g00=g00+a(i)/(al2ar(i)+t)
       gp0=gp0+sqt*b(i)/(be2ar(i)+t)
       gpm=gpm+t*c(i)/(ga2ar(i)+t)
      enddo

      gc=gd2e*( (1d0-c2i3*eta)*g00+4d0*c2i3*sq2e*gp0
     +     +c2i3*(2d0*eta-1d0)*gpm)
      gm=gd2e*(2d0*g00+2d0*(2d0*eta-1d0)/sq2e*gp0-2d0*gpm)
      gq=gd2e*(-g00+2d0/sq2e*gp0-(1d0+1d0/eta)*gpm)

      end

C **********************************************************************
C
C     FFHE3
C
C **********************************************************************

      subroutine ffhe3(t,ge,gm)
c
c   l.i.shiff phys. rev. 133(1964)3b,802
c
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"

      tf=t/chbar**2
      qf=sqrt(tf)
      a=.675
      b=.366
      c=.836
      am=.654
      bm=.456
      cm=.821
      d=-6.78d-3
      p=.9
      q0=3.98
      f0=ddexp(-a**2*qf**2) - b**2*qf**2*ddexp(-c**2*qf**2)
      fm=ddexp(-am**2*qf**2) - bm**2*qf**2*ddexp(-cm**2*qf**2)
      df=d*ddexp(-((qf-q0)/p)**2)
      ge=f0+df
      gm=fm*rtara/rtarz * (-2.13)

      end
C
      subroutine ffhe4(dq2,ge,gm)

C charge form factor by SOG(Sum-of-Gaussian)
C Helium-4,Carbon-12 :    I.Sick, Phys.Lett.116B(1982)212.
C For the Method see also I.Sick, Nucl.Phys.A218(1974)509.
C                      or M.Arneodo, thesis Princeton 1991 page 270
C Data:  H. De Vries et al., Atomic Data and Nuclear Data Tables, 36 (1987) 495
C Q2 IN 1/FERMI**2
C DGAM2=(RP/SQRT(3/2))**2 = 0.66666666            for He4
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"

      DIMENSION DR(20),DQ(20)
      DATA DR/.2D0,.6D0,.9D0,1.4D0,1.9D0,2.3D0,2.6D0
     :     ,3.1D0,3.5D0,4.2D0,4.9D0,9*0.D0/
      DATA DQ/.034724D0,.430761D0,.203166D0,.192986D0,.083866D0
     : ,.033007D0,.014201D0,0.D0,.00686D0,0.D0,.000438D0,9*0.D0/
      DATA DGAM2/.666666666667D0/
      DATA NSOG /11/

*     write(6,*) 'q2=',dq2

      DQ2F=DQ2/CHBAR**2
      DQ1F=DSQRT(DQ2F)
      DELAST=0.D0
      DO I=1,NSOG
        D2RG=2.D0*DR(I)**2/DGAM2
        DQRI=DQ1F*DR(I)
        DJ0=1.D0
        IF(DQRI.GE.1.D-6) DJ0=DSIN(DQRI)/DQRI
        DAI=DQ(I)/(1.D0+D2RG)*(DCOS(DQRI)+D2RG*DJ0)
        DELAST=DELAST+DAI
      enddo
      DELAST=DEXP(-DQ2F*DGAM2/4.D0) * DABS(DELAST)
      ge=delast
      gm=0.0

      end


      REAL FUNCTION FORGETp(Q2)
C     GET FORM FACTOR BY INTERPOLATION OUT OF TABLE FFNUC

      implicit none
      include "tailcom.inc"
      integer iq
      real q2

      FORGETp=0.
C     LOWER BIN OF Q2
      IQ=MAX0(INT(Q2/Q2BIN)+1,1)
      IF(IQ+1.LT.NQBIN) FORGETp=FFNUC(IQ)+
     1     (FFNUC(IQ+1)-FFNUC(IQ))*(Q2/Q2BIN-FLOAT(IQ-1))
c     if (iq.lt.20) print *,iq,ffnuc(iq),forgetp
c     print *,iq,q2,ffnuc(iq),forgetp
      RETURN
      END

      SUBROUTINE FORDOp
C     produce form factor table by fourier-transformation
C     of the charge distribution

      implicit none
      include "tailcom.inc"
      real q2max,rmin,rmax,eps,hbc
      real q2
      real forintp,gauss
      external forintp
      integer nqbin1,iq

      DATA Q2MAX,RMIN,RMAX,EPS,NQBIN1/.3,1.E-6,20.,1.E-6,600/
      DATA HBC /197.3234/
C
      NQBIN=NQBIN1
C
C       loop over q-bins
      Q2BIN=Q2MAX/FLOAT(NQBIN-1)
C
      DO 10 IQ=1,NQBIN
C               first bin is q2=0.    last is q2=q2max
          Q2=Q2BIN*FLOAT(IQ-1)
C               calc 3-momentum forq (transform from gev to fermi?)
          QFOR=SQRT(Q2)*1000./HBC
C               perform fourier-transform by integ with CERN-Gauss routine
          FFNUC(IQ)=GAUSS(FORINTp,RMIN,RMAX,EPS)
C               normalise form factor in a way, that FFNUC(Q=0)=1.
          IF(IQ.GT.1)FFNUC(IQ)=FFNUC(IQ)/FFNUC(1)
c        if (iq.lt.20) print *,iq,q2,ffnuc(iq)
10    CONTINUE
      FFNUC(1)=1.
C
      RETURN
      END


        REAL FUNCTION FORINTp(R)
C       INTEGRAND OF FOURIER TRANSFORMATION INTEGRAL
C       THREE-DIMENSIONAL INTEGRATION IS REDUCED TO
C       AN INTEGRATION OVER THE RADIUS BY MULTIPLYING
C       WITH A BESSEL FKT.

      implicit none
      include "tailcom.inc"
      real r,qr,forrhop

        IF (QFOR.GT.1.E-5) THEN
           QR=QFOR*R
           FORINTp=FORRHOp(R)*SIN(QR)*R**2/QR
        ELSE
           FORINTp=FORRHOp(R)*R**2
        ENDIF
c        print *,'forintp=',forintp
        RETURN
        END


        REAL FUNCTION FORRHOp(R)
C       charge distribution as fkt. of radius
C     generalized 2-parameter fermi for not too light nuclei
C     Ref: T.de Forest and J.D.Walecka, Adv.in Phys.15(1966)1  -- page 21
C          they cite R.Hofstadter, Rev Mod Phys 28 (1963) 214
C         z = 2.4 / (4*ln(3))
 
      implicit none
      include "tailcom.inc"
      include "cmpcom.inc"
       real r,zpara,cpara,wpara

       if (ire.eq.14) then
C - 3 parameter Fermi charge distr of 14N
C H. de Vries et al, Atomic Data and Nuclear Tables 36, 495-536 (1987)
         CPARA=2.57
         ZPARA=0.5052
         WPARA=-0.180
         FORRHOp=(1.0+WPARA*(R/CPARA)**2)/
     1       (1.0+EXP((R-CPARA)/ZPARA))
       elseif (ire.eq.20) then
C - 3 parameter Fermi charge distr of 20 Ne
C H. de Vries et al, Atomic Data and Nuclear Tables 36, 495-536 (1987)
         CPARA=2.791
         ZPARA=0.698
         WPARA=-0.168
         FORRHOp=(1.0+WPARA*(R/CPARA)**2)/
     1       (1.0+EXP((R-CPARA)/ZPARA))
       elseif (ire.eq.84) then
C - RHO2BOMO (TERADNMC.CAR)
C     generalized 2-parameter fermi
C     Ref: Bohr & Mottelson, Nuclear Stucture (Benjamin, New York 1969)
C          cited in Date et al., PRL 52 (1984) 2344
         ZPARA=0.54
         CPARA=1.12*rTARA**(1./3.)-0.86/rTARA**(1./3.)
         FORRHOp=1./(1.+EXP((R-CPARA)/ZPARA))
       else
C   generalized 2-parameter fermi for not too light nuclei
C   Ref: T.de Forest and J.D.Walecka, Adv.in Phys.15(1966)1  -- page 21
C        they cite R.Hofstadter, Rev Mod Phys 28 (1963) 214
C        z = 2.4 / (4*ln(3))
         ZPARA=0.546144
         CPARA=1.07*rTARA**(1./3.)
         FORRHOp=1./(1.+EXP((R-CPARA)/ZPARA))
       endif
       RETURN
       END


C **********************************************************************
C
C     FFQUAS
C
C **********************************************************************

      subroutine ffquas(t,geun,gmun,gepo,gmpo)
      implicit real*8(a-h,o-z)

      include "cmpcom.inc"
      include "tailcom.inc"

      call ffpro(t,gep,gmp)
      tf=t/chbar**2
      tau=t/4./amp**2
      tau1=1.+tau
c
c   t. de forest and d.j. valecka adv. in phys. 15(1966) no.57
c
        if(ire.eq.3)then
         supele=1.
         supmag=1.
         if(t.le.(2.d0*fermom)**2)then
          sqrat=dsqrt(t)/fermom
          supele=0.75*sqrat-sqrat**3/16.
          supmag=supele
         endif
        else
         supele=supst(t) 
c            qbold=sqrt(t*tau1) 
c            supele=1.-dsbern(qbold)**2
             supmag=supele 
        endif
        geun=gep*dsqrt(supele*rtarz)
        rtarn=rtara-rtarz
        gmun=gep*dsqrt(supmag*(rtarz*amm**2+rtarn*amn**2))
        if(ire.eq.2)then
         gepo=geun
         gmpo=gmun
        else
         gepo=0.
         rtarn=rtara-rtarz
         gmpo=gep*dsqrt(supmag*(rtarn*amn**2))
        endif
        end

********************** supst ************************************
 
      double precision function supst(t)
      implicit real*8(a-h,o-z)
      data chbar/.197328d0/
c
c     tf is t in fermi**(-2)
c
      tf=t/chbar**2
c
c   s.stein et al phys. rev. 12(1975)1884 (appendix 1)
c
               sqtf=dsqrt(tf)
               delff=(datan(sqtf/.93d0)-2.*datan(sqtf/3.19d0)+
     +    datan(sqtf/5.45d0))*1.580d0/sqtf
               delff=dmax1(0.d0,delff)
               supele=1.-delff**2
               supst=dmax1(0.d0,supele)
       return
      end 

C **********************************************************************
C
C     DDEXP
C
C **********************************************************************

      double precision function ddexp(x)
      implicit real*8(a-h,o-z)

      ddexp=0.
      if(x.gt.-50.)ddexp=exp(x)

      return
      end

C **********************************************************************
C
C     DQN32
C
C **********************************************************************

      subroutine dqn32(xl,xu,fct,y)
c
c
      external fct
      double precision xl,xu,y,a,b,c,fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))

      return
      end

C **********************************************************************
C
C     INTBTK2
C
C **********************************************************************

      subroutine intbtk2(dbtk,nbtk,dtkmax)
      implicit real*8 (a-h,o-z)
c
c ***  choose 7 bins for integration over the theta k angle
      dimension  dbtk(8)

      include "sxycom.inc"
      include "cmpcom.inc"
      include "ppicom.inc"
      include "radgen.inc"

      data dc /3.5d0/
c
c      bins are chosen according to the s peak and p peak
c      width of the peaks:

      dwids=dsqrt(ap*aml/s)/dc
      dwidp=dsqrt(ap*aml/x)/dc
* avoid values above 1.0 (for x values below 1e-07, abr 19.02.04)
*     dts=acos( (s*sx+ap2*y)/sqls/sqly )
*     dtp=acos( (x*sx-ap2*y)/sqlx/sqly )
      dts=acos( min(1.0d0,(s*sx+ap2*y)/sqls/sqly) )
      dtp=acos( min(1.0d0,(x*sx-ap2*y)/sqlx/sqly) )

      dpi=pi

c
c      define bin boundaries
      dbtk(1)=0.d0
      dbtk(2)=dmin1(4.d0*dc*dwids,dts/3.0d0)
      dbtk(3)=dmax1(dts-dwids,dts/1.5d0)
      dbtk(4)=dmin1(dts+dwids,dts+(dtp-dts)/3.0d0)
      dbtk(5)=dmax1(dtp-dwidp,dts+(dtp-dts)/1.5d0)
      dbtk(6)=dmin1(dtp+dwidp,dtp+(dpi-dtp)/3.0d0)
      dbtk(7)=.98d0*dbtk(6)+.02d0*dpi
      dbtk(8)=dpi
c
c      check cut in theta which follows from infra red cut:
c
      nbtk=0
      do 10 ka=2,8
        nbtk=nbtk+1
        if(dbtk(nbtk+1).ge.dtkmax) goto 20
 10   continue
      goto 30
 20   dbtk(nbtk+1)=dmin1(dbtk(nbtk+1),dtkmax)
 30   continue
c
c      nbtk=number of valid bins
c
      return
      end

C **********************************************************************
C
C     INTTK2
C
C **********************************************************************

      subroutine inttk2(isumtk,dbtk1,dbtk2,dsumtk,derrtk)
      implicit real*8 (a-h,o-z)

c      integrate over ntk bins of brems-gamma theta angle dtk
      include "radgen.inc"
      include "sxycom.inc"
      include "cmpcom.inc"
      include "tailcom.inc"
      include "ppicom.inc"

      ddtk=(dbtk2-dbtk1)/ntk
      dtk=dbtk1-ddtk*.5d0
      dsum=0.d0
      derr1=0.d0
c     loop over all bins
      do 10 itk=1,ntk
        isumtk=isumtk+1
        dtk=dtk+ddtk
        itkcur=isumtk
c     calculate integrand
C       call intsai(dtk,dtsai)
C       dsum=dsum+dtsai
        ta=(sx-sqly*cos(dtk))/ap2
        dsigma=an*alfa/pi*rv2(ta) * sin(dtk)*sqly/ap2
        dsum=dsum+dsigma
C write(9,*)dsigma,ta,dtk
        if(ita.le.3)then
          dsitkm(isumtk,ita) = dsumtk + ddtk * dsum
          dcvtkm(isumtk,ita) = dtk
          ddetkm(isumtk,ita) = ddtk
        endif
c      add up errors:
c      errors of theta integration:
        if    (itk.gt.2.and.itk.lt.ntk) then
          derr1=derr1+dabs(dsigma-dold)
        elseif(itk.eq.2.or.itk.eq.ntk) then
          derr1=derr1+1.5d0*dabs(dsigma-dold)
        endif
c     
        dold=dsigma
c
 10   continue
c
c      integral on big bin is:
C     print *,dsumtk,ddtk,dsum
      dsumtk=dsumtk+ddtk*dsum
      derrtk=derrtk+ddtk*derr1
c
      return
      end


C **********************************************************************
C
C     ITAFUN
C
C **********************************************************************

      subroutine itafun(ys,xs,pl,pn,ixytest,itagen)

      include "cmpcom.inc"
      include "radgen.inc"
      include "density.inc"
      include "xytabcom.inc"

      dimension xym(2),na(2),a(2*nt)

      do nny=nty,1,-1
C        print *, nny, nty, nt, ys, y(nny)
        if (ys.gt.y(nny)) goto 10
      enddo
 5    print *,' ys=',ys,' out of the region'
      stop
 10   if (nny.eq.nt) goto 5

      if(ixytest.eq.1)then
        stop ' ixytest.eq.1 is not supported in the version'
      else
         xym(1)=xs
         xym(2)=ys
         na(1)=ntx
         na(2)=nty
        do i=1,ntx
          a(i)=x(i)
          a(i+ntx)=y(i)
        enddo
* what nonsense if this ? - abr 19.02.04
* none of the polarised quantities (tbor_p etc) is ever different from zero !
* use final quantities (tbor etc) immediately
* different array sizes for DIS and phoitoproduction kinematics
* also include other quantities as redfac to be available with LUT
      if (ntx.eq.ntdis) then
        tbor=fint(2,xym,na,a,tbor_udis)
        sigrad=fint(2,xym,na,a,sigrad_udis)
        sig1g=fint(2,xym,na,a,sig1g_udis)
        vac=fint(2,xym,na,a,vac_udis)
        vertex=fint(2,xym,na,a,vertex_udis)
        small=fint(2,xym,na,a,small_udis)
        redfac=fint(2,xym,na,a,redfac_udis)
      else
        tbor=fint(2,xym,na,a,tbor_upho)
        sigrad=fint(2,xym,na,a,sigrad_upho)
        sig1g=fint(2,xym,na,a,sig1g_upho)
        vac=fint(2,xym,na,a,vac_upho)
        vertex=fint(2,xym,na,a,vertex_upho)
        small=fint(2,xym,na,a,small_upho)
        redfac=fint(2,xym,na,a,redfac_upho)
      endif
*       tbor=tboru+pl*pn*tborp
*       sigrad=sigradu+pl*pn*sigradp
*       sig1g=sig1gu+pl*pn*sig1gp

      endif

      sigcor=sigrad/sig1g
!|bs>
cc      write(*,*)'radgen:',xs,ys,sigcor,sigrad,sig1g,sig1gu,pl,pn,sig1gp

      r1=rlu(0)
      rr1=r1*sigrad

* same nonsense here (no need for tineu and tinep, use tine immediately)
        if (ntx.eq.ntdis) then
          tine=fint(2,xym,na,a,tine_udis)
          tnuc=fint(2,xym,na,a,tnuc_udis)
          tpro=fint(2,xym,na,a,tpro_udis)
        else
          tine=fint(2,xym,na,a,tine_upho)
          tnuc=fint(2,xym,na,a,tnuc_upho)
          tpro=fint(2,xym,na,a,tpro_upho)
        endif

      if(rr1.gt.sigrad-tbor)then
        itagen=0
        return
      endif

C       write(*,'(6g13.5)')sigrad,tbor,tine,tpro,tnuc,rr1
c      pause

      if(rr1.gt.(tpro+tnuc)) then
        itagen=1
c        scgen=rr1-tpro-tnuc

      elseif(rr1.gt.tnuc)then
        itagen=3
c       scgen=rr1-tnuc
      else
c       scgen=rr1
        itagen=2
      endif

! Speed optimization by Joe Seele (seele@uiuc.edu) June 4, 2002

C #undef OLDCODE
      do ii=1,7*ntk
C #ifdef OLDCODE
C         do ix=1,ntx
C           do iy=1,nty
C             if (ntx.le.ntdis) then
C               wrk(ix,iy)=densdis(ntx+1-ix,iy,ii,itagen)
C             else
C               wrk(ix,iy)=denspho(ntx+1-ix,iy,ii,itagen)
C             endif
C           enddo
C         enddo
C         dsitkm(ii,itagen)=fint(2,xym,na,a,wrk)
C #else
      if (ntx.le.ntdis) then
        dsitkm(ii,itagen)=fint(2,xym,na,a,densdis(1,1,ii,itagen))
      else
        dsitkm(ii,itagen)=fint(2,xym,na,a,denspho(1,1,ii,itagen))
      endif
C #endif
      enddo
      
      do ii=1,7
C #ifdef OLDCODE
C         do ix=1,ntx
C           do iy=1,nty
C             if (ntx.le.ntdis) then
C               wrk(ix,iy)=widdis(ntx+1-ix,iy,ii,itagen)
C             else
C               wrk(ix,iy)=widpho(ntx+1-ix,iy,ii,itagen)
C             endif
C           enddo
C         enddo
C         ddetkm(ii,itagen)=fint(2,xym,na,a,wrk)
C #else
        if (ntx.le.ntdis) then
          ddetkm(ii,itagen)=fint(2,xym,na,a,widdis(1,1,ii,itagen))
        else
          ddetkm(ii,itagen)=fint(2,xym,na,a,widpho(1,1,ii,itagen))
        endif
C #endif
      enddo

      do iittkk=1,7*ntk
        ssuumm=0.
        do itoo=1,iittkk/ntk
          ssuumm=ssuumm+ddetkm(itoo,itagen)*ntk
        enddo
        dcvtkm(iittkk,itagen)=ssuumm+ddetkm((iittkk-1)/ntk+1,itagen)
     +       *(mod(iittkk,ntk)-0.5)
      enddo

      end



C **********************************************************************
C
C     XYTABL
C
C **********************************************************************

      subroutine xytabl(tname,e1,plrun,pnrun,ixytest,ire)

      include "cmpcom.inc"
      include "radgen.inc"
      include "density.inc"
      include "xytabcom.inc"
      include "mc_set.inc"
      
      character*256 tname

      data lun/66/
      open (unit=lun,file=tname)

      if(ixytest.eq.0)then

        do iy=1,nty
          y(iy)=radgen_ymin+(radgen_ymax-radgen_ymin)/(nty-1)*(iy-1)
        enddo

*abr      step=(radgen_xmax-radgen_xmin)/(nt-1)
* abr 19.02.04
        do ix=1,ntx
* start from high x, then bring in ascending order
          if (ix.le.18) then
            step=0.05
            x(ix)=1.0-step*ix
          elseif (ix.le.27) then
            step=0.01
            x(ix)=x(18)-step*(ix-18)
          elseif (ix.le.36) then
            step=0.001
            x(ix)=x(27)-step*(ix-27)
          elseif (ix.le.48) then
            if (mod(ix,2).gt.0) then
              x(ix)=x(ix-1)/2
            else
              x(ix)=x(ix-1)/5
            endif
          endif
        enddo
        CALL FLPSOR(x,ntx)

        do iy=1,nty
          do ix=1,ntx
*abr        x(ix,iy)=radgen_xmin+step*(ix-1)
            bmp=amp
            bmc2=1.16
            ylim=(bmc2-bmp**2)/(2.*bmp*e1*(1d0-x(ix)))
            ys=y(iy)
            if(ys.lt.ylim ) ys=ylim
            call mpolrad(e1,ys,x(ix),1.,plrun,pnrun,-1)
            tbor_u(ix,iy)=tbor
            sig1g_u(ix,iy)=sig1g
            sigrad_u(ix,iy)=sigrad
            tine_u(ix,iy)=tine
            tnuc_u(ix,iy)=tnuc
            tpro_u(ix,iy)=tpro
            vac_u(ix,iy)=vac
            vertex_u(ix,iy)=vertex
            small_u(ix,iy)=small
            redfac_u(ix,iy)=redfac
C #ifdef TEST
C             write(*,'(a3,6g12.4)')' u ',
C      +           y(iy),x(ix),sig1g,tnuc,tpro,sigrad
C #endif
            write(lun,'(1x,14g14.4)')x(ix),y(iy)
     +           ,sig1g_u(ix,iy),sigrad_u(ix,iy),tbor_u(ix,iy)
     +           ,tine_u(ix,iy),tnuc_u(ix,iy),tpro_u(ix,iy)
     +           ,vac_u(ix,iy),vertex_u(ix,iy),small_u(ix,iy)
     +           ,redfac_u(ix,iy)

            do itkm=0,6
              width(ix,iy,itkm+1,1)=ddetkm(itkm*ntk+1,1)
              width(ix,iy,itkm+1,2)=ddetkm(itkm*ntk+1,2)
              if(ire.ne.1)
     +             width(ix,iy,itkm+1,3)=ddetkm(itkm*ntk+1,3)
            enddo
            
            ittmax=3
            if(ire.eq.1)ittmax=2
            write(lun,*)((width(ix,iy,itkm+1,itt)
     +           ,itkm=0,6),itt=1,ittmax),ndxtkm
            
            do itkm=1,ntk*7
* avoid values smaller than E-24, happens in case of N14 (abr)
* similar problem seems to there at high x values
* is this really the solution ???
* at least it seems to work....
*             denstk(ix,iy,itkm,1)=dsitkm(itkm,1)/dsitkm(ntk*7,1)
              if (dsitkm(itkm,1).lt.1.E-24 .or.
     &            dsitkm(ntk*7,1).lt.1E-24 ) then
                denstk(ix,iy,itkm,1)=1.0
              else
                denstk(ix,iy,itkm,1)=dsitkm(itkm,1)/dsitkm(ntk*7,1)
              endif

              if (dsitkm(itkm,2).lt.1.E-24 .or.
     &            dsitkm(ntk*7,2).lt.1E-24 ) then
                denstk(ix,iy,itkm,2)=1.0
              else
                denstk(ix,iy,itkm,2)=dsitkm(itkm,2)/dsitkm(ntk*7,2)
              endif

c             if (ix.gt.20) then
c              write(6,*) dsitkm(itkm,2),dsitkm(ntk*7,2),
c    +                    denstk(ix,iy,itkm,2)
c             endif
* 28.01.04 - better include safety limits here as well
*             if(ire.ne.1)
*    +             denstk(ix,iy,itkm,3)=dsitkm(itkm,3)/dsitkm(ntk*7,3)
              if(ire.ne.1) then
              if (dsitkm(itkm,3).lt.1.E-24 .or.
     &            dsitkm(ntk*7,3).lt.1E-24 ) then
                denstk(ix,iy,itkm,3)=1.0
              else
                denstk(ix,iy,itkm,3)=dsitkm(itkm,3)/dsitkm(ntk*7,3)
              endif
              endif
            enddo
            write(lun,'(35g14.4)')(denstk(ix,iy,itkm,1)
     +           ,itkm=1,7*ntk)
            write(lun,'(35g14.4)')(denstk(ix,iy,itkm,2)
     +           ,itkm=1,7*ntk)
            if(ire.ne.1)
     +           write(lun,'(35g14.4)')(denstk(ix,iy,itkm,3)
     +           ,itkm=1,7*ntk)
          enddo
        enddo
        close(lun)
        
      elseif(ixytest.gt.0)then
        do iy=1,nty
          do ix=1,ntx
C #ifdef TEST
C             print *,iy,ix
C #endif
            if (ntx.eq.ntdis) then    ! DIS range (xmin=0.002)
              read(lun,*)x(ix),y(iy)
     +           ,sig1g_udis(ix,iy),sigrad_udis(ix,iy),tbor_udis(ix,iy)
     +           ,tine_udis(ix,iy),tnuc_udis(ix,iy),tpro_udis(ix,iy)
     +           ,vac_udis(ix,iy),vertex_udis(ix,iy)
     +           ,small_udis(ix,iy),redfac_udis(ix,iy)
c             print *,ix,iy
c    +           ,sig1g_u(ix,iy),sigrad_u(ix,iy),tbor_u(ix,iy)
c    +           ,tine_u(ix,iy),tnuc_u(ix,iy),tpro_u(ix,iy)

              ittmax=3
              if(ire.eq.1)ittmax=2
              read(lun,*)
     +           ((widdis(ix,iy,ite,itt),ite=1,7),itt=1,ittmax),ndxtkm
              read(lun,*)(densdis(ix,iy,itkm,1)
     +           ,itkm=1,7*ntk)
              read(lun,*)(densdis(ix,iy,itkm,2)
     +           ,itkm=1,7*ntk)
              if(ire.ne.1)
     +           read(lun,*)(densdis(ix,iy,itkm,3)
     +           ,itkm=1,7*ntk)
            else                    ! photoproduction (xmin=1e-07)
              read(lun,*)x(ix),y(iy)
     +           ,sig1g_upho(ix,iy),sigrad_upho(ix,iy),tbor_upho(ix,iy)
     +           ,tine_upho(ix,iy),tnuc_upho(ix,iy),tpro_upho(ix,iy)
     +           ,vac_upho(ix,iy),vertex_upho(ix,iy)
     +           ,small_upho(ix,iy),redfac_upho(ix,iy)

              ittmax=3
              if(ire.eq.1)ittmax=2
              read(lun,*)
     +           ((widpho(ix,iy,ite,itt),ite=1,7),itt=1,ittmax),ndxtkm
              read(lun,*)(denspho(ix,iy,itkm,1)
     +           ,itkm=1,7*ntk)
              read(lun,*)(denspho(ix,iy,itkm,2)
     +           ,itkm=1,7*ntk)
              if(ire.ne.1)
     +           read(lun,*)(denspho(ix,iy,itkm,3)
     +           ,itkm=1,7*ntk)
            endif

          enddo
        enddo
      else
       stop 'xytabl -- ixytest < 0'
      endif
      close(lun)

      end

C **********************************************************************
C
C
C
C **********************************************************************
