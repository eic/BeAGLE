
      subroutine radgen_event

      implicit none

       include "mc_set.inc"
       include "mconsp.inc"
       include "phiout.inc"
       include "tailcom.inc"
       include "cmpcom.inc"
       include "radgen.inc"
       include "mcRadCor.inc"

      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      INTEGER MSTP,MSTI
      REAL PARP,PARI
      SAVE /PYPARS/

      real PhRAD(4),q2true,nutrue,radweight
      real nu,q2,phi,yys,xxs

! calculate radiative corrections
      nu=sngl(gennu)
      q2=sngl(genq2)
      phi=sngl(genphi)
      yys=sngl(geny)
      xxs=sngl(genx)
      call RADGEN(mcSet_EneBeam,q2,nu,yys,xxs,phi,PhRAD,q2true,nutrue,
     +     radweight)

! fill mcRadCor WCB with ...

      mcRadCor_ID=1

! ... true kinematics

* by definition we calculate xbj using the proton mass
* ---> xbj for elastic events is A
*     mcRadCor_XTrue=q2true/(2.0d0*amp*nutrue)
      mcRadCor_NuTrue=nutrue
      mcRadCor_Q2True=q2true
      mcRadCor_YTrue=nutrue/mcSet_EneBeam
      mcRadCor_XTrue=q2true/mcRadCor_YTrue/(4.*ebeam*pbeam)
      mcRadCor_W2True=amp2 + (q2true*(1./mcRadCor_XTrue-1.))
!      mcRadCor_XTrue=q2true/(2.0d0*0.938272d0*nutrue)
!      mcRadCor_W2True=amp2 - q2true + 2.*amp*nutrue

! ... kinematics of real photon
      mcRadCor_EBrems=phrad(4)
      mcRadCor_ThetaBrems=0.
      if(phrad(4).gt.0.)
     +    mcRadCor_ThetaBrems = acos(phrad(3)/phrad(4))
      mcRadCor_PhiBrems=0.
      if (.not.(phrad(1).eq.0..and.phrad(2).eq.0.)) then
        mcRadCor_PhiBrems = atan2(phrad(2),phrad(1))
        if (mcRadCor_PhiBrems.lt.0.)
     +       mcRadCor_PhiBrems = mcRadCor_PhiBrems + twopi
      endif

c...if we would like to have the TSAI system angles for the real gamma than
c     mcRadCor_ThetaBrems = dthg
c     mcRadCor_PhiBrems = dphig

! ... radiative contributions

      mcRadCor_Sigrad=sigrad
      mcRadCor_Sigcor=sigcor
      mcRadCor_Sigcorerr=0.
      mcRadCor_TailIne=tine
      mcRadCor_TailEla=tpro
      mcRadCor_TailCoh=tnuc
      mcRadCor_Vacuum=vac
      mcRadCor_Vertex=vertex
      mcRadCor_Small=small
      mcRadCor_Redfac=redfac

! ... radiative correction type

      if (ita.eq.2) then
        mcRadCor_cType='elas'
      else if (ita.eq.3) then
        mcRadCor_cType='qela'
      else
        mcRadCor_cType='inel'
      endif

      end
