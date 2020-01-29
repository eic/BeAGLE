* These are interfaces to rapgap needed by dpmjet
*
* Written by Liang Zheng and Mark Baker
*
*=====dt_rginitep=========================================================
*...Initialize rapgap when using dpmjet
* Adapted from RAPGAP main program rgmain.F 
*=======================================================================
      SUBROUTINE DT_RGINITEP(PBEAML,PBEAMH,Q2MIN,Q2MAX,YMIN,YMAX,
     &              THMIN,THMAX,INPUT)
*     input:  (all energy units GeV. c=1.)
*           PBEAML   lepton beam momentum in lab frame 
*           PBEAMH   hadron beam momentum in lab frame
*           Q2MIN    Q2 cut low
*           Q2MAX    Q2 cut high
*           YMIN     Y cut low
*           YMAX     Y cut high
*           THMIN    Scattered electron angle cut low (in radians)   
*           THMAX    Scattered electron angle cut high (in radians)
*           INPUT    Input file   
C Note: For now all but INPUT are dummy variable. Later we will use them.

      Implicit None
      include 'rgluco.inc'  !  For IFPS in /INPU/
      include 'beagle.inc'

C...Pythia event counter (since we keep PYINITing) Mark 2017-01-31
      COMMON /PYCNTR/ MYNGEN
      INTEGER MYNGEN
      SAVE /PYCNTR/

      CHARACTER*8 INPUT
      integer LINP
      parameter ( LINP=28 )
      CHARACTER*256 inputfilename

      real timeleft
      Integer Minuts
      External Minuts
      External pydata

      Integer I

      Double precision PBEAML, PBEAMH, Q2MIN, Q2MAX, YMIN, YMAX
      Double precision THMIN, THMAX

C...Switches for nuclear correction
      COMMON /PYNUCL/ INUMOD,CHANUM,ORDER,genShd
      SAVE /PYNUCL/
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER,genShd

* properties of interacting particles
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5) 
      INTEGER IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,ITMMOD
      INTEGER MODHYP,NHYPER,IDHYP 
C...Parameters and switch for energy loss
      DOUBLE PRECISION QHAT
      INTEGER QSWITCH
      COMMON /QUENCH/ QHAT, QSWITCH

C...Pythia eA shadowing common block from Mark 2017-06-30
      COMMON /PYSHAD/ NKNOTS,RDUMMY,RAVAL,XKNOT(100),YKNOT(100),
     &BKNOT(100),CKNOT(100),DKNOT(100),EKNOT(100),FKNOT(100),SHDFAC
      SAVE /PYSHAD/
      DOUBLE PRECISION XKNOT,YKNOT,BKNOT,CKNOT,DKNOT,EKNOT,FKNOT,RAVAL,
     &SHDFAC
      INTEGER NKNOTS,RDUMMY

C... Locals
      INTEGER IFSEED, initseed, idum1, idum2
      DOUBLE PRECISION Mneut, Mprot, Mlept, Mnucl, pbeamE, ebeamE
      DOUBLE PRECISION ebeamEnucl, epznucl,sqrts,altpbeam,altpbeamE
      DOUBLE PRECISION altsqrts

C... External
      DOUBLE PRECISION PYMASS, AZMASS
      EXTERNAL PYMASS, AZMASS

      WRITE (*,*)'Inputs: ',PBEAML, PBEAMH, Q2MIN, Q2MAX, YMIN, YMAX,
     &  THMIN, THMAX, "Filename: ", INPUT

C - TEMPORARY - Hard-code parameters we will want to read in later
      genShd = 1
      ORDER = 101
      QSWITCH = 1
      QHAT = 0.0D0
      SHDFAC = 1.0D0
C      IFSEED = 1234567
      IFSEED = 44788029
      INUMOD=IT
      CHANUM=ITZ
      initseed = IFSEED
      call rndmq (idum1,idum2,initseed,' ')


C     Fill /LABTONR/ & /BEAEVT/ variables for use by DT_PYEVNTEP & DT_PYOUTEP
C      COMMON /LABTONR/ MAscl,pgamma,pbeta,pbeamP,pbeamN,pbeamdbl,
C     &                 idNucPY,idNucBAM, lName
       pbeamP=PBEAMH
       pbeamN=PBEAMH
       pbeamdbl=PBEAMH
       if( IJPROJ.eq.3) then
          ltype = 11
       elseif( IJPROJ.eq.4 ) then
          ltype = -11
       elseif( IJPROJ.eq.10 ) then
          ltype = -13
       elseif( IJPROJ.eq.11 ) then
          ltype = 13
       endif

      !set up nucleon beam
      Mneut=PYMASS(2112)
      Mprot=PYMASS(2212)
      Mlept=PYMASS(ltype)

      if( IJTARG.eq.1 ) then 
         Mnucl=Mprot
         idNucPY=2212
         idNucBAM=1
      elseif ( IJTARG.eq.8 ) then
         Mnucl=Mneut
         idNucPY=2112
         idNucBAM=8
      else
         write (*,*) "Unsupported target. IJTARG= ", IJTARG
         STOP "Aborting"
      endif

C MDB 2016-11-10 For A>1, match rapidities for nucleon and nucleus 
C Note: massp means nucleon mass as seen in radgen.
      IF (IT.GT.1) THEN
         MAscl=AZMASS(IT,ITZ,ITMMOD)/INUMOD 
         pbeamdbl=PBEAMH*Mnucl/MAscl
         pbeamP=(PBEAMH*Mprot/MAscl)
         pbeamN=(PBEAMH*Mneut/MAscl)
      ENDIF

C     proton (neutron) is defined in positive z and as target
C      P(2,1)=0.0  
C      P(2,2)=0.0
C      P(2,3)= pbeamdbl
      PZTARG = PBEAMH
      PZNUCL = pbeamdbl
C     lepton is defined in negative z and as beam
C      P(1,1)=0.0  
C      P(1,2)=0.0  
C      P(1,3)=-PBEAML
      PZLEP = -PBEAML

      pbeamE=sqrt(pbeamdbl**2+Mnucl**2)
      pbeta=pbeamdbl/pbeamE
      pgamma=pbeamE/Mnucl
      ebeamE=sqrt(PBEAML**2+Mlept**2)
      ebeamEnucl=pgamma*ebeamE-pgamma*pbeta*(-PBEAML)
      epznucl=-pgamma*pbeta*(ebeamE)+pgamma*(-PBEAML)
      write(*,*) ebeamEnucl, ebeamE, epznucl, -PBEAML

      sqrts=sqrt((pbeamE+ebeamE)**2-(pbeamdbl-PBEAML)**2)
      altpbeam=PBEAMH
      altpbeamE=sqrt(PBEAMH*PBEAMH+PYMASS(2212)**2)
      altsqrts=sqrt((altpbeamE+ebeamE)**2-(altpbeam-PBEAML)**2)
      write(*,*) '*********************************************'
      IF(IT.GT.1) THEN
         write(*,*) 'Nucleus: Z, A: ', ITZ, ' ', IT
         write(*,*) 'Per nucleon beam momentum:', PBEAMH, 'GeV'
         write(*,*) 'sqrts for momentum matched ep: ',altsqrts,'GeV'
         write(*,*) 'Rapidity-matched neutron beam momentum: ', 
     &               pbeamN, 'GeV'
         write(*,*) 'Rapidity-matched proton beam momentum: ', 
     &               pbeamP, 'GeV'
      ELSEIF (ITZ.GE.1) THEN
         write(*,*) 'proton beam momentum: ', pbeamdbl, 'GeV'
      ELSE
         write(*,*) 'neutron beam momentum: ', pbeamdbl, 'GeV'
      ENDIF
      write(*,*) 'lepton beam momentum: ', PBEAML, 'GeV'
      write(*,*) 'sqrts for eN system: ', sqrts, 'GeV'
      write(*,*) 'lepton beam energy in TRF: ',ebeamEnucl,'GeV'
      write(*,*) '*********************************************'

C     GenNucDens is used even without quenching.
      call GenNucDens(ITZ, IT)

C      IF(QSWITCH.EQ.1) THEN
C         print*,'Quenching requested, qhat=',QHAT
C         print*,'Using recoil factor:',PQRECF
C         print*,'      Ptmodel iPtF: ',PYQ_IPTF
C         print*,'      EmitGluon iEg:',PYQ_IEG
C         print*,'      SupFactor:    ',PYQ_SUPF
C         print*,'      Heavy Quarks: ',PYQ_HQ
C         print*,'      Energy Thres: ',PYQ_IET
Cc...when quenching is used switch off internal parton shower
C         MSTP(61)=0
C         MSTP(71)=0
C      ENDIF

C     Force weak decays to happen outside the nucleus
      CALL DT_PYDECY(1)
      MYNGEN=0

C---initialise ARIADNE parameters, now done via block data
C---initialise PYTHIA 6 parameters, via pythia block data 
C     initialize random number generator
C      ISEED = 213123
C      ISEED = Iabs(MINUTS())
C
C NOTE: This was a duplication of the rndmq call above
C       since it calls RLUXGO too.
C      ISEED = 44788029
C      LUX = 4
C      K1=0
C      K2=0
C      CALL RLUXGO(LUX,ISEED,K1,K2)
C---initialise RAPGAP parameters
      CALL GRAINI
C-- read in parameters from file 
      inputfilename=INPUT
      write(*,*) 'the input file is: ', inputfilename
      open(LINP, file=inputfilename,STATUS='UNKNOWN')
      Call rgsteer
C-- change standard parameters of RAPGAP 	
      Call rapcha
C-- change standard parameters of HERACLES	
      Call hercha
C-- change standard parameters of JETSET
      Call pytcha
      IF(IFPS.EQ.10) then
C Initialize ARIADNE
         CALL ARINIT('RAPGAP')
C-- change standard parameters of ARIADNE
         Call aricha
      endif

C---  CALCULATE X SECTION
      CALL PTIME(' rapgap  ',1,0)
      CALL RAPGAP
C TEMPTEMP
         write(*,*) 'In DT_RGINITEP. Ran RAPGAP.'
C---  print x section
      CALL RAEND(1)

      RETURN
      END

*=====dt_rgoutep=========================================================
*used for the output of pythia event list and statistics information
* Adapted from B. Page's version of rgmain.F and analys.F along with 
* DT_PYOUTEP
      SUBROUTINE DT_RGOUTEP(MODE)     

*     input:
*           MODE: 1:reject statistics - not really used
*                 2:event output
*                 3:total statistics print  
*=======================================================================

      Implicit none 
      INTEGER MODE

      include 'rgrapgki.inc'
      include 'rgluco.inc'
      include 'beagle.inc'

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 
      INTEGER NEVENT,ICASCA

C. Could also use PARU(2)
      double precision twopi
      parameter (twopi=6.283185307179586d0)

C...Pythia Commons. Spell out to avoid IMPLICIT NONE clash w/ pythia.inc. 
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      INTEGER MINT,MSTP,MSTI
      DOUBLE PRECISION VINT,PARP,PARI

* added by liang to store the output event variables 1/20/12
      COMMON /EVTOUT/ XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT
      DOUBLE PRECISION XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT

C...output file name definition
      COMMON /OUNAME/ outname
      CHARACTER*256 outname

*KEEP, RGPART. change DBCMS(4)->DBCMSS(4)
      DOUBLE PRECISION SSS,CM,DBCMSS
      COMMON /PARTON/ SSS,CM(4),DBCMSS(4)

*From pythia.inc - straight copy.
      INTEGER N,NPAD,K
      Double Precision P,V
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)

*KEEP, RGPARA1. change q2--> q2r
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2R,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2R,Q2Q,PCM(4,18)

*KEEP, RGPARAS. change IPRO --> IPROR
      DOUBLE PRECISION Q2START,W_Q2,OMEG2
      INTEGER IRUNA,IQ2,IRUNAEM
      INTEGER IPROR
      COMMON/RAPA /IPROR,IRUNA,IQ2,IRUNAEM,Q2START,W_Q2,OMEG2
      DOUBLE PRECISION SCALFA
      COMMON/SCALF/ SCALFA

C*KEEP, RGDISDIF. change IDIR --> IDIRR
      INTEGER IDIRR,IDIRINT,IDISDIF
      COMMON/DISDIF/ IDIRR,IDIRINT,IDISDIF

C*KEEP, RGEFFIC. change NOUT --> NOUTT
      DOUBLE PRECISION AVGI,SD
      INTEGER NIN,NOUTT
      COMMON/EFFIC/AVGI,SD,NIN,NOUTT

      DOUBLE PRECISION DBCMS(4)
      Integer NOUT,NOU110,NOU210,NOU1100
      COMMON/MYOUT/NOUT,NOU110,NOU210,NOU1100
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN

C     Local
      DOUBLE PRECISION  DOT
      real pt2_gen,phi_gen
      Double Precision PYP,pyangl
      Double Precision sphi,stheta,phit1,phit2
      Integer ncall,naf1,naf2,l,ll
c  ntuple gluon
      Real x_bj,y,xgam,Q2,x_pom,t,sigm
      Integer ipro,idir
      Real s_h,p_t,p_t1,p_t2,p_t3,phi1,phi2,phi3,eta1,eta2,
     + eta3,xg_p,xm2_tot
      common/cwnco1/x_bj,y,xgam,Q2,x_pom,t,ipro,idir,sigm
      common/cwnco2/s_h,p_t,p_t1,p_t2,p_t3,phi1,phi2,phi3,eta1,eta2,
     + eta3,xg_p,xm2_tot

      real phi_e,phi_p
      Integer isprot,i,ipom
      double precision bochck
      Integer ievent, genevent, lastgenevent

      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
      EXTERNAL DOT

C definitions for outputfile and so on
      integer txtLun
      parameter (txtLun=29)
      character outputfilename*80
      parameter (outputfilename='rapgap.txt')
      outname = outputfilename

      IF (MODE.LT.1 .OR. MODE.GT.3) THEN
         WRITE(*,*) 'DT_RGOUTEP: Illegal mode: ', MODE
         STOP
      ENDIF

      GOTO (1,2,3) MODE

c...mode 1 is used to update the reject statistics - not implemented
1     CONTINUE      
c      write(99,*),'event ',ievent,' rejected,',' proces=',
c     & msti(1),', X=',XBJOUT,' Q2=',Q2OUT 
         
      RETURN

c...mode 2 is used to output the event list
2     CONTINUE

      WRITE (*,*) 'Temporary PYLIST inside of DT_RGOUTEP'
      CALL PYLIST(2)
      IF (.NOT. OLDOUT) THEN
C     Fix a few missing things in RAPGAP
         K(1,4) = 3
         K(1,5) = 4
C     Writeout RAPGAP in BeAGLE format
         CALL DT_PYOUTEP(MODE)     
      ELSE
        x_bj = SNGL(Q2R/DBLE(YY)/SSS)
        y = yy
        xgam = xel/yy
        Q2 = SNGL(Q2R)
        x_pom = XFGKI
        t = T2GKI
        ipro = ipror
        idir = idirr
        sigm = sngl(avgi)
        s_h = shh
        if(ipror.eq.12) shh = 0.
        p_t = sqrt(pt2h)
        phi_e = pyangl(P(4,1),P(4,2))
        do I=1,n
           if(k(i,2).eq.2212.and.k(i,1).eq.1) isprot =i
        enddo
        phi_p=pyangl(P(isprot,1),P(isprot,2))
C     
C TEMPTEMP
C
        write(*,*)'x_bj(M=m=0), x_bj(Q2,y,s,M,m), XBJOUT',x_bj,' ',
     &        Q2R/DBLE(YY)/(SSS-VINT(4)*VINT(4)-VINT(303)),' ',XBJOUT
        write(*,*)'y,Q2: ',x_bj,' ',y,' ',Q2
        write(*,*)'xgam,x_pom,t,p_t: ',xgam,' ',x_pom,' ',t,' ',p_t
         
C   Writeout Event in RAPGAP format.
        IF (FIRST) THEN
          open(txtLun, file=outname,STATUS='UNKNOWN')
          write(txtLun,*)' RAPGAP EVENT FILE ' 
          write(txtLun,*)'============================================'
          write(txtLun,40)
 40       format('I, ievent, genevent, subprocess, idir, idisdif, c.s.,'
     &         ,'sigma(c.s.), s, q2, y, xgam, xpr, Pt_hat, pt2_hat,',  
     &         's_hat, t, x_pom, s_hat, z, x, phi, nrTracks')
          write(txtLun,*)'============================================'
          write(txtLun,*)' I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)',
     &         '  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)',
     &         '  V(I,1)  V(I,2)  V(I,3)'
          write(txtLun,*)'============================================'
          NEVENT = 0
          NOUT = 0
          FIRST=.FALSE.
        ENDIF
       
c     write(6,*) ' analys: phi_e/p ',phi_e,phi_p
c     looking for phi asymmetries
C     boost in gamma proton system
C     and rotate
        IF(IDIRR.EQ.1) THEN
           DBCMS(1)= DBLE(P(NIA1,1) + P(2,1))
           DBCMS(2)= DBLE(P(NIA1,2) + P(2,2))
           DBCMS(3)= DBLE(P(NIA1,3) + P(2,3))
           DBCMS(4)= DBLE(P(NIA1,4) + P(2,4))
           xfgki=1.
        ELSE
c     find pomeron
           ipom = 2
           do i=1,n
              if(k(i,2).eq.100) ipom=i
           enddo
c     boost to gamma pomeron system
           DBCMS(1)= DBLE(P(NIA1,1) + P(ipom,1))
           DBCMS(2)= DBLE(P(NIA1,2) + P(ipom,2))
           DBCMS(3)= DBLE(P(NIA1,3) + P(ipom,3))
           DBCMS(4)= DBLE(P(NIA1,4) + P(ipom,4))
        ENDIF
c     write(6,*) ' analys nia1 ',nia1
        BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2
     +       + (DBCMS(3)/DBCMS(4))**2
        BOCHCK = DSQRT(BOCHCK)
        IF(BOCHCK.GT.0.99999999D0) then
           write(6,*) BOCHCK,nia1,ipom
           call pylist(1)
        endif
        CALL PYROBO(0,N,0.d0,0.d0,-DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +       -DBCMS(3)/DBCMS(4))
        SPHI = pyangl(P(nia1,1),P(nia1,2))
        CALL PYROBO(0,0,0.d0,-sphi,0.d0,0.d0,0.d0)
        STHETA = pyangl(P(nia1,3),P(nia1,1))
        CALL PYROBO(0,0,-STHETA,0.d0,0.d0,0.d0,0.d0)
        IF(P(NF1,4).GT.P(NF2,4)) THEN
           phit1 = pyangl(P(NF1,1),P(NF1,2))
           phit2 = pyangl(P(NF2,1),P(NF2,2))
        ELSE
           phit2 = pyangl(P(NF1,1),P(NF1,2))
           phit1 = pyangl(P(NF2,1),P(NF2,2))
        ENDIF
        if(phit1.lt.0.) phit1 = twopi + phit1
        if(phit2.lt.0.) phit2 = twopi + phit2
        
        if(nf1.ne.nf2) then
           naf1=nf1
           naf2=nf2
        else
           do 20  l=nf1+1,n-1
              if(k(l,1).eq.12.or.k(l,1).eq.2) THEN
                 do 10 ll = l,n
                    IF(k(ll,1).eq.11.or.k(ll,1).eq.1) then
                       naf1=l
                       naf2=ll
                       goto 30
                    endif
 10              continue
              endif
 20        continue
 30        continue
        endif
c     call pylist(1)
c     write(6,*) ' nf1,nf2,naf1,naf2 ',nf1,nf2,naf1,naf2
        p_t1 = PYP(naf1,10)
        p_t2 = PYP(naf2,10)
        p_t3 = PYP(naf1+1,10)
        phi1 = pyangl(P(NaF1,1),P(NaF1,2))
        phi2 = pyangl(P(NaF2,1),P(NaF2,2))
        phi3 = pyangl(P(NaF1+1,1),P(NaF1+1,2))
        eta1 = PYP(naf1,19)
        eta2 = PYP(naf2,19)
        eta3 = PYP(naf1+1,19)
        xg_p = (shh+ sngl(Q2R))/yy/sngl(sss)
        xg_p = xpr
c     xg_p = 3.14159
        xm2_tot = -sngl(Q2R) + yy*xfgki*sngl(sss) + t2gki
        
        CALL PYROBO(0,0,STHETA,0.d0,0.d0,0.d0,0.d0)
        CALL PYROBO(0,0,0.d0,sphi,0.d0,0.d0,0.d0)
        CALL PYROBO(0,N,0.d0,0.d0,DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +       DBCMS(3)/DBCMS(4))
        
        if(pt2gen.ne.0.0) then
           pt2_gen = sngl(pt2gen)
           phi_gen = sngl(phigen)
        else
           pt2_gen=PT2H
           phi_gen=phitgki
        endif
        ncall = ncall + 1
c     CALL HFNT(100)
        
        NEVENT = NEVENT + 1
        NOUT=NEVENT
        
C     Lets generate a nice outputfile
        
        ievent=NOUTT
        genevent=NIN-lastgenevent

        write(txtLun,50) 0,ievent,genevent,ipro,idir,idisdif,avgi,sd,
     &       sss,q2,yy,xgam,xpr,p_t,pt2h,shh,t2gki,xfgki,shat,zqgki,
     &       xpgki,phitgki,N
 50     format((I4,1x,$),(I10,1x,$),4(I4,1x,$),16(f12.4,1x,$),I12,/)
        write(txtLun,*)'============================================'

        DO I=1,n
           write(txtLun,52) I,K(I,1),K(I,2),K(I,3),K(I,4),K(I,5),
     &          P(I,1),P(I,2),P(I,3),P(I,4),P(I,5),
     &          V(I,1),V(I,2),V(I,3)
        ENDDO
 52     format(2(I6,1x,$),I10,1x,$,3(I6,1x,$),5(f16.9,1x,$),
     &       3(f15.6,1x,$)/)
        write(txtLun,*)'=============== Event finished ==============='

        lastgenevent=NIN
      ENDIF
      RETURN

c...mode 3 is used to print the whole statistics information
3     CONTINUE
      CALL RAEND(20)
C      CALL PYSTAT(1)
C      CALL PYSTAT(4)
      CALL PTIME(' rapgap  ',2,0)
      CALL PTIME('        ',2,99)

      RETURN
      END

*=====dt_rgevntep=========================================================
*...Generate a rapgap event when using dpmjet
* Adapted from RAPGAP main program rgmain.F 
*=======================================================================
      SUBROUTINE DT_RGEVNTEP(Q2evt,Yevt,MODE,IREJ)
*     input:   MODE  1: generate a RAPGAP event get its Q2 and Y
*                    2: copy this event to dpmjet data common block
*     Output:  Q2evt, Yevt - event kinematics
*              IREJ  reject flag for current event (used in MODE 2)
      Implicit None

      double precision Q2evt, Yevt
      integer IREJ, MODE
      include 'rgrapgki.inc'  ! REAL YY etc.
      include 'rgpara1.inc'   ! DOUBLE PRECISION Q2
      include 'rgpart.inc'    ! DOUBLE SSS, PBEAM(2,5)...

C...Pythia Commons. Spell out to avoid IMPLICIT NONE clash w/ pythia.inc. 
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      INTEGER MINT,MSTP,MSTI
      DOUBLE PRECISION VINT,PARP,PARI

* added by liang to store the output event variables 1/20/12
      COMMON /EVTOUT/ XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT
      DOUBLE PRECISION XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT

      IF (MODE.EQ.1) THEN
         CALL PTIME(' event  ',1,0)
         CALL EVENT
C TEMPTEMP
         write(*,*) 'In DT_RGEVNTEP. List event:'
         CALL PYLIST(2)
         Q2evt = Q2
         Yevt = DBLE(YY)
         Q2OUT = Q2evt
         YYOUT = Yevt
         XBJOUT = XPR
         write(*,*) 'Rapgap Q2, y:',Q2,YY
         write(*,*) 'Rapgap s, xel, xpr',SSS,XEL,XPR
C Fill info. into Pythia common blocks?
         VINT(302)=SSS
         VINT(307)=Q2evt
         VINT(309)=Yevt
         write(*,*) 'PBEAM(1,*):',PBEAM(1,1),PBEAM(1,2),PBEAM(1,3),
     &        PBEAM(1,4),PBEAM(1,5)
         write(*,*) 'PBEAM(2,*):',PBEAM(2,1),PBEAM(2,2),PBEAM(2,3),
     &        PBEAM(2,4),PBEAM(2,5)
         write(*,*) 'XEL, XPR, PT2H, SHH, SHAT:',XEL,XPR,PT2H,SHH,SHAT
         PARI(33)=DBLE(XEL)
         PARI(34)=DBLE(XPR)
         PARI(18)=DBLE(PT2H)
         PARI(14)=DBLE(SHAT)
         VINT(4)=PBEAM(2,5)
         VINT(303)=PBEAM(1,5)
         W2OUT = Q2evt*(1.0d0/XBJOUT-1.0d0)+VINT(4)*VINT(4)
         NUOUT = Q2evt/(2.0d0*VINT(4)*XBJOUT)
         VINT(2) = W2OUT
         VINT(1) = SQRT(VINT(2))
         CALL PTIME(' event  ',2,0)
         IREJ=0
      ELSEIF (MODE.EQ.2) THEN
C  Swap particles 3 & 4 to match Pythia order: e p e' Z
         CALL RGTOPY
         CALL DT_PYEVNTEP(Q2evt,Yevt,MODE,IREJ)
      ELSE
         STOP ("DT_RGEVNTEP: FATAL ERROR. Called with invalid MODE")
      ENDIF

      RETURN
      END
