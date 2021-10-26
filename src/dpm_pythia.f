*These are the interfaces to pythia needed by dpmjet
*
* Written by Liang Zheng and Mark Baker
*
C...intialize pythia when using dpmjet
*======================================================================
      subroutine DT_PYINITEP(EPN,PPN,Q2MIN,Q2MAX,YMIN,YMAX,INPUT)
*     input:
*           EPN      electron beam momentum in lab frame
*           PPN      proton beam momentum in lab frame
*           Q2MIN    Q2 cut low
*           Q2MAX    Q2 cut high
*           YMIN     Y cut low
*           YMAX     Y cut high

      include 'pythia.inc'              ! All PYTHIA commons blocks
      include "mc_set.inc"
      include "py6strf.inc"
      include "mcRadCor.inc"
      include "radgen.inc"
      include "phiout.inc"
      include "beagle.inc"
      include "bea_pyqm.inc"

* properties of interacting particles
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5) 
c...target/proj mass, charge and projectile internal ID
      integer IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,ITMMOD,
     &     MODHYP,NHYPER,IDHYP 

      double precision EPN,PPN
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY10=1.0D-10,MAXNCL=260)
C...Added by liang 1/6/12
C...Switches for nuclear correction
      COMMON /PYNUCL/ INUMOD,CHANUM,ORDER,genShd
      SAVE /PYNUCL/
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER,genShd

C...Pythia event counter (since we keep PYINITing) Mark 2017-01-31
      COMMON /PYCNTR/ MYNGEN
      INTEGER MYNGEN
      SAVE /PYCNTR/

C...Pythia eA shadowing common block from Mark 2017-06-30
      COMMON /PYSHAD/ NKNOTS,RDUMMY,RAVAL,XKNOT(100),YKNOT(100),
     &BKNOT(100),CKNOT(100),DKNOT(100),EKNOT(100),FKNOT(100),SHDFAC
      SAVE /PYSHAD/
      DOUBLE PRECISION XKNOT,YKNOT,BKNOT,CKNOT,DKNOT,EKNOT,FKNOT,RAVAL,
     &SHDFAC
      INTEGER NKNOTS,RDUMMY
C... Locals by Mark 7/15/16
      INTEGER IKNOT
      DOUBLE PRECISION XKTMP, YKTMP, BKTMP, CKTMP, DKTMP, EKTMP, FKTMP
      DOUBLE PRECISION RATMP

C...Parameters and switch for energy loss
      DOUBLE PRECISION QHAT
      INTEGER QSWITCH
      COMMON /QUENCH/ QHAT, QSWITCH

C...output file name definition
      COMMON /OUNAME/ outname

      integer NEV, NPRT, ievent, genevent, I, tracknr 
      integer lastgenevent, idum1, idum2, initseed, nrtrack
      REAL trueX, trueW2, trueNu
      DOUBLE PRECISION sqrts, radgamE, radgamp, radgamEnucl
      DOUBLE PRECISION pbeamE, ebeamE, epznucl 
      DOUBLE PRECISION altpbeamE, altpbeam, altsqrts
      CHARACTER PARAM*100
      LOGICAL UseLut, GenLut

C MDB 2016-11-09 Reorganize to conserve charge & momentum for eA collisions 
C MDB 2017-05-25 Allow 4 cases according to ITMODE
C                n/p handling: 0 (D) = sequential n, then p      
C                              1 = all en collisions            
C                              2 = all ep collisions            
C                              3 = random mix (missing proc<91!)
C
C The nucleon and nuclear target rest frame should all be the same or
C our Pythia treatment of the hard collision will be inconsistent with
C our treatment of the nuclear remnant. So protons, neutrons, and A should 
C all have the same beta, gamma, rapidity, in the lab. So match the p/m.
C
C We are using the standard Pythia approach of treating the struck nucleon
C as being on mass shell (free nucleon mass).
C
C The input beam momentum: PPN=p_A/A
C For struck neutron: pbeamN = PPN*(A*Mn)/M_A
C For struck proton:  pbeamP = PPN*(A*Mp)/M_A
C
C Note: massp is nucleon mass, not necessarily proton
C
C MDB 2016-10-22 variables for nucleus A and it's nucleons P,N
C MAscl = M_A/A, 
C
C massp, masse, ebeam, pbeam, mcSet_EneBeam needed in mc_set.inc
C used in radgen routines in: pythia_radgen_extras.f, pythia_xsec.f,
C                           radgen_event.f, radgen.f, radgen_init.f
C PYTHIA routines: pydiff, pygaga, pysigh 
C and also: DT_PYEVNTEP, DT_PYOUTEP along with current DT_PYINITEP
C
      DOUBLE PRECISION Mprot,Mneut,Mnucl,Mlept
      ! beam type
      CHARACTER*10 tName
c ---------------------------------------------------------------------
c     Run parameter
c ---------------------------------------------------------------------
      integer*4 today(3), now(3)
c---------------------------------------------------------------------
c     ASCII output file and input file
c ---------------------------------------------------------------------
      CHARACTER*8 INPUT
      CHARACTER*256 outname

      integer LINP
      parameter ( LINP=28 )
      CHARACTER*256 inputfilename

* 14-Dec-2016 MDB Event variable to find datafiles
      CHARACTER*255 ENVDIR
      COMMON /ENVCOM/ ENVDIR
      CHARACTER*255 FILNAM

* Locals
      INTEGER IFSEED

c---------------------------------------------------------------------
! ... force block data modules to be read
C       external pydata
c ---------------------------------------------------------------------
       iModel=0
       pbeam=real(PPN)
       pbeamP=PPN
       pbeamN=PPN
       pbeamdbl=PPN
       ebeam=real(EPN) 
       ievent=0
       genevent=0
       lastgenevent=0
       tracknr=0
c ---------------------------------------------------------------------
c     Open ascii input file
c ---------------------------------------------------------------------
       inputfilename=INPUT
       open(LINP, file=inputfilename,STATUS='UNKNOWN')
       write(*,*) 'the input file is: ', inputfilename

C...Read output file name
       READ(LINP,*) outname
C...Read min/max x of radgen lookup table
       READ(LINP,*) mcSet_XMin, mcSet_XMax
C...Read information for cross section used in radgen
       READ(LINP,*) genSet_FStruct, genSet_R
C...Read parameters of radcorr: do radcorr (1), generate look-up table (2)
       READ(LINP,*) qedrad
       IF (qedrad.GT.0) WRITE(*,*) 'Warning: radcorr untested in BeAGLE'
C...Read parameters for PYTHIA-Model = which generation is done     
       READ(LINP,*) iModel
C...Read target type mass and charge
       READ(LINP,*) mcSet_TarA, mcSet_TarZ
       IF (mcSet_TarA.NE.IT .OR. mcSet_TarZ.NE.ITZ) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'WARNING: Pythia control file A,Z unused'
          WRITE(*,*) 'Beagle will use DPMJET A,Z: ',IT,', ',ITZ
       ENDIF
C...Read switch for shadowing
       READ(LINP,*) genShd
C...Read nuclear pdf correction order
       READ(LINP,*) ORDER
C...Read the switch for quenching
       READ(LINP,*) QSWITCH
C...Read q hat
       READ(LINP,*) QHAT
       IF (genShd.GT.1 .AND. QHAT.GT.1D-06) 
     &    STOP 'Nonzero qhat incompatible with genShd>1. Aborting!' 
C...Read information for cross section used in radgen
C...Mark allow fixed seed 08/14/2016 and shadowing factor 06/30/2017
       IFSEED=0
       SHDFAC=1.0
  100  READ(LINP,'(A)',END=200) PARAM
       IF (PARAM(1:5).EQ."FSEED") THEN
          READ(PARAM(7:16),'(I10)') IFSEED
          WRITE(*,*)'Using fixed PYR seed: ', IFSEED
       ELSEIF (PARAM(1:6).EQ."SHDFAC") THEN
          READ(PARAM(8:17),*) SHDFAC
          WRITE(*,*)'Using Shadowing Factor: ', SHDFAC 
       ELSE
          CALL PYGIVE(PARAM)
       ENDIF
       GOTO 100
c ---------------------------------------------------------------------
C...Initialize PYTHIA.      
c ---------------------------------------------------------------------
  200  CLOSE(LINP)
C...Added by Mark 2016-07-15. Hard-coded input file for now.
       IF (genShd.GE.2) THEN
          IF (IT.EQ.1) THEN
             WRITE(*,*)'ERROR: Shadowing corrections requested (',
     &                 'genShd>1) for A=1. genShd set to 1 (off).'
             genShd=1
          ELSE
             IF (IT.NE.197.OR.ITZ.NE.79) 
     &           WRITE(*,*)'WARNING: Shadowing corrections for gamma*',
     &             '+Au are being used. Not fully accurate for A=',IT
             CHECKR=0
             IKNOT=0
             IF (IT.EQ.197.AND.ITZ.EQ.79) THEN
                FILNAM = TRIM(ENVDIR)//'/shadowmapAu.dat'
             ELSEIF (IT.EQ.40.AND.ITZ.EQ.20) THEN
                FILNAM = TRIM(ENVDIR)//'/shadowmapCa.dat'
             ELSEIF(IT.LE.98) THEN
                FILNAM = TRIM(ENVDIR)//'/shadowmapCa.dat'
                WRITE(*,*)'WARNING: Shadowing corrections for gamma*',
     &               '+Cu are being used. Not fully accurate for A=',IT
             ELSE
                FILNAM = TRIM(ENVDIR)//'/shadowmapAu.dat'
                WRITE(*,*)'WARNING: Shadowing corrections for gamma*',
     &               '+Au are being used. Not fully accurate for A=',IT
             ENDIF
             WRITE(*,*)'Opening file: ',FILNAM
             OPEN(LINP, file=FILNAM,STATUS='OLD')
 300         READ(LINP,*,END=400)XKTMP,YKTMP,BKTMP,CKTMP,DKTMP,EKTMP,
     &                           FKTMP
             IKNOT=IKNOT+1
             XKNOT(IKNOT)=XKTMP
             YKNOT(IKNOT)=YKTMP
             BKNOT(IKNOT)=BKTMP
             CKNOT(IKNOT)=CKTMP
             DKNOT(IKNOT)=DKTMP
             EKNOT(IKNOT)=EKTMP
             FKNOT(IKNOT)=FKTMP
             GOTO 300
 400         CLOSE(LINP)
             NKNOTS=IKNOT
          ENDIF ! (IT.EQ.1)
       ENDIF ! (genShd.GT.2)
c...read parameters from dpmjet       
       INUMOD=IT
       CHANUM=ITZ
       mcSet_YMin=real(YMIN)
       mcSet_YMax=real(YMAX)
       mcSet_Q2Min=real(Q2MIN)
       mcSet_Q2Max=real(Q2MAX)

       write(*,*) '*********************************************'
       write(*,*) 'NOW all parameters are read by PYTHIA'
       write(*,*) '*********************************************'
       write(*,*) 'the output file is: ', outname
       if (IOULEV(4).GT.1) THEN
          call PYLIST(11)
          call PYLIST(12)
       endif
      
       print*,'kinematics cut read by PYTHIA:'
       print*,YMIN,' < y < ',YMAX,', ',Q2MIN,' < Q2 < ',Q2MAX

C...The random number for pythia 
C...Mark - fixed seed if requested in input for debugging 08/14/2016
      IF (IFSEED.NE.0) THEN
         initseed=IFSEED
      ELSE
         call getseed(initseed,.TRUE.)
      ENDIF
      write(6,*) 'SEED = ', initseed
      call rndmq (idum1,idum2,initseed,' ')
        
      !default ltype, lName and tName
      ltype = 11
      lName = 'gamma/e-'
      tName = 'p+'

      !set up lepton beam
      if( IJPROJ.eq.3 ) then
         lName = 'gamma/e-'
         ltype = 11
      elseif( IJPROJ.eq.4 ) then
         lName = 'gamma/e+'
         ltype = -11
      elseif( IJPROJ.eq.10 ) then
         lName = 'gamma/mu+'
         ltype = -13
      elseif( IJPROJ.eq.11 ) then
         lName = 'gamma/mu-'
         ltype = 13
      endif

      !set up nucleon beam
      Mneut=PYMASS(2112)
      Mprot=PYMASS(2212)
      Mlept=PYMASS(ltype)
      masse=real(Mlept)
C Mark 2015-10-21 This is for the first PYINIT only.
C                 For eA we may have to re-PYINIT event by event,
C                 but start with a neutron since that's the most likely
C Mark 2018-01-24 Use double precision here.
      if( IJTARG.eq.1 ) then 
         tName = 'p+'
         mcSet_TarA=1
         mcSet_TarZ=1
         massp=real(Mprot)
         Mnucl=Mprot
         idNucPY=2212
         idNucBAM=1
      elseif ( IJTARG.eq.8 ) then
         tName = 'n0'
         mcSet_TarA=1
         mcSet_TarZ=0
         massp=real(Mneut)
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
         pbeamdbl=PPN*Mnucl/MAscl
         pbeam=real(pbeamdbl)
         pbeamP=(PPN*Mprot/MAscl)
         pbeamN=(PPN*Mneut/MAscl)
      ENDIF
C     proton (neutron) is defined in positive z and as target
      P(2,1)=0.0  
      P(2,2)=0.0
      P(2,3)= pbeamdbl
      PZTARG = PPN
      PZNUCL = pbeamdbl
C     lepton is defined in negative z and as beam
      P(1,1)=0.0  
      P(1,2)=0.0  
      P(1,3)=-EPN
      PZLEP = -EPN

      pbeamE=sqrt(pbeamdbl**2+Mnucl**2)
      pbeta=pbeamdbl/pbeamE
      pgamma=pbeamE/Mnucl
      ebeamE=sqrt(EPN**2+Mlept**2)
      ebeamEnucl=pgamma*ebeamE-pgamma*pbeta*(-EPN)
      epznucl=-pgamma*pbeta*(ebeamE)+pgamma*(-EPN)
      write(*,*) ebeamEnucl, ebeamE, epznucl, -EPN
      mcSet_EneBeam=sngl(ebeamEnucl)

      sqrts=sqrt((pbeamE+ebeamE)**2-(pbeamdbl-EPN)**2)
      altpbeam=PPN
      altpbeamE=sqrt(PPN*PPN+PYMASS(2212)**2)
      altsqrts=sqrt((altpbeamE+ebeamE)**2-(altpbeam-EPN)**2)
      write(*,*) '*********************************************'
      IF(IT.GT.1) THEN
         write(*,*) 'Nucleus: Z, A: ', ITZ, ' ', IT
         write(*,*) 'Per nucleon beam momentum:', PPN, 'GeV'
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
      write(*,*) 'lepton beam momentum: ', EPN, 'GeV'
      write(*,*) 'sqrts for eN system: ', sqrts, 'GeV'
      write(*,*) 'lepton beam energy in TRF: ',ebeamEnucl,'GeV'
      write(*,*) '*********************************************'

      if (iModel.eq.0) then
         UseLUT=.false.
         GenLUT=.false.
         qedrad=0
         MSTP(199)=0
         mcRadCor_EBrems=0.
      elseif (iModel.eq.1) then
         if (qedrad.eq.0) then
            mcRadCor_EBrems=0.
            UseLUT=.false.
            GenLUT=.false.
            MSTP(199)=1
         elseif (qedrad.eq.1) then
            mcRadCor_EBrems=0.
            UseLUT=.true.
            GenLUT=.false.
            MSTP(199)=1
            call radgen_init(UseLUT,GenLUT)
            write(*,*) 'I have initialized radgen'
         elseif (qedrad.eq.2) then
            write(*,*) 'radgen lookup table will be generated'
            mcRadCor_EBrems=0.
            UseLUT=.true.
            GenLUT=.true.
            MSTP(199)=1
            call radgen_init(UseLUT,GenLUT)
            goto 500
         endif
      endif

C     GenNucDens is used even without quenching.
      call GenNucDens(ITZ, IT)
      IF(QSWITCH.EQ.1) THEN
         print*,'Quenching requested, qhat=',QHAT
         print*,'Using recoil factor:',PQRECF
         print*,'      Ptmodel iPtF: ',PYQ_IPTF
         print*,'      EmitGluon iEg:',PYQ_IEG
         print*,'      SupFactor:    ',PYQ_SUPF
         print*,'      Heavy Quarks: ',PYQ_HQ
         print*,'      Energy Thres: ',PYQ_IET

c...when quenching is used switch off internal parton shower
C... 2020-03-24 MDB - No don't!
C         MSTP(61)=0
C         MSTP(71)=0
      ENDIF

      !changed by liang to fit with INC particle status
      CALL DT_PYDECY(1)

      !initialize pythia
      MYNGEN=0
C-TEMP-TEMP-TEMP
C      write(*,*) '1st call to pyinit, MYNGEN=',MYNGEN
      call pyinit('3MOM', lName, tName, WIN)
  500  if (qedrad.eq.2) then
         write(*,*) 'lookup table is generated;'
         write(*,*) 'to run now pythia change parameter qedrad to 1'
       endif

      RETURN
      END
*=====dt_pydecy=====================================================
*used to setup the particle decay channels in pythia
*since we don't want E-W decay particles before intra-nuclear
*cascade process (those particles decay way out of nuclear
*environment), we will switch them off in pythia initialzation
*but turn on for decay after all the particle generation is done
* input: MODE 1: switch off decay for the list particles
*             2: Restore original decay state for the list particles
* MDB 2019-09-11 Major change. Mode 2 used to switch on decay, now it
*     restores the status quo ante. Also the routine will printout 
*     the initial state of the MDCY array.
*
      SUBROUTINE DT_PYDECY(MODE)

      include 'pythia.inc'              ! All PYTHIA commons blocks

      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
      SAVE LFIRST
      CHARACTER PNAME*16

*c... Extended to allow jpsi, phi to decay outside nuclei
      DIMENSION IDXSTA(40),MDCYINI(40)
      DATA IDXSTA
*          K0s   pi0  lam   alam  sig+  asig+ sig-  asig- tet0  atet0
     &  /  310,  111, 3122,-3122, 3222,-3222, 3112,-3112, 3322,-3322,
*          tet- atet-  om-  aom-   D+    D-    D0    aD0   Ds+   aDs+
     &    3312,-3312, 3334,-3334,  411, -411,  421, -421,  431, -431,
*          etac lamc+alamc+sigc++ sigc+ sigc0asigc++asigc+asigc0 Ksic+
     &     441, 4122,-4122, 4222, 4212, 4112,-4222,-4212,-4112, 4232,
*         Ksic0 aKsic+aKsic0 sig0 asig0 jpsi phi
     &    4132,-4232,-4132, 3212,-3212, 443, 333, 3*0/

      INTEGER MODE

c...Store & printout initial decay status on first call. 
      IF (LFIRST) THEN
         WRITE(*,*)
         WRITE(*,*) "WEAK DECAY TABLE"
         WRITE(*,*) "Particle, PID, MDCY value"
         DO I=1,37
            MDCYINI(I) = MDCY(PYCOMP(IDXSTA(I)),1)
            CALL PYNAME(IDXSTA(I),PNAME)
            WRITE(*,*) PNAME, IDXSTA(I), MDCYINI(I)
         ENDDO
         WRITE(*,*)
         LFIRST=.FALSE.
      ENDIF

      GOTO(1,2) MODE

c...Switch off decay 
1     CONTINUE

      DO I=1,37
         MDCY(PYCOMP(IDXSTA(I)),1) = 0
      ENDDO

      RETURN

c...Restore initial decay settings
2     CONTINUE

      DO I=1,37
C         MDCY(PYCOMP(IDXSTA(I)),1) = 1
         MDCY(PYCOMP(IDXSTA(I)),1) = MDCYINI(I)
      ENDDO

      RETURN

      END

*=====dt_decall=====================================================
*used to perform all long time decay in pythia at the end of run
*liang-2021-10 
*
      SUBROUTINE DT_DECALL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
* event history

      PARAMETER (NMXHKK=200000)

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHIST(2,NMXHKK)

      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      PARAMETER (MAXLND=4000)
      COMMON/PYJETS/N,NPAD,K(MAXLND,5),P(MAXLND,5),V(MAXLND,5)

      INTEGER PYCOMP,PYK

      DIMENSION IHISMO(NMXHKK),P1(4)

*c...data array added by liang to do test whether 1 particle is supposed
*c...decay. Extended to allow jpsi, phi to decay outside nuclei
      DIMENSION IDXSTA(40)
      DATA IDXSTA
*          K0s   pi0  lam   alam  sig+  asig+ sig-  asig- tet0  atet0
     &  /  310,  111, 3122,-3122, 3222,-3222, 3112,-3112, 3322,-3322,
*          tet- atet-  om-  aom-   D+    D-    D0    aD0   Ds+   aDs+
     &    3312,-3312, 3334,-3334,  411, -411,  421, -421,  431, -431,
*          etac lamc+alamc+sigc++ sigc+ sigc0asigc++asigc+asigc0 Ksic+
     &     441, 4122,-4122, 4222, 4212, 4112,-4222,-4212,-4112, 4232,
*         Ksic0 aKsic+aKsic0 sig0 asig0 jpsi phi
     &    4132,-4232,-4132, 3212,-3212, 443, 333, 3*0/

      !allow long life-time particle decay
      CALL DT_PYDECY(2)

      NN  = 0
      DO 1 I=1,NHKK
        ISTAB = 1 ! status to label a stable particle, 1-stable,0-decay
        DO J=1,37
           IF(IDHKK(I).EQ.IDXSTA(J)) ISTAB=0
        ENDDO

        IF ( (ISTHKK(I).EQ.1).AND.ISTAB.EQ.0 ) THEN
           NN = NN+1

           !liang-test-2021-10
c           print*,'id=',IDHKK(I),' px=',PHKK(1,I),' py=',PHKK(2,I),
c     &  ' pz=',PHKK(3,I),' x=',VHKK(1,I),' y=',VHKK(2,I),' z=',VHKK(3,I)

           !copy wanted decay hadrons to PYJETS
           K(NN,1)=2
           K(NN,2)=IDHKK(I)
           K(NN,3)=0
           P(NN,5)=PHKK(5,I)
           DO JDX=1,4
           P(NN,JDX)=PHKK(JDX,I)
           V(NN,JDX)=VHKK(JDX,I)
           ENDDO
          
           IHISMO(NN)= I

        ENDIF

    1 CONTINUE

      IF (NN.GT.0) THEN
      !found NN number of particles to be decayed
         N=NN !need to set the PYJETS size
         CALL PYEXEC
         !liang-test-2021-10
c         call pylist(3)

         NHKOLD=NHKK !record DTEVT size before add decay

         !start from the first decay daughter
         DO 2 II=1,N
           IF(II.le.NN) THEN
           !decay mothers
             IDX=IHISMO(II)
             ISTHKK(IDX)=2
             !NN must be subtracted converting PYJETS index to 
             !DTEVT line index
             JDAHKK(1,IDX)=K(II,4)-NN+NHKOLD
             JDAHKK(2,IDX)=K(II,5)-NN+NHKOLD
           ELSE
           !decay daughters
             NHKK=NHKK+1
             PHKK(5,NHKK)=P(II,5)
             IDHKK(NHKK)=K(II,2)
             ISTHKK(NHKK)=K(II,1)
             !unstable daughter
             IF(ISTHKK(NHKK).ne.1) THEN
               ISTHKK(NHKK)=2
               JDAHKK(1,NHKK)=K(II,4)-NN+NHKOLD
               JDAHKK(2,NHKK)=K(II,5)-NN+NHKOLD
             ENDIF
             IF(K(II,3).le.NN) THEN
               !primary decay
               JMOHKK(1,NHKK)=IHISMO(K(II,3))
             ELSE
               !secondary decay
               JMOHKK(1,NHKK)=K(II,3)-NN+NHKOLD
             ENDIF
             DO JDX=1,4
               PHKK(JDX,NHKK)=P(II,JDX)
               VHKK(JDX,NHKK)=V(II,JDX)
             ENDDO
             !inherit mother nobam
             !as mother-daughter already sorted in JMOHKK, can easily
             !use the JMOHKK nobam value
             NOBAM(NHKK)=NOBAM(JMOHKK(1,NHKK))
           ENDIF
    2    CONTINUE

      ENDIF

      !get back to INC particle status
      CALL DT_PYDECY(1)


      RETURN
      END


*=====dt_pyevnt========================================================
*used to the sample a pythia event and dump this event to DTEVT1
*Some important thing to be mentioned here is the output of pyhia
*data is in the lab frame, while the calculation in dpmjet is made
*in photon-proton c.m.s frame. So when we are copying the data common
*to dpmjet, we must make the right Lorentz transformation. Plus,
*the virtual photon has to be directed to the z+ direction. A rotation
*for the reference is also necessary.

      SUBROUTINE DT_PYEVNTEP(Q2,YY,MODE,IREJ)
 
*     input:
*           MODE  1: generate a pythia event get its Q2 and Y
*                 2: copy this event to dpmjet data common block
*     output:
*           Q2    Q2 of this current event (used in MODE 1)
*           YY    Y of this current event  (used in MODE 1)
*           IREJ  reject flag for current event (used in MODE 2)

      include 'pythia.inc'              ! All PYTHIA commons blocks
      include "mc_set.inc"
      include "py6strf.inc"
      include "mcRadCor.inc"
      include "radgen.inc"
      include "phiout.inc"
      include "beagle.inc"
      include "bea_pyqm.inc"

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 

* event history

      PARAMETER (NMXHKK=200000)

      PARAMETER (FM2MM=1.0D-12)

* properties of interacting particles
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5) 
c...target/proj mass, charge and projectile internal ID
      integer IT, ITZ, IP, IPZ, IJPROJ, IBPROJ, IJTARG, IBTARG, ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP 

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history - IHIST renamed because it's in pythia.inc
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHISTDPM(2,NMXHKK)

      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY10=1.0D-10,MAXNCL=260)

C...Parameters and switch for energy loss from PyQM
      DOUBLE PRECISION QHAT
      INTEGER QSWITCH
      COMMON /QUENCH/ QHAT, QSWITCH

* Lorentz-parameters of the current interaction from DPMJET
      COMMON /DTLTRA/ GACMS(2),BGCMS(2),GALAB,BGLAB,BLAB,
     &                UMO,PPCM,EPROJ,PPROJ

* lorentz transformation parameter 
      COMMON /LTPARA/ BGTA(4), GAMM(5), eveBETA(4), GAA
      COMMON /ROTATE/ COF,COD,SIF,SID

* properties of photon/lepton projectiles from DPMJET
      COMMON /DTGPRO/ VIRT,PGAMM(4),PLEPT0(4),PLEPT1(4),PNUCL(4),IDIREC

* kinematics at lepton-gamma vertex from DPMJET
      COMMON /DTLGVX/ PPL0(4),PPL1(4),PPG(4),PPA(4)

C* position of interacted nucleon   ! part of bea_pyqm.inc
C      DOUBLE PRECISION PosNuc
C      COMMON /NPARINT/ PosNuc(4)
* Added by Mark 2016-08-18
      INTEGER MXINTS,IIGA
      PARAMETER (MXINTS=300)
      INTEGER IINTER(MXINTS),IPARTN,IMAXZ,NPOS(MXINTS),IMAIN
      DOUBLE PRECISION PHIGAM,THEGAM,PTGAM,PTK1,PTK2,PZMAX
      DOUBLE PRECISION PosAlt(MXINTS,4)
      INTEGER MomAlt(MXINTS)
      INTEGER IDAlt(MXINTS)
      DOUBLE PRECISION ZMIN, ZTRY
      DOUBLE PRECISION EEFIX,PZEFIX,BETAFX,GAMMFX

Cc...added by liang & Mark to include pythia energy loss datas
      double precision PAUX, DPF
      COMMON /PFAUX/ PAUX(4), DPF(4)

      DOUBLE PRECISION P1, P2, P3, P4, PTOT, PP1, PP2, PP3, PP4, QCHG

      DOUBLE PRECISION UPVP,DNVP,USEAP,DSEAP,STRP,CHMP,BOTP,TOPP,GLUP
      DOUBLE PRECISION UPVA,DNVA,USEAA,DSEAA,STRA,CHMA,BOTA,TOPA,GLUA

C      DOUBLE PRECISION USEAP1000, USEAP1001, USEAP1690

      INTEGER MODE, IREJ

      DOUBLE PRECISION Q2, YY, XX, XXALT

C...Shadowing Map & Switches for nuclear correction - Mark 7/22/16
C...Added by liang 1/6/12
C...Switches for nuclear correction
      COMMON /PYNUCL/ INUMOD,CHANUM,ORDER,genShd
      SAVE /PYNUCL/
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER,genShd

      EXTERNAL ANEAR
      DOUBLE PRECISION ATEMP
      CHARACTER*8 chpset
      CHARACTER*20 Nprm

C...Pythia event counter (since we keep PYINITing) Mark 2017-01-31
      COMMON /PYCNTR/ MYNGEN
      INTEGER MYNGEN
      SAVE /PYCNTR/

C...Pythia eA shadowing common block from Mark 2017-06-30
      COMMON /PYSHAD/ NKNOTS,RDUMMY,RAVAL,XKNOT(100),YKNOT(100),
     &BKNOT(100),CKNOT(100),DKNOT(100),EKNOT(100),FKNOT(100),SHDFAC
      SAVE /PYSHAD/
      DOUBLE PRECISION XKNOT,YKNOT,BKNOT,CKNOT,DKNOT,EKNOT,FKNOT,RAVAL,
     &SHDFAC
      INTEGER NKNOTS,RDUMMY

C  Local

      LOGICAL LFIRST
      INTEGER IIMAIN, IIMAINN
      DOUBLE PRECISION PDEUT(5), PSPEC(5), MNUCL

      SAVE LFIRST

      GOTO(1,2) MODE

c...generate a pythia event
1     CONTINUE      
c...MDB 2016-11-11 For eA, decide whether we are hitting a proton or neutron
c...               and reinitialize Pythia (and radgen) if needed
c...               Note: we assume p,n equally likely. Not strictly correct
c...               especially at low energies. 
c...
      if (ITMODE.eq.3) then
         if (NINT(INUMOD).GT.1) then
C     WRITE(*,*)'Prob(proton) = Z/A = ',CHANUM/INUMOD
            if (DT_RNDM(RDUMMY).LE.CHANUM/INUMOD) then
c...  proton - reinit if it was a neutron last event
C     WRITE(*,*)'Choice for this event: proton'
               if(idNucPY.NE.2212) then
CC     WRITE(*,*)'Need to reinitialize for proton'
                  CALL REINIT(2212)
CC     -TEMP-TEMP-TEMP
CC     write(*,*) 'pre-pyinit, MYNGEN=',MYNGEN
C                  call pyinit('3MOM', lName, 'p+', WIN)
               endif
            else
c...  neutron - reinit if it was a proton last event
C     WRITE(*,*)'Choice for this event: neutron'
               if(idNucPY.NE.2112) then
CC     WRITE(*,*)'Need to reinitialize for neutron'
                  CALL REINIT(2112)
CC     -TEMP-TEMP-TEMP
CC     write(*,*) 'pre-pyinit, MYNGEN=',MYNGEN
C                  call pyinit('3MOM', lName, 'n0', WIN)
               endif
            endif 
         endif
      endif

C MDB - Save fragmentation & particle decay for later 
      MSTJ(1) = 0
      MSTJ(21) = 0
C      print*,'in event:', NEVENT,' MYNGEN=',MYNGEN,
C     &       ' calling PYEVNT'
C     Don't decay particles yet if we are going to apply pF or ptkick
999   CALL PYEVNT
      MYNGEN = MYNGEN+1
C-TEMP-TEMP-TEMP
C      print*,'After PYEVNT: MYNGEN=',MYNGEN
      IF(MSTI(61).EQ.1) THEN
         WRITE(*,*) 'go back to PYEVNT call'
         GOTO 999
      ENDIF

      Q2=VINT(307)
      YY=VINT(309)
      XX = Q2/YY/(VINT(302)-VINT(4)**2-VINT(303)**2)

      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
         print*,'in event:', NEVENT
         print*,'Q2 from pythia:',Q2
         print*,'Y from pythia:',YY
         print*,'X Bj from pythia:',XX
C         print*,'X Bj alt Q2/(W2+Q2-M2):',Q2/(VINT(2)+Q2-VINT(4)**2)
         print*,'W2 from VINT:',VINT(2)
         print*,'M from VINT(4):',VINT(4)
         print*,'me2 from VINT:',VINT(303)**2
         print*,'s from VINT:',VINT(302)
         print*,'subprocess: ',MINT(1)
C         CALL PYLIST(2)
      endif

      RAVAL = 1.0
      IF (genShd.GE.2) THEN
C Simplified method. Use sea RA(xBj, Q2) for all processes
C x = Q2/y/(s-M2-m2) matches x from PYPDFU better than  Q2/(W2+Q2-M2)

         IF(XX.LE.0.1) THEN
C            CALL STRUCTM(XX,1.00,UPVP,DNVP,USEAP,DSEAP,STRP,CHMP,BOTP
C     &                  ,TOPP,GLUP)
C            WRITE(*,*)"UseaP(X=",XX,",Q2=1.00 GeV2 from STRUCTM: ",USEAP
C            USEAP1000=USEAP
C            CALL STRUCTM(XX,1.001,UPVP,DNVP,USEAP,DSEAP,STRP,CHMP,BOTP
C     &                  ,TOPP,GLUP)
C            WRITE(*,*)"UseaP(X=",XX,",Q2=1.001 GeV2 from STRUCTM: ",
C     &                   USEAP
C            USEAP1001=USEAP
C            CALL STRUCTM(XX,1.69,UPVP,DNVP,USEAP,DSEAP,STRP,CHMP,BOTP
C     &                  ,TOPP,GLUP)
C            WRITE(*,*)"UseaP(X=",XX,",Q2=1.69 GeV2 from STRUCTM: ",USEAP
C            USEAP1690=USEAP
            CALL STRUCTM(XX,Q2,UPVP,DNVP,USEAP,DSEAP,STRP,CHMP,BOTP,TOPP
     &                  ,GLUP)
C            WRITE(*,*)"UseaP(X=",XX,",Q2=",Q2," GeV2 from STRUCTM: ",
C     &                   USEAP
            ATEMP = INUMOD
            IF((CHANUM.GT.1D0).AND.(ORDER/100.GT.0)) THEN
               write(chpset,'(I2)') MOD(ORDER,100)
               IF(ORDER/100.EQ.1) THEN
C                 print*,'Now run LO nuclear pdf'
                  Nprm='EPS09LO,'//chpset
                  CALL SETLHAPARM(Nprm)
C     MDB - Use nearest A value as advised by H. Paukkunen 10/3/17
                  ATEMP = ANEAR(INUMOD)
               ELSE IF(ORDER/100.EQ.2) THEN
C                 print*,'Now run NLO nuclear pdf'
                  Nprm='EPS09NLO,'//chpset
                  CALL SETLHAPARM(Nprm)
               ENDIF
            ENDIF
            CALL STRUCTA(XX,Q2,ATEMP,UPVA,DNVA,USEAA,DSEAA,STRA,CHMA,
     &                   BOTA,TOPA,GLUA)
C
C    Simplified eA shadowing implementation - Mark 2016-08-18   
C    Use  u-sea R(A) for xBj for all processes including diffraction
C   and 95 etc.  For R(A)>1 or x>0.1, set/leave R(A)=1
C Note: For PGF (and QCDC) the x_SF is not the same as x_Bj
C       Proper flavor dependence will take a bit of care...
            RAVAL = USEAA/USEAP
C            print*,'Used: ',Nprm
C            print*,'Raw RA value: ',RAVAL
            IF (RAVAL.GE.1.D0) RAVAL=1.D0
         ENDIF
C         print*,'Final RA value: ',RAVAL
C         print*
      ENDIF

C     TEMP-TEMP-TEMP - comment these in instead of out
C      WRITE(*,*) 'End of DT_PYEVNTEP mode 1, PYLIST(2):'
C      CALL PYLIST(2)
      LFIRST=.FALSE.

      RETURN

c...transform this event to certain reference and copy it to DPMJET
2     CONTINUE
c...guess the position of final state particle
c...NPOINT(4), the position where final particle starts is very important!!    

      IREJ=0

      IF (IOULEV(4).GE.2) THEN
         WRITE(*,*) 'Start of DT_PYEVNTEP(mode=2)'
         CALL PYLIST(2) 
      ENDIF

C     MDB 2017-02-22
C     DT_FOZOCA, DT_SCN4BA need NPOINT(1) pointing to the last nucleon
C     DT_FOZOCA, DT_RESNCL, DT_SCN4BA: NPOINT(4) points to the 1st produced particle
C     DT_CHASTA: NPOINT(3) points to the 1st produced particle
C                                   
      NPOINT(1)=NHKK
      NPOINT(2)=NPOINT(1)
      NPOINT(3)=NPOINT(1)+5
      NPOINT(4)=NPOINT(1)+5

********find the position of wounded nucleon
c...Don't assume just one! Mark 2016-08-18
      NINTS=0
      DO J=1,MXINTS
         IMAIN=0
         IINTER(J)=0
         NPOS(J)=0
         DO KK=1,4
            PosAlt(J,KK)=0.0D0
         ENDDO
      ENDDO
      DO J=1,NHKK      
         IF(ISTHKK(J).EQ.12) THEN
            NINTS=NINTS+1
            IINTER(NINTS)=J
         ENDIF
      ENDDO

C...  Pick main (Pythia) interaction
      IF (NINTS.LT.1) THEN
         WRITE(*,*) 'DT_PYEVNTEP ERROR: No Interactions!'
         STOP
      ELSEIF (genShd.LE.1) THEN
C...  For genShd=1 should only be one.
         IMAIN = 1
         IF (NINTS.GT.1) THEN
            write(*,*) 'DT_PYEVNTEP ERROR: Found ',NINTS,
     &                 ' interactions. Expected 1.'
         ENDIF
      ELSEIF (genShd.EQ.2) THEN
C...  For genShd=2, pick the main Pythia interaction at random
         IMAIN=INT(NINTS*DT_RNDM(RDUMMY))+1
      ELSEIF (genShd.GE.3) Then
C...  The first interaction is the main Pythia interaction
C...  Note: VHKK is always in TRF
         ZMIN=10000.0E0
         IMAIN=0
         DO III=1,NINTS       
            ZTRY=VHKK(3,IINTER(III))
            IF (ZTRY.LT.ZMIN) THEN
               ZMIN=ZTRY
               IMAIN=III
            ENDIF
         ENDDO
      ENDIF
      IIMAIN = IINTER(IMAIN)
      PosNuc(1)=VHKK(1,IIMAIN)
      PosNuc(2)=VHKK(2,IIMAIN)
      PosNuc(3)=VHKK(3,IIMAIN)
      PosNuc(4)=VHKK(4,IIMAIN)
      PXF = PHKK(1,IIMAIN)
      PYF = PHKK(2,IIMAIN)
      PZF = PHKK(3,IIMAIN)
      EKF = PHKK(4,IIMAIN)-PHKK(5,IIMAIN)
C Thickness is twice the integral from 0 to infinity
      THKB = 2.0D0*DCALC(0.0D0)
      IDUM = 0
      THKSCL = THKB*SCLFAC(IDUM)
C     Calculate distance travelled in the nucleus
      IF (NINTS.EQ.1) THEN
         DAVG = DCALC(VHKK(3,IIMAIN)/FM2MM)
         DFIRST = DAVG
      ELSE
         ZMIN = 10000.0
         DO III=1,NINTS       
            ZTRY=VHKK(3,IINTER(III))
            IF (ZTRY.LT.ZMIN) THEN
               ZMIN=ZTRY
            ENDIF
         ENDDO
         DFIRST = DCALC(ZMIN/FM2MM)
         DAVG = DCALC(VHKK(3,IIMAIN)/FM2MM)
      ENDIF
c     WRITE(*,*) 'Calculated D1st, Davg: ',DFIRST,' ',DAVG

c      write(*,*) 'DPMJET position',PosNuc(1),PosNuc(2),PosNuc(3)
      PYQREC(1)=0.0D0
      PYQREC(2)=0.0D0
      PYQREC(3)=0.0D0
      PYQREC(4)=0.0D0
      IF (USERSET.EQ.1 .OR. USERSET.EQ.2) THEN
         USER1=0.0D0
         USER2=0.0D0
         USER3=0.0D0
      ENDIF

c...MDB 2016-11-14 Make sure the nucleon is the correct flavor
c...               If not, swap flavors with a random nucleon
c...MDB 2017-01-21 Swap momenta as well as flavor - neutrons have a
c...               higher Fermi momentum than protons for A/Z>2.
c...               Maybe this is destabilitzing the nucleus.
C      WRITE(*,*) 'Output current event before possible swap'
C      CALL DT_PYOUTEP(4)
      IF (IDHKK(IIMAIN).NE.idNucPY) THEN
C         WRITE(*,*) 'SWAP! Interaction# is: ',IIMAIN,' ID=',
C     &    IDHKK(IIMAIN),' idNucPY=',idNucPY,' idNucBAM=',idNucBAM
c... Pick a # between 2 and NHKK and then look for a nucleon to swap flavors
         ISWAP=INT((NHKK-1)*DT_RNDM(RDUMMY))+2
         ISWTRY=0
         DO WHILE (ISWTRY.LT.NHKK-1)            
            IF (IDHKK(ISWAP).EQ.idNucPY) GOTO 25
            ISWTRY=ISWTRY+1
            ISWAP = ISWAP+1
            IF (ISWAP.GT.NHKK) ISWAP=2
         ENDDO
         WRITE(*,*)"ERROR: Can't find a nucleon to swap!"
         ISWAP=IIMAIN
 25      CONTINUE
C         WRITE(*,*)'Swap nucleon# is: ', ISWAP,'ID= ',IDHKK(ISWAP)
C         WRITE(*,*)'ISTHKK, , PHKK(1-5), VHKK(1-3), IDRES, IDXRES,' 
C     &          ,'NOBAM, IDBAM, IDCH'     
C         WRITE(*,*)'BEFORE SWAP:'
C         WRITE(*,*)'Main: ',ISTHKK(IIMAIN),IDHKK(IIMAIN), 
C     &      PHKK(1,IIMAIN), PHKK(2,IIMAIN), PHKK(3,IIMAIN), 
C     &      PHKK(4,IIMAIN), PHKK(5,IIMAIN), VHKK(1,IIMAIN), 
C     &      VHKK(2,IIMAIN), VHKK(3,IIMAIN), IDRES(IIMAIN), 
C     &      IDXRES(IIMAIN), NOBAM(IIMAIN), IDBAM(IIMAIN), IDCH(IIMAIN)
C         WRITE(*,*)'Swap: ',ISTHKK(ISWAP),IDHKK(ISWAP), PHKK(1,ISWAP),
C     &      PHKK(2,ISWAP), PHKK(3,ISWAP), PHKK(4,ISWAP), 
C     &      PHKK(5,ISWAP), VHKK(1,ISWAP), VHKK(2,ISWAP),
C     &      VHKK(3,ISWAP), IDRES(ISWAP), IDXRES(ISWAP),
C     &      NOBAM(ISWAP), IDBAM(ISWAP), IDCH(ISWAP)
         IF (ISWAP.NE.IIMAIN) THEN
            IDHKK(ISWAP)=IDHKK(IIMAIN)
            IDHKK(IIMAIN)=idNucPY
            DO I5=1,5 
               DTEMP = PHKK(I5,IIMAIN)
               PHKK(I5,IIMAIN)=PHKK(I5,ISWAP)
               PHKK(I5,ISWAP) = DTEMP
            ENDDO
            ITEMP = NOBAM(IIMAIN)
            NOBAM(IIMAIN)=NOBAM(ISWAP)
            NOBAM(ISWAP)=ITEMP
            ITEMP = IDBAM(IIMAIN)
            IDBAM(IIMAIN)=IDBAM(ISWAP)
            IDBAM(ISWAP)=ITEMP
            ITEMP = IDCH(IIMAIN)
            IDCH(IIMAIN)=IDCH(ISWAP)
            IDCH(ISWAP)=ITEMP
C     Use new nucleon to specify PF
            PXF = PHKK(1,IIMAIN)
            PYF = PHKK(2,IIMAIN)
            PZF = PHKK(3,IIMAIN)
            EKF = PHKK(4,IIMAIN)-PHKK(5,IIMAIN)
         ENDIF
C         WRITE(*,*)'AFTER SWAP:'
C         WRITE(*,*)'Main: ',ISTHKK(IIMAIN),IDHKK(IIMAIN), 
C     &      PHKK(1,IIMAIN), PHKK(2,IIMAIN), PHKK(3,IIMAIN), 
C     &      PHKK(4,IIMAIN), PHKK(5,IIMAIN), VHKK(1,IIMAIN), 
C     &      VHKK(2,IIMAIN), VHKK(3,IIMAIN), IDRES(IIMAIN), 
C     &      IDXRES(IIMAIN), NOBAM(IIMAIN), IDBAM(IIMAIN), IDCH(IIMAIN)
C         WRITE(*,*)'Swap: ',ISTHKK(ISWAP),IDHKK(ISWAP), PHKK(1,ISWAP),
C     &      PHKK(2,ISWAP), PHKK(3,ISWAP), PHKK(4,ISWAP), 
C     &      PHKK(5,ISWAP), VHKK(1,ISWAP), VHKK(2,ISWAP),
C     &      VHKK(3,ISWAP), IDRES(ISWAP), IDXRES(ISWAP),
C     &      NOBAM(ISWAP), IDBAM(ISWAP), IDCH(ISWAP)
c      ELSE
C         WRITE(*,*) 'KEEP! nucleon# is: ',IIMAIN,' ID=',
C     &    IDHKK(IIMAIN),' idNucPY=',idNucPY,' idNucBAM=',idNucBAM
C         WRITE(*,*)'ISTHKK, , PHKK(1-5), VHKK(1-3), IDRES, IDXRES,' 
C     &          ,'NOBAM, IDBAM, IDCH'     
C         WRITE(*,*)'Main: ',ISTHKK(IIMAIN),IDHKK(IIMAIN), 
C     &      PHKK(1,IIMAIN), PHKK(2,IIMAIN), PHKK(3,IIMAIN), 
C     &      PHKK(4,IIMAIN), PHKK(5,IIMAIN), VHKK(1,IIMAIN), 
C     &      VHKK(2,IIMAIN), VHKK(3,IIMAIN), IDRES(IIMAIN), 
C     &      IDXRES(IIMAIN), NOBAM(IIMAIN), IDBAM(IIMAIN), IDCH(IIMAIN)
      ENDIF
C     Calculate kinematics of actual eN collision w/o Fermi motion
CC      MAscl = AZMASS(NINT(INUMOD),NINT(CHANUM))/INUMOD ! Done already
C      WRITE(*,*) 'ebeamEnucl, INUMOD, CHANUM, AZMASS, AZMASS/A:',
C     &            ebeamEnucl, ' ', INUMOD, ' ', CHANUM, ' ', 
C     &            AZMASS(NINT(INUMOD),NINT(CHANUM)),' ',MAscl
C      WRITE(*,*) 'sqrts(e+ A/A)  = ',sqrt(MAscl*MAscl+masse*masse+
C     &                                    2.0*MAscl*ebeamEnucl)
C      WRITE(*,*) 'masse = ',masse
C      WRITE(*,*) 'Naive sqrtss(eN)= ',sqrt(
C     &     ((2.0*ebeamEnucl+PHKK(5,IIMAIN))*PHKK(5,IIMAIN))+masse*masse)

C      WRITE(*,*) 'Output current event after possible swap'
C      CALL DT_PYOUTEP(4)

      CALL DT_PICKSRC(PHKK,VHKK,NINT(INUMOD),IIMAIN,IIMAINN)
* initialize IIMAINN to be -1 in DT_PICKSRC if not in SRC. 
* otherwise, it would be the index for the SRC partner
      
      !using IFMPOST to switch on applying LF kinematics to E and pz
      if(IFMPOST .EQ. 2) then 
        if(NINT(INUMOD).NE.2) then
          STOP "FATAL: CAN ONLY DO DEUTERON KINEMATICS"    
        else
          CALL DT_SPECTRALFUNC(PHKK,NINT(INUMOD),IIMAIN)
        endif
      endif

      IF( IIMAINN.NE.-1 ) THEN
        ISTHKK(IIMAINN)=-12
        NINTS=NINTS+1
        IINTER(NINTS)=IIMAINN
      ENDIF

C... Note: DPF(mu) = P(mu)_true - P(mu)_naive is a 4-momentum too.   
C    DPF is the name in the HCMS
C    PXF,PYF,PZF,EKF in the TRF
      DPF(1) = PXF
      DPF(2) = PYF
      CALL DT_LTNUC(PZF,EKF,DPF(3),DPF(4),3)

C  ROBO Pythia event into same frame as PHKK (TRF g*=z,e' px>0,py-0) 
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         write(*,*) "DT_PYEVNTEP: Lab, z=p"
         CALL PYLIST(2)
         write(*,*) ""
         write(*,*) "P(mu)_true - P(mu)_naive (z along gamma*):"
         write(*,*) "TRF(1-4): ",PXF,PYF,PZF,EKF
         write(*,*) "HCMS(1-4): ",DPF(1),DPF(2),DPF(3),DPF(4)
         write(*,*) ""
      endif
      CALL PYROBO(0,0,0.0D0,0.0D0,0.0D0,0.0D0,-pbeta)
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         write(*,*) "DT_PYEVNTEP: TRF e=-z"
         CALL PYLIST(2)
      endif
C...Rotate so that the gamma* (or Z0) is along +z            
      DO ITRK=1,N
         IF ((K(ITRK,2).EQ.22 .OR. K(ITRK,2).EQ.23) 
     &        .AND. K(ITRK,1).EQ.21 .AND. K(ITRK,3).EQ.1) THEN
            PHIGAM=PYANGL(P(ITRK,1),P(ITRK,2))
            PTGAM=DSQRT(P(ITRK,1)**2+P(ITRK,2)**2)
            THEGAM=PYANGL(P(ITRK,3),PTGAM)
            GOTO 30
         ENDIF
      ENDDO
 30   CONTINUE
C            write(*,*) "TRF e=-z, g* theta, phi, COD, SID, COF, SIF: ",
C     &        THEGAM," ",PHIGAM," ",DCOS(THEGAM)," ",DSIN(THEGAM),
C     &                              DCOS(PHIGAM)," ",DSIN(PHIGAM)
C... Rotate to TRF, gamma* along z, e' px>0, py=0
      CALL PYROBO(0,0,0.0D0,PARU(1)-PHIGAM,0.D0,0.D0,0.D0)
      CALL PYROBO(0,0,THEGAM,0.0D0,0.D0,0.D0,0.D0)
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         write(*,*) "DT_PYEVNTEP: TRF g*=z, e' px>0 py=0"
         CALL PYLIST(2)
      endif
C
C If needed:
C Calculate kinematics for Fermi momentum correction.
C We are in the nuclear TRF with gamma* along z. 
c VINT(1-4) is W, W2, -SQRT(Q2), M_N

C     Find and flag the spectator nucleon for the Deuteron
      IF (ITZ.EQ.1 .AND. IT.EQ.2) THEN
         IF(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
            WRITE(*,*)"I  ISTHKK(I)  IDHKK(I) PHKK(5,I)"
            DO IHKK=1,NHKK
               WRITE(*,*)IHKK,ISTHKK(IHKK),IDHKK(IHKK),PHKK(5,IHKK)
            ENDDO
         ENDIF
         DO IHKK=1,NHKK 
            IF (ISTHKK(IHKK).EQ.14) THEN
               IIMAINN = IHKK
               GOTO 31
            ENDIF
         ENDDO
         STOP 'FATAL: Could not find a spectator nucleon in deuteron.'
         IF(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
            WRITE(*,*)"I=IIMAINN  ISTHKK(I)  IDHKK(I) PHKK(5,I)"
            WRITE(*,*)IIMAINN,ISTHKK(IIMAINN),IDHKK(IIMAINN),
     &           PHKK(5,IIMAINN)
         ENDIF
 31      CONTINUE
      ENDIF

      if (USERSET.EQ.3 .OR. IFERPY.GT.1) then
         VALNU = 0.5D0*VINT(309)*(VINT(302)-VINT(4)**2-VINT(303)**2)
     &        /VINT(4)
C     For now, the deuteron is at rest in the IRF (may change for SRC)
C     Not used unless IFERPY>2
         PDEUT(1)=0.0D0
         PDEUT(2)=0.0D0
         PDEUT(3)=0.0D0
         PDEUT(5)=AZMASS(IT,ITZ,ITMMOD)
         PDEUT(4)=PDEUT(5)
         MNUCL=VINT(4)
         DO IDIM=1,5
            PSPEC(IDIM)=0.0D0
         ENDDO
         if (IFERPY.EQ.2) then
            CALL PFSHIFT(VINT(1),VINT(2),VINT(307),VALNU,MNUCL,PDEUT,
     &           PSPEC,1,IREJ)
         elseif (IFERPY.EQ.3) then
            if (IIMAINN.NE.-1) then
C     Here MNUCL is the SPECTATOR (not struck) nucleon mass from D or SRC-pair:
               MNUCL=PHKK(5,IIMAINN)
               DO IDIM=1,5
                  PSPEC(IDIM)=PHKK(IDIM,IIMAINN)
               ENDDO
               PDEUT(1)=PHKK(1,IIMAIN)+PHKK(1,IIMAINN)
               PDEUT(2)=PHKK(2,IIMAIN)+PHKK(2,IIMAINN)
               PDEUT(3)=PHKK(3,IIMAIN)+PHKK(3,IIMAINN)
               PDEUT(4)=SQRT(PDEUT(1)*PDEUT(1)+PDEUT(2)*PDEUT(2)+
     &              PDEUT(3)*PDEUT(3)+PDEUT(5)*PDEUT(5))
               CALL PFSHIFT(VINT(1),VINT(2),VINT(307),VALNU,MNUCL,PDEUT,
     &              PSPEC,2,IREJ)
            else
               CALL PFSHIFT(VINT(1),VINT(2),VINT(307),VALNU,MNUCL,PDEUT,
     &              PSEPC,1,IREJ)
            endif
         else
            CALL PFSHIFT(VINT(1),VINT(2),VINT(307),VALNU,MNUCL,PDEUT,
     &           PSPEC,0,IREJ)
         endif
         if (IREJ.NE.0) goto 9999
         if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
            write(*,*) "DT_PYEVNTEP: TRF g*=z, after pf post-fix"
            CALL PYLIST(2)
         endif
      endif

C... If requested, fix the e+D event kinematics.
      ! change to .EQ. 1 instead of .GE. 1
      if (IFMPOST.EQ.1) then
         if (ITZ.NE.1 .OR. IT.NE.2) 
     &        STOP "FATAL: CAN ONLY POST-FIX DEUTERON KINEMATICS"
         if (IFERPY.GT.2) 
     &        STOP "FATAL: IFERPY>2 and IFMPOST>0 are incompatible"
C...     Copy the spectator nucleon to the PYTHIA event record 
C...     and flag it as having been involved in the interaction 
C...     but not needing a ptkick using temporary value -12
C     nu = y(s-M2-m2)/2M  where M=M_nucleon and m=m_lepton
         NINTS=NINTS+1
         IINTER(NINTS)=IIMAINN
         ISTHKK(IIMAINN)=-12
         N = N + 1
         P(N,1) = PHKK(1,IINTER(NINTS))
         P(N,2) = PHKK(2,IINTER(NINTS))
         P(N,3) = PHKK(3,IINTER(NINTS))
         P(N,4) = PHKK(4,IINTER(NINTS))
         P(N,5) = PHKK(5,IINTER(NINTS))
         K(N,1) = 1
         K(N,2) = IDHKK(IINTER(NINTS))
         K(N,3) = IINTER(NINTS)
         K(N,4) = 0
         K(N,5) = 0
C...  Note: Can't fill V(N,I) yet or PYROBO will boost it around.
         NPOS(NINTS) = N
         PosAlt(NINTS,1) = VHKK(1,IINTER(NINTS))
         PosAlt(NINTS,2) = VHKK(2,IINTER(NINTS))
         PosAlt(NINTS,3) = VHKK(3,IINTER(NINTS))
         PosAlt(NINTS,4) = VHKK(4,IINTER(NINTS))
         MomAlt(NINTS) = IINTER(NINTS)
         VALNU = 0.5D0*VINT(309)*(VINT(302)-VINT(4)**2-VINT(303)**2)
     &        /VINT(4)
         QQ=VINT(307)
         CALL DEUTFIX(VALNU,QQ,AZMASS(2,1,1))
      endif

C
C...Struck "parton" is the (non e') particle with highest pz (along g*)
      PZMAX=0.D0
      IPARTN=0
      DO ITRK=1,N
         IF ( (1.LE.K(ITRK,1) .AND. K(ITRK,1).LE.3) .AND.
     &       .NOT.(11.LE.ABS(K(ITRK,2)).AND.ABS(K(ITRK,2)).LE.18) ) THEN
            IF (P(ITRK,3).GT.PZMAX) THEN
               PZMAX=P(ITRK,3)
               IPARTN=ITRK
            ENDIF
         ENDIF
      ENDDO
c...  For the non-main interactions, just do a ptkick & recoil
C...  Copy the recoiling nucleon to the PYTHIA event record 
C...  and flag it as special.
C...  If it is a spectator nucleon from e+D, leave it alone.
C...  If it is a spectator nucleon from e+A, put in the event record, but 
C...  don't ptkick it...
      DO IINT=1,NINTS
         IF (ISTHKK(IINTER(IINT)).EQ.-12) THEN
            ISTHKK(IINTER(IINT))=12
            IF (IT.GT.2) THEN
               N = N + 1
               P(N,1) = PHKK(1,IINTER(IINT))
               P(N,2) = PHKK(2,IINTER(IINT))
               P(N,3) = PHKK(3,IINTER(IINT))
               P(N,4) = PHKK(4,IINTER(IINT))
               P(N,5) = PHKK(5,IINTER(IINT))
               K(N,1) = 1
               K(N,2) = IDHKK(IINTER(IINT)) 
               K(N,3) = IINTER(IINT)
               K(N,4) = 0
               K(N,5) = 0
C...  Note: Can't fill V(N,I) yet or PYROBO will boost it around.
               NPOS(IINT) = N
               PosAlt(IINT,1) = VHKK(1,IINTER(IINT))
               PosAlt(IINT,2) = VHKK(2,IINTER(IINT))
               PosAlt(IINT,3) = VHKK(3,IINTER(IINT))
               PosAlt(IINT,4) = VHKK(4,IINTER(IINT))
               MomAlt(IINT) = IINTER(IINT)
            ENDIF
         ELSEIF (IINT.NE.IMAIN) THEN
            N = N + 1
            P(N,1) = PHKK(1,IINTER(IINT))
            P(N,2) = PHKK(2,IINTER(IINT))
            P(N,3) = PHKK(3,IINTER(IINT))
            P(N,4) = PHKK(4,IINTER(IINT))
            P(N,5) = PHKK(5,IINTER(IINT))
            K(N,1) = 1
            K(N,2) = IDHKK(IINTER(IINT)) 
            K(N,3) = IINTER(IINT)
            K(N,4) = 0
            K(N,5) = 0
C...  Note: Can't fill V(N,I) yet or PYROBO will boost it around.
            NPOS(IINT) = N
            PosAlt(IINT,1) = VHKK(1,IINTER(IINT))
            PosAlt(IINT,2) = VHKK(2,IINTER(IINT))
            PosAlt(IINT,3) = VHKK(3,IINTER(IINT))
            PosAlt(IINT,4) = VHKK(4,IINTER(IINT))
            MomAlt(IINT) = IINTER(IINT)
            CALL KICKIT(P(IPARTN,1), P(IPARTN,2), P(IPARTN,3), 
     &           P(IPARTN,4), P(IPARTN,5),
     &           P(N,1), P(N,2), P(N,3), P(N,4), P(N,5))
         ENDIF
      ENDDO            
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         write(*,*) "DT_PYEVNTEP: TRF g*=z, post-ptkick pre-ApplyQW"
         CALL PYLIST(2)
         WRITE(*,*) "PYQREC=0: ",PYQREC(1),PYQREC(2),PYQREC(3),
     &        PYQREC(4)
      endif
C Decay particles (like a J/psi) from the event skeleton
      MSTJ(21)=2
      CALL PYEXEC
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*)"PYLIST: After PYEXEC"
         CALL PYLIST(2)
      endif
      if(IOULEV(6).GE.1) print*, 'Event ', NEVENT
c...  do quenching to scattered partons if requested      
      IF (QSWITCH.EQ.1) THEN
         call InterPos          ! Set the interaction position in nucleus
         if(IOULEV(6).GE.1) then
            print*, 'before pyqm'
            call PYLIST(1)
         endif
         call ApplyQW(QHAT)     ! Compute QW (also fill PYQREQ(mu))
         if(IOULEV(6).GE.1) then
            print*,'after pyqm'
            call PYLIST(1)
         endif 
C     2017-08-26 MDB For Userset1, use PYQREC in TRF z along gamma*
         IF (USERSET.EQ.1 .OR. USERSET.EQ.2) THEN
            USER3 = PYQREC(4)
            IF (USERSET.EQ.1) THEN
               USER1 = PYQREC(1)*PYQREC(1)+PYQREC(2)*PYQREC(2)
               USER2 = DSQRT(USER1 +PYQREC(3)*PYQREC(3))
               USER1 = DSQRT(USER1)
            ENDIF
         ENDIF
      ENDIF
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         write(*,*) 'DT_PYEVNTEP: TRF g*=z, After ApplyQW'
         CALL PYLIST(2)
      endif
      CALL PYROBO(0,0,-THEGAM,PHIGAM-PARU(1),0.0D0,0.0D0,0.0D0)
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         write(*,*) 'DT_PYEVNTEP: Back to TRF e=-z'
         CALL PYLIST(2)
      endif
      CALL PYROBO(0,0,0.0D0,0.0D0,0.0D0,0.0D0,pbeta)
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         write(*,*) 'DT_PYEVNTEP: Back in lab frame'
         CALL PYLIST(2)
      endif
C We are done with the partons. Hadronize the event now.
      MSTJ(1) = 1
      CALL PYEXEC
      if(NEVENT.LE.IOULEV(5)) then
         if(IOULEV(4).GE.2) then
            write(*,*) 'DT_PYEVNTEP: After  PYEXEC'
            CALL PYLIST(2)
         endif
         if(IOULEV(4).GE.2) then
            WRITE(*,*) 'Recoil PYQREC(1-4) TRF: ',PYQREC(1),PYQREC(2),
     &           PYQREC(3),PYQREC(4)
         endif
      endif
C MDB 2017-08-07 Rename PF as PYQREC.
C      Boost directly from TRF z=g* to HCMS z=g* in one shot...
C      WRITE(*,*) 'Boost from TRF to HCMS. gamma, betagamma:',
C     &      GACMS(2),BGCMS(2)
      P3=PYQREC(3)
      P4=PYQREC(4)
      PYQREC(3)=GACMS(2)*P3-BGCMS(2)*P4
      PYQREC(4)=GACMS(2)*P4-BGCMS(2)*P3
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*) 'Mycalc Recoil PYQREC(1-4) HCMS: ',PYQREC(1),
     &        PYQREC(2),PYQREC(3),PYQREC(4)
      endif
      CALL DT_LTNUC(P3,P4,PYQREC(3),PYQREC(4),3)
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*) 'DT_LTNUC Recoil PYQREC(1-4) HCMS: ',PYQREC(1),
     &        PYQREC(2),PYQREC(3),PYQREC(4)
      endif
c...Translate PYJETS into HEPEVT event record
      CALL PYHEPC(1)

c...Output part
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*) 'Listing HEPEVT as we go. Before rotation which',
     &        'takes x->-x and z->-z.'
         WRITE(*,*) 'J, ISTHEP(J), IDHEP(J), JMOHEP(1-2,J), ',
     &        'JDAHEP(1-2,J)','PHEP(1-5,J)'
      endif
c...First loop to find exchanged boson
      DO  J=1,NHEP
         if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
c...output of event list from pythia for a temp check      
            WRITE(*,998) J,ISTHEP(J),IDHEP(J),JMOHEP(1,J),
     &           JMOHEP(2,J),JDAHEP(1,J),JDAHEP(2,J),
     &           PHEP(1,J),PHEP(2,J),PHEP(3,J),PHEP(4,J),
     &           PHEP(5,J)
 998        FORMAT(I5,I5,I8,4I5,5F17.5)
         endif
c...find the virtual photon to do LT from lab to gamma c.m.s      
         PHEP(3,J)=-PHEP(3,J)
c... Mark 2016-09-14 Flip px also to make it a rotation
         PHEP(1,J)=-PHEP(1,J)  
         IF((IDHEP(J).EQ.22 .OR. IDHEP(J).EQ.23) .AND.
     &        (ISTHEP(J).EQ.3).AND. (JMOHEP(1,J).EQ.1)) THEN
            GAMM(1)=PHEP(1,J)
            GAMM(2)=PHEP(2,J)
            GAMM(3)=PHEP(3,J)
            GAMM(4)=PHEP(4,J)
            GAMM(5)=PHEP(5,J)
         ENDIF
      ENDDO
c...loop to find exchanged boson end

**************copy entries from pythia to DTEVT1*************
*********transform pythia entries from lab to c.m.s of gamma+p********
      eveBETA(1)=-(PHEP(1,2)+GAMM(1))/(PHEP(4,2)+GAMM(4))
      eveBETA(2)=-(PHEP(2,2)+GAMM(2))/(PHEP(4,2)+GAMM(4))

**!!!!!!!!remember to change the sign of proton beam in LT*****      
      eveBETA(3)=-(PHEP(3,2)+GAMM(3))/(PHEP(4,2)+GAMM(4))

      eveBETA(4)=eveBETA(1)*eveBETA(1)+eveBETA(2)*eveBETA(2)+
     & eveBETA(3)*eveBETA(3)

      GAA=1./SQRT(1-eveBETA(4))
      
      eveBETA(1)=eveBETA(1)*GAA
      eveBETA(2)=eveBETA(2)*GAA
      eveBETA(3)=eveBETA(3)*GAA
*  in the n rest frame rotate virtual photon  angles to +z axis
*...COD=cos(theta) SID=sin(theta) COF=cos(phi) SIF=sin(phi)
      call DT_DALTRA(GAA,eveBETA(1),eveBETA(2),eveBETA(3),
     &GAMM(1),GAMM(2),GAMM(3),GAMM(4),PTOT,P1,P2,P3,P4)
      PTOT=SQRT(P1**2+P2**2+P3**2)
      COD = P3/PTOT
      PPT = SQRT(P1**2+P2**2)
      SID = PPT/PTOT
      IF(P1.GT.ZERO) THEN
         COF = ONE
      ELSE
         COF = -ONE
      ENDIF
      SIF = ZERO      
      IF (PTOT*SID.GT.TINY10) THEN
         COF = P1/(SID*PTOT)
         SIF = P2/(SID*PTOT)
         ANORF = SQRT(COF*COF+SIF*SIF)
         COF = COF/ANORF
         SIF = SIF/ANORF
      ENDIF
C
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*)'Transformation from Lab to HCMS'
         WRITE(*,*)'GAA,eveBETA(1-3): ',GAA,' ',eveBETA(1),' ',
     &        eveBETA(2),' ',eveBETA(3)
         WRITE(*,*)'Lab GAMM(1-4): ',GAMM(1),' ',GAMM(2),' ',GAMM(3),
     &        GAMM(4)
         WRITE(*,*)'HCMS P1,P2,P3,P4: ',P1,' ',P2,' ',P3,' ',P4
         WRITE(*,*)'Photon angles in HCMS'
         WRITE(*,*)'COD,SID,COF,SIF: ',COD,' ',SID,' ',COF,' ',SIF
      endif
c...Collector to get the sum of final state momentums other than 
c...the scattered electron in gamma nucleon c.m.s frame
      PP1=0.
      PP2=0.
      PP3=0.
      PP4=0.
C...Mark also charge 2016-09-10
      QCHG=0.
***************from HEPEVT to HKKEVT***************************
C...Mark - 2017-02-22 Change defn. of NPOINT(1)
      DO J=1,NHEP
         I=J+NPOINT(1)
         ISTHKK(I)=ISTHEP(J)
         IDHKK(I)=IDHEP(J)
         IF(JMOHEP(1,J).GE.1) JMOHKK(1,I)=JMOHEP(1,J)+NPOINT(1)
         IF(JMOHEP(2,J).GE.1) JMOHKK(2,I)=JMOHEP(2,J)+NPOINT(1)
         IF(JDAHEP(1,J).GE.1) JDAHKK(1,I)=JDAHEP(1,J)+NPOINT(1)
         IF(JDAHEP(2,J).GE.1) JDAHKK(2,I)=JDAHEP(2,J)+NPOINT(1)
**********rotate in nucleon c.m.s frame***********************         
         call DT_DALTRA(GAA,eveBETA(1),eveBETA(2),eveBETA(3),
     &   PHEP(1,J),PHEP(2,J),PHEP(3,J),PHEP(4,J),PTOT,P1,P2,P3,P4)
         PHKK(1,I)=COD*(COF*P1+SIF*P2)-SID*P3
C  Mark 2016-09-15 Make this a rotation, not a flip
C        theta=phi=0 should be the identity transformation.
C        Old line: PHKK(2,I)=SIF*P1-COF*P2
         PHKK(2,I)=COF*P2-SIF*P1
         PHKK(3,I)=SID*(COF*P1+SIF*P2)+COD*P3
         PHKK(4,I)=P4
c         call DT_LTNUC(P3,P4,PHKK(3,I),PHKK(4,I),3)
*********LT from nucleon rest c.m.s**************************
c         call DT_LTRANS(P1,P2,P3,P4,PHKK(1,I),PHKK(2,I),
c     &    PHKK(3,I),PHKK(4,I),1,3)
******rotate to virtual photon directing z+ ******************         
         PHKK(5,I)=PHEP(5,J)
c         PHKK(1,I)=P1
c         PHKK(2,I)=P2
C
      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*)'Particle I, ISTHKK, IDHKK: ',I,' ',ISTHKK(I),
     &             ' ', IDHKK(I)
         WRITE(*,*)'Lab frame'
         WRITE(*,*),'PHEP(1-4): ',PHEP(1,J),' ',PHEP(2,J),' ',
     &                            PHEP(3,J),' ',PHEP(4,J)
         WRITE(*,*)'HCMS: boosted from lab'
         WRITE(*,*)'P1,P2,P3,P4: ',P1,' ',P2,' ',P3,' ',P4
         WRITE(*,*)'HCMS: z along gamma*'
         WRITE(*,*),'PHKK(1-4): ',PHKK(1,I),' ',PHKK(2,I),' ',
     &                            PHKK(3,I),' ',PHKK(4,I)
      endif
********get the position of particles in nucleon rest frame***         
c... we simply use the position of involved nucleon for all the
c... particles         
c... Mark - 2017-02-28 treat recoiling extra nucleons differently
         JPOINT = 0
         DO MM=1,NINTS
            IF (NPOS(MM).EQ.J) THEN
               JPOINT = MM
            ENDIF
         ENDDO
         IF (JPOINT.EQ.0) THEN         
            DO M=1,4
               VHKK(M,I)=PosNuc(M)
            ENDDO
         ELSE
c... Set VHKK for recoiling nucleons.
            DO M=1,4
               VHKK(M,I)=PosAlt(JPOINT,M)
            ENDDO
            JMOHKK(1,I)=MomAlt(JPOINT)
         ENDIF

c...set BAM ID for the particles         
         IDBAM(I)=IDT_ICIHAD(IDHKK(I))
c...change the IS of scattered lepton from 1 to 99 in order to avoid its 
c...interaction in cascade. Mother is 3 for Pythia, 1 for RAPGAP
         IF( (ISTHEP(J).EQ.1).AND.
     &        (JMOHEP(1,J).EQ.3 .OR. JMOHEP(1,J).EQ.1).AND.
     &        (ABS(IDHEP(J)).EQ.11.OR.ABS(IDHEP(J)).EQ.13) ) THEN
            ISTHKK(I)=99
C            JMOHKK(1,I)=JMOHEP(1,J)
         ENDIF
         NHKK=NHKK+1

c...collect final state momentum in gamma*p cms frame
         IF( ISTHKK(I).EQ.1 ) THEN
            PP1=PHKK(1,I)+PP1
            PP2=PHKK(2,I)+PP2
            PP3=PHKK(3,I)+PP3
            PP4=PHKK(4,I)+PP4
            QCHG=QCHG+PYCHGE(IDHKK(I))/3D0
C            print*,'I=',I,' ID=',IDHKK(I),' IS=',ISTHKK(I)
C            print*,'BAM=',IDBAM(I),' P1=',PHKK(1,I),' P2=',PHKK(2,I)
C            print*,'NOBAM=',NOBAM(I),' P3=',PHKK(3,I),' P4=',PHKK(4,I)
         ENDIF
C      WRITE(89,997) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),
C     &          JMOHKK(2,I),JDAHKK(1,I),JDAHKK(2,I),
C     &          PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
C     &              PHKK(5,I)
C  997      FORMAT(I5,I5,I8,4I5,5F17.5)
      ENDDO

C
      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
         print*,'4-momentum & charge totals for Status=1:'
         print*,'PP1=',PP1,' PP2=',PP2,' PP3=',PP3,' PP4=',PP4,
     &        ' QCHG=',QCHG
      endif
c      DO M=1,4
c         print*,'P1(',M,')=',PHKK(M,NPOINT(1)+5),'
c     &    P2(',M,')=',PHKK(M,NPOINT(1)+6)
c      ENDDO
c      print*,'NPOINT1=',NPOINT(1)
c      print*,'ID1=',IDHKK(NPOINT(1)+5),' ID2=',IDHKK(NPOINT(1)+6)

C MDB 2017-02-28 This method can't be used for genShd>1
      IF (genShd.EQ.1) THEN
      !!!PAUX is used to balance the PFSP used to estimate 
      !!!total final state particle information
      !!!if energy loss, the lost momenta must be saved for 
      !!!to correctly calculate the nucleus remnant mass
      !!!in the DPMJET routines, where only evaporation energy
      !!!should be considered
         PAUX(1)=PHKK(1,NPOINT(1)+4)+PHKK(1,NPOINT(1)+5)-PP1
         PAUX(2)=PHKK(2,NPOINT(1)+4)+PHKK(2,NPOINT(1)+5)-PP2
         PAUX(3)=PHKK(3,NPOINT(1)+4)+PHKK(3,NPOINT(1)+5)-PP3
         PAUX(4)=PHKK(4,NPOINT(1)+4)+PHKK(4,NPOINT(1)+5)-PP4
      ELSE
         PAUX(1)=0.0
         PAUX(2)=0.0
         PAUX(3)=0.0
         PAUX(4)=0.0
      ENDIF

      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*) 'End of DT_PYEVNTEP - mode 2'
         CALL DT_PYOUTEP(4)
      endif

      LFIRST=.TRUE.
9999  CONTINUE
      IF (IOULEV(1).GE.1) WRITE(*,*) 'DT_PYEVNTEP IREJ= ', IREJ
      RETURN

      END


*=====dt_pyout=========================================================
*used for the output of pythia event list and statistics information
      SUBROUTINE DT_PYOUTEP(MODE)     
 
*     input:
*           MODE: 1:reject statistics - not really used
*                 2:event output
*                 3:total statistics print
*                 4:event output to screen (for debugging) Mark 08/17/2016

      include 'pythia.inc'              ! All PYTHIA commons blocks
      include "mc_set.inc"
      include "py6strf.inc"
      include "mcRadCor.inc"
      include "radgen.inc"
      include "phiout.inc"
      include "beagle.inc"

      EXTERNAL PYCHGE, NBARY, PYMASS
      INTEGER  NBARY

* event history
      PARAMETER (NMXHKK=200000)
   
      LOGICAL ISHADR
      EXTERNAL ISHADR

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history - IHIST renamed because it's in pythia.inc
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHISTDPM(2,NMXHKK)

      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY10=1.0D-10,MAXNCL=260)

* Glauber formalism: collision properties
      COMMON /DTGLCP/ RPROJ,RTARG,BIMPAC,
     &                NWTSAM,NWASAM,NWBSAM,NWTACC,NWAACC,NWBACC

C...Pythia event counter (since we keep PYINITing) Mark 2017-01-31
      COMMON /PYCNTR/ MYNGEN
      INTEGER MYNGEN
      SAVE /PYCNTR/

* lorentz transformation parameter
      COMMON /LTPARA/ BGTA(4), GAMM(5), eveBETA(4), GAA
      COMMON /ROTATE/ COF,COD,SIF,SID

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 

* added by liang to store the output event variables 1/20/12
      COMMON /EVTOUT/ XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT
      DOUBLE PRECISION XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT

C...Added by liang 1/6/12
C...Switches for nuclear correction
      COMMON /PYNUCL/ INUMOD,CHANUM,ORDER,genShd
      SAVE /PYNUCL/
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER,genShd

* treatment of residual nuclei: 4-momenta
      LOGICAL LRCLPR,LRCLTA
      COMMON /DTRNU1/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5),
     &                TRCLPR(5),TRCLTA(5),LRCLPR,LRCLTA

* treatment of residual nuclei: properties of residual nuclei
      COMMON /DTRNU2/ AMRCL0(2),EEXC(2),EEXCFI(2),
     &                NTOT(2),NPRO(2),NN(2),NH(2),NHPOS(2),NQ(2),
     &                NTOTFI(2),NPROFI(2)

* treatment of residual nuclei: wounded nucleons
      COMMON /DTWOUN/ NPW,NPW0,NPCW,NTW,NTW0,NTCW,IPW(210),ITW(210)

* properties of interacting particles
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5) 
c...target/proj mass, charge and projectile internal ID
      integer IT, ITZ, IP, IPZ, IJPROJ, IBPROJ, IJTARG, IBTARG, ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP 

* added by liang to check the photon flux 12/28/11
      COMMON /FLCHK/ PFXCHK
      DOUBLE PRECISION PFXCHK

C...output file name definition
      COMMON /OUNAME/ outname

      CHARACTER*256 outname

      DOUBLE PRECISION P1, P2, P3, P4, PTOT, PP1, PP2, PP3, PP4

      integer NEV, NPRT, ievent, genevent, I, tracknr
      integer lastgenevent, idum1, idum2, initseed, nrtrack
      DOUBLE PRECISION radgamE, radgamp, radgamEnucl
      DOUBLE PRECISION epznucl
      CHARACTER PARAM*100
      LOGICAL UseLut, GenLut
      INTEGER MODE

C     Local
      INTEGER IDOUT
      DOUBLE PRECISION BETA2,P2TEMP,XCALC
      DOUBLE PRECISION P5SUM(5)
      LOGICAL VERBOSE

      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./

      if (mcRadCor_EBrems.gt.0.) then
         radgamEnucl=sqrt(dplabg(1)**2+dplabg(2)**2+dplabg(3)**2)
         radgamE=pgamma*radgamEnucl-pgamma*pbeta*dplabg(3)
         radgamp=-pgamma*pbeta*radgamEnucl+pgamma*dplabg(3)
C         write(*,*) radgamEnucl, radgamE, dplabg(3), radgamp
      else
        radgamEnucl=0D0
        radgamE=0D0
        radgamp=0D0 
      endif

      if ((msti(1).ge.91).and.(msti(1).le.94)) msti(16)=0

      ievent = NEVENT

      GOTO (1,2,3,4) MODE

c...mode 1 is used to update the reject statistics in pythia
1     CONTINUE      
c      write(99,*),'event ',ievent,' rejected,',' proces=',
c     & msti(1),', X=',XBJOUT,' Q2=',Q2OUT 
         
      RETURN

c...mode 2 is used to output the event list
2     CONTINUE

      IF(FIRST) open(29, file=outname,STATUS='UNKNOWN')

      AREMN = 0
      NNEVAP = 0
      NPEVAP = 0
      NINC = 0
      NINCCH = 0
**************check HKKEVT***********************************
      DO  J=1,NHKK
C Flag INC particles here. Count them below.
         IF (IDCH(J).GT.0) NOBAM(J)=10+IDCH(J)
********rotate back from the gamma* +z direction***********
         IF(J.GT.NPOINT(1)) THEN
            P1=COF*(COD*PHKK(1,J)+SID*PHKK(3,J))+SIF*PHKK(2,J)
            P2=SIF*(COD*PHKK(1,J)+SID*PHKK(3,J))-COF*PHKK(2,J)
            P3=COD*PHKK(3,J)-SID*PHKK(1,J)
            P4=PHKK(4,J)
c**************transform back to the lab frame**************      
            call DT_DALTRA(GAA,-eveBETA(1),-eveBETA(2),-eveBETA(3),
     &           P1,P2,P3,P4,PTOT,PP1,PP2,PP3,PP4)
c     WRITE(89,996) J,ISTHKK(J),IDHKK(J),JMOHKK(1,J),
c     &          JMOHKK(2,J),JDAHKK(1,J),JDAHKK(2,J),
c     &          PP1,PP2,-PP3,PP4, !remember to change the sign of z back
c     &              PHKK(5,J)
c  996      FORMAT(I5,I5,I8,4I5,5F17.5)
            PHKK(1,J)=PP1
            PHKK(2,J)=PP2      
            PHKK(3,J)=-PP3
            PHKK(4,J)=PP4
c     PHKK(3,J)=pgamma*(P3+pbeta*P4)
c     PHKK(4,J)=pgamma*(P4+pbeta*P3)
c...find the exchanged boson and out e- to make it fit root tree making rules
c...in the following steps
            IF((ISTHKK(J).EQ.3).AND.(IDHKK(J).EQ.22 .OR. IDHKK(J).EQ.23)
     &           .AND. (JMOHKK(1,J).EQ.(NPOINT(1)+1))) THEN
               IBOSON=J
            ELSEIF(ISTHKK(J).EQ.99) THEN
               ISTHKK(J)=1
C               JMOHKK(1,J)=3   ! Do this later. Don't change internal record.
               ILEPT=J
c...2017-01-02 MDB Fill some new event variables
            ELSEIF (.NOT. OLDOUT) THEN
               IF (ISTHKK(J).EQ.1001) THEN
                  AREMN = IDRES(J)
                  IF (USERSET.EQ.5) USER1=IDXRES(J)
               ELSEIF (ISTHKK(J).EQ.-1) THEN
                  IF (IDHKK(J).EQ.2212) THEN
                     NPEVAP=NPEVAP+1
                  ELSEIF (IDHKK(J).EQ.2112) THEN
                     NNEVAP=NNEVAP+1
                  ENDIF
               ELSEIF (ISTHKK(J).EQ.1.AND.NOBAM(J).GT.10) THEN
                  NINC = NINC + 1
                  IF (IDHKK(J).EQ.80000) THEN
                     IF (IDXRES(J).NE.0) NINCCH = NINCCH + 1
                  ELSE
                     IF (PYCHGE(IDHKK(J)).NE.0) NINCCH = NINCCH + 1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF ! (J.GT.(NPOINT(1)+1))
      ENDDO ! J=1,NHKK
*************check HKKEVT end***************************      
      !liang-2021-10 perform long-lifetime decay
      CALL DT_DECALL

      genevent=MYNGEN-lastgenevent
      tracknr = NHKK
      if (mcRadCor_EBrems.gt.0.) then
         nrtrack=tracknr+1
      else
         nrtrack=tracknr
      endif

c...print a title for the event file
      If (FIRST) then
C        write(29,*)' PYTHIA EVENT FILE '
        write(29,*)' BEAGLE EVENT FILE '
        write(29,*)'============================================'
        if (OLDOUT) then
           write(29,30)
        else
           write(29,31) 
        endif
C OLD OLD
C 30     format('I, ievent, genevent, subprocess, nucleon,
C     &  targetparton, xtargparton, beamparton, xbeamparton,
C     &  thetabeamprtn, truey, trueQ2, truex, trueW2, trueNu, leptonphi, 
C     &  s_hat, t_hat, u_hat, pt2_hat, Q2_hat, F2, F1, R, sigma_rad, 
C     &  SigRadCor, EBrems, photonflux, nrTracks')
 30     format('I, ievent, genevent, lepton, Atarg, Ztarg, pzlep,pztarg,
     &  pznucl, subprocess, 
     &  nucleon, targetparton, xtargparton, beamparton, xbeamparton,
     &  thetabeamprtn, truey, trueQ2, truex, trueW2, trueNu, leptonphi, 
     &  s_hat, t_hat, u_hat, pt2_hat, Q2_hat, F2, F1, R, sigma_rad, 
     &  SigRadCor, EBrems, photonflux, b, Phib, Thickness, ThickScl, 
     &  Ncollt, Ncolli, Nnevap, Npevap, Aremn, NINC, NINCch, d1st, davg,
     &  pxf, pyf, pzf, Eexc, RAevt, User1, User2, User3, nrTracks')
 31     format('I, ievent, genevent, lepton, Atarg, Ztarg, pzlep,pztarg,
     &  pznucl, crang, crori, subprocess, nucleon, targetparton,
     &  xtargparton, beamparton, xbeamparton, thetabeamprtn, truey, 
     &  trueQ2, truex, trueW2, trueNu, leptonphi, s_hat, t_hat, u_hat,
     &  pt2_hat, Q2_hat, F2, F1, R, sigma_rad, SigRadCor, EBrems, 
     &  photonflux, b, Phib, Thickness, ThickScl, Ncollt, Ncolli,
     &  Nwound, Nwdch, Nnevap, Npevap, Aremn, NINC, NINCch, d1st, davg,
     &  pxf, pyf, pzf, Eexc, RAevt, User1, User2, User3, nrTracks')
        write(29,*)'============================================'

c...similar to the dpmjet track wide title 
      write(29,*)'I  ISTHKK(I)  IDHKK(I)  JMOHKK(2,I)  JMOHKK(1,I)
     & JDAHKK(1,I)  JDAHKK(2,I)  PHKK(1,I)  PHKK(2,I)  PHKK(3,I)
     & PHKK(4,I)  PHKK(5,I)  VHKK(1,I) VHKK(2,I) VHKK(3,I) IDRES(I)
     & IDXRES(I) NOBAM(I)'

c        write(29,*)' I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)
c     &  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)'
        write(29,*)'============================================'
         FIRST=.FALSE.
      endif

      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
         XCALC = VINT(307)/VINT(309)/(VINT(302)-VINT(4)**2-VINT(303)**2)
         write(*,*) "MN VINT(4): ",VINT(4),"massp: ",massp
         write(*,*) "MN VINT(304): ",VINT(304),"MN2 VINT(308)",VINT(308)
         write(*,*) "X(VINT-calc): ",XCALC," XBJOUT: ",XBJOUT
         write(*,*) "beamx VINT(305): ",VINT(305),
     &        "targx VINT(306): ",VINT(306)
         write(*,*) "Y VINT(309): ",VINT(309)," YYOUT:  ",YYOUT
         write(*,*) "W2 VINT(2): ",VINT(2)," W2OUT: ",W2OUT
         write(*,*) "Nu(from X): ",VINT(307)/(2.*VINT(4)*XCALC),
     &        " NUOUT: ",NUOUT
         write(*,*) "Q2 VINT(307): ",VINT(307)," Q2OUT: ",Q2OUT
         write(*,*) "Recalc. W2::",VINT(4)**2+VINT(307)*(1./XCALC-1.)
         write(*,*) "Recalc. nu from W2:",
     &              (VINT(2)+VINT(307)-VINT(4)**2)/(2.*VINT(4))
         write(*,*) ""
      endif

***************standard output for event info***************************
C      For USERSET 0:
C      USER1 = sigma_dipole
C      USER2 = <Q_T>
C      USER3 = Ngrey 
C
C      Already filled with PYQREC (PyQM Recoil) if USERSET1:
C      USER1 = PYQREC_T
C      USER2 = |PYQREC|
C      USER3 = PYQREC(4)  all in TRF with z along gamma* 
C
C      USERSET 2:
C      USER1 = Q2SPLAT - Q2 of last aborted event if retry>0
C      USER2 = YYSPLAT - y of last aborted event if retry>0
C      USER3 = PYQREC(4)  in TRF 
C
C      USERSET 3: 
C      USER1 = Fermi-corrected W2
C      USER2 = W2 after correction/W2F - 1 
C      USER3 = # of iterations needed for correction
C
C      USERSET 4 (already filled with):
C      USER1-3 = x-z of Nuclear longitudinal axis for 3d nuclei
C
C      USERSET 5: USER1=Zremn, USER2=sigma_dipole,USER3 not used
C
C      USERSET 6: 
C      USER1 = D2-corrected W2
C      USER2 = W2 after correction/W2F - 1 
C      USER3 = # of iterations needed for correction
C
C      USERSET 7: Calculated elsewhere
C
C      USERSET 8: Particle 4-momenta sums in lab (collider) frame
C      USER1 = EOUT
C      USER2 = PZOUT (+z along ion_in)
C      USER3 = PTOUT
C
C      USERSET 9-11: Particle 4-momenta sums in ion rest frame
C      USER1 = EOUT
C      USER2 = PZOUT (-z along e_in)
C      USER3 = PTOUT
C
C      MDB 2017-07-01 Count wounded nucleons, n_g, <Q_T> by hand.
C      n_g only counted for Fixed target. 0.3 < beta < 0.7
C      <Q_T> assumes all ID=80000 are absorbed in target
C
      NWND = 0
      NWDCH  = 0
      IF (tracknr.LT.IT+2) STOP "DT_PYOUT FATAL ERROR: Too few tracks."
      DO I=2,IT+1
         IF(IDHKK(I).NE.2112 .AND. IDHKK(I).NE.2212) THEN
            WRITE(*,*) 'DT_PYOUT ERROR: IN-NUCLEON W/ IDHKK=',IDHKK(I)
         ELSEIF (ISTHKK(I).EQ.18 .OR. ISTHKK(I).EQ.12) THEN
            NWND = NWND + 1
            IF (IDHKK(I).EQ.2212) NWDCH = NWDCH + 1
         ELSEIF (ISTHKK(I).NE.14) THEN
            WRITE(*,*) 'DT_PYOUT ERROR: IN-NUCLEON W/ ISTHKK=',ISTHKK(I)
         ENDIF
      ENDDO
      IF (USERSET.EQ.0) THEN
         USER1 = SIGEFF
         USER2 = 0.0D0
         USER3 = 0.0D0
         DO I=IT+2,tracknr
            IF (ISTHKK(I).EQ.1) THEN
               USER2 = USER2 + DBLE(PYCHGE(IDHKK(I))/3)
               IF (ABS(PZTARG).LT.0.001.AND.PYCHGE(IDHKK(I)).EQ.3) THEN
                  P2TEMP = PHKK(4,I)*PHKK(4,I) - PHKK(5,I)*PHKK(5,I)
                  BETA2 = P2TEMP/PHKK(4,I)*PHKK(4,I)
                  IF (0.04.LT.P2TEMP .AND. P2TEMP.LT.0.36 .AND. 
     &                 BETA2.LT.0.49) USER3=USER3+1.0D0
               ENDIF
            ENDIF
         ENDDO
      ELSEIF (USERSET.EQ.2) THEN
         USER1 = Q2SPLAT
         USER2 = YYSPLAT
      ELSEIF (USERSET.EQ.5) THEN         
         USER2 = SIGEFF
      ELSEIF (USERSET.GE.8.AND.USERSET.LT.14) THEN         
         VERBOSE = (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5))
         AMASS1 = AZMASS(IT,ITZ,1)
         AMASS2 = AZMASS(IT,ITZ,2)
         AMASS3 = AZMASS(IT,ITZ,3)
         EINA1 = SQRT(PZTARG*PZTARG*IT*IT + AMASS1*AMASS1)
         EINA2 = SQRT(PZTARG*PZTARG*IT*IT + AMASS2*AMASS2)
         EINA3 = SQRT(PZTARG*PZTARG*IT*IT + AMASS3*AMASS3)
         IF (VERBOSE) THEN
            PZIN=PZLEP+IT*PZTARG
            WRITE(*,*) 'Check electron mass:',PYMASS(ltype)
            EINE = SQRT(PZLEP*PZLEP+PYMASS(ltype)*PYMASS(ltype))
            EIN1 = EINE+EINA1
            EIN2 = EINE+EINA2
            EIN3 = EINE+EINA3
            WRITE(*,*)'Lab frame quantities:'
            WRITE(*,*)'Naive AZMASS lookup: Ein, pzin:',EIN1,PZIN
            WRITE(*,*)'Corr. AZMASS lookup: Ein, pzin:',EIN2,PZIN
            WRITE(*,*)'DPMJET-F EXMSAZ:     Ein, pzin:',EIN3,PZIN
            WRITE(*,*)'Lab sums of stable particles:'
         ENDIF
         CALL EVTSUM(P5SUM,IZSUM,IASUM,VERBOSE,.FALSE.)
         USER3 = SQRT(P5SUM(1)*P5SUM(1)+P5SUM(2)*P5SUM(2))
         IF (USERSET.EQ.8) THEN
            USER1 = P5SUM(4)
            USER2 = P5SUM(3)
         ENDIF
         IF (USERSET.GT.8 .OR. VERBOSE) THEN
            GAMTMP1 = EINA1/AMASS1
            GAMTMP2 = EINA2/AMASS2
            GAMTMP3 = EINA3/AMASS3
            BGAMTMP1 = PZTARG*IT/AMASS1
            BGAMTMP2 = PZTARG*IT/AMASS2
            BGAMTMP3 = PZTARG*IT/AMASS3
            ESUM1 = GAMTMP1*P5SUM(4)-BGAMTMP1*P5SUM(3)
            ESUM2 = GAMTMP2*P5SUM(4)-BGAMTMP2*P5SUM(3)
            ESUM3 = GAMTMP3*P5SUM(4)-BGAMTMP3*P5SUM(3)
            PZSUM1 = GAMTMP1*P5SUM(3)-BGAMTMP1*P5SUM(4)
            PZSUM2 = GAMTMP2*P5SUM(3)-BGAMTMP2*P5SUM(4)
            PZSUM3 = GAMTMP3*P5SUM(3)-BGAMTMP3*P5SUM(4)
            IF (USERSET.EQ.9) THEN
               USER1=ESUM1
               USER2=PZSUM1
            ELSEIF (USERSET.EQ.10) THEN
               USER1=ESUM2
               USER2=PZSUM2
            ELSEIF (USERSET.EQ.11) THEN
               USER1=ESUM3
               USER2=PZSUM3
            ELSEIF (USERSET.EQ.12) THEN
               USER1=ESUM3
               USER2=IZSUM
               USER3=IASUM
            ELSEIF (USERSET.EQ.13) THEN
               USER1=NHYPER
               USER2=IZSUM
               IF (NHYPER.GE.1) THEN
                  USER3=IDHYP(1)
               ELSE
                  USER3=0
               ENDIF
            ENDIF
            IF (VERBOSE) THEN
               WRITE(*,*)'IRF frame (e- defines -z) quantities:'
               WRITE(*,*)'Naive AZMASS lookup Ein, pzin:  ',
     &            GAMTMP1*EIN1-BGAMTMP1*PZIN,GAMTMP1*PZIN-BGAMTMP1*EIN1
               WRITE(*,*)'Naive AZMASS lookup Eout, pzout:',ESUM1,PZSUM1
               WRITE(*,*)'Corr. AZMASS lookup Ein, pzin:  ',
     &            GAMTMP2*EIN2-BGAMTMP2*PZIN,GAMTMP2*PZIN-BGAMTMP2*EIN2
               WRITE(*,*)'Corr. AZMASS lookup Eout, pzout:',ESUM2,PZSUM2
               WRITE(*,*)'DPMJET-F EXMSAZ Ein, pzin:      ',
     &            GAMTMP3*EIN3-BGAMTMP3*PZIN,GAMTMP3*PZIN-BGAMTMP3*EIN3
               WRITE(*,*)'DPMJET-F EXMSAZ Eout, pzout:    ',ESUM3,PZSUM3
            ENDIF
         ENDIF
      ENDIF

      IF (NCOLLT.NE.NTW0)
     &  WRITE(*,*)'PYOUTEP ERROR: NCOLLT .NE. NTW0: ',NCOLLT,' ',NTW0
C      WRITE(*,*)
C      WRITE(*,*) 'Output XBJ= ',XBJOUT
C      WRITE(*,*) 'Output RAVAL= ',RAVAL
C      WRITE(*,*) 'Output RAEVT= ',RAEVT
C      WRITE(*,*)
      if (OLDOUT) then
C        Note: nrtrack is still # of tracks-4. So add 4.
         write(29,32) 0, ievent, genevent, ltype, it, itz, pzlep, 
     &        pztarg, pznucl, msti(1), msti(12), msti(16), pari(34), 
     &        msti(15), pari(33), pari(53), 
     &        YYOUT, Q2OUT, XBJOUT, W2OUT, NUOUT,
     &        VINT(313), pari(14), pari(15), pari(16), 
     &        pari(18),  pari(22), sngl(py6f2), sngl(py6f1), 
     &        py6r, mcRadCor_Sigrad, mcRadCor_sigcor, radgamEnucl,
     &        VINT(319), BBEA, PHIB, THKB, THKSCL, NCOLLT, NCOLLI,
     &        NNEVAP, NPEVAP, AREMN, NINC, NINCCH, DFIRST, DAVG,
     &        PXF, PYF, PZF, EEXC(2), RAEVT, USER1, USER2, USER3, 
     &        nrtrack+4 
      else
         write(29,33) 0, ievent, genevent, ltype, it, itz, pzlep, 
     &        pztarg, pznucl, crang, crori, msti(1), msti(12),
     &        msti(16), pari(34), msti(15), pari(33), pari(53), 
     &        YYOUT, Q2OUT, XBJOUT, W2OUT, NUOUT,
     &        VINT(313), pari(14), pari(15), pari(16), 
     &        pari(18),  pari(22), sngl(py6f2), sngl(py6f1), 
     &        py6r, mcRadCor_Sigrad, mcRadCor_sigcor, radgamEnucl,
     &        VINT(319), BBEA, PHIB, THKB, THKSCL, NCOLLT, NCOLLI,
     &        NWND, NWDCH,
     &        NNEVAP, NPEVAP, AREMN, NINC, NINCCH, DFIRST, DAVG,
     &        PXF, PYF, PZF, EEXC(2), RAEVT, USER1, USER2, USER3, 
     &        nrtrack+4 
      endif
C 32   format((I4,1x,$),(I10,1x,$),3(I4,1x,$),(I10,1x,$),f9.6,1x,$,
C     &     I12,1x,$,2(f12.6,1x,$),7(f18.11,3x,$),11(f19.9,3x,$),I12,/)
 32   format((I4,1x,$),(I10,1x,$),4(I4,1x,$),3(f12.6,1x,$),2(I4,1x,$),
     &     I6,1x,$,f9.6,1x,$,I6,1x,$,2(f12.6,1x,$),7(f18.11,3x,$),
     &     11(f19.9,3x,$),4(f10.6,1x,$),7(I5,1x,$),2(f10.6,1x,$),
     &     3(f15.6,1x,$),2(f9.6,1x,$),3(e17.8,1x,$),I6,/)
 33   format((I4,1x,$),(I10,1x,$),4(I4,1x,$),4(f12.6,1x,$),3(I4,1x,$),
     &     I6,1x,$,f9.6,1x,$,I6,1x,$,2(f12.6,1x,$),7(f18.11,3x,$),
     &     11(f19.9,3x,$),4(f10.6,1x,$),9(I5,1x,$),2(f10.6,1x,$),
     &     3(f15.6,1x,$),f12.6,1x,$,f9.6,1x,$,3(e17.8,1x,$),I6,/)
      write(29,*)'============================================'

***************standard output for particle info************************
c...add 2 beam information at first to fit into root tree making rule      
c... MDB Change these lines to use the correct pid for lepton+nucleon!
      I=NPOINT(1)+1   
C      write(29,34) 1,21,11,0,0,I+4,0,
      write(29,34) 1,21,IDHKK(I),0,0,I+4,0,
     &     PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),PHKK(5,I),
     &     VHKK(1,I),VHKK(2,I),VHKK(3,I)
     &     ,0,0,0
c...nuclear beam (if nucleus) needs to be made with modification
      PP1=0
      PP2=0
c         IF(INUMOD.GT.1) THEN
c            PP5 = PINITA(5)
c            call DT_LTNUC(PINITA(3),PINITA(4),PP3,PP4,-3)
c            !!Lorentz Transformation
c            !!E = gama*(E-beta*pz); pz = gamma*(pz-beta*E)
c            P4 = pgamma*PP4 - pgamma*(-pbeta)*PP3
c            P3 = -pgamma*(-pbeta)*PP4 + pgamma*PP3
c         ELSE
      PP5 = PHKK(5,I+1)
      P4 = PHKK(4,I+1)
      P3 = PHKK(3,I+1)
c         ENDIF
      IF (OLDOUT) THEN
         AOUT = 0
         ZOUT = 0
      ELSE 
         AOUT = 1
C         ZOUT = 1 
         ZOUT = PYCHGE(IDHKK(I+1))/3
      ENDIF
C      write(29,34) 2,21,2212,0,0,I+5,0,PP1,PP2,P3,P4,PP5,
      write(29,34) 2,21,IDHKK(I+1),0,0,I+5,0,PP1,PP2,P3,P4,PP5,
     &        VHKK(1,I+1),VHKK(2,I+1),VHKK(3,I+1),AOUT,ZOUT,0
c...add the lepton 
C      write(29,34) 3,21,11,0,1,ILEPT+4,0,
      write(29,34) 3,21,IDHKK(ILEPT),0,1,ILEPT+4,0,
     &     PHKK(1,ILEPT),PHKK(2,ILEPT),PHKK(3,ILEPT),
     &     PHKK(4,ILEPT),PHKK(5,ILEPT),VHKK(1,ILEPT),
     &     VHKK(2,ILEPT),VHKK(3,ILEPT)
     &     ,0,0,0
c...add the exchanged boson  
C      IDEXBO = 22
C      IF (MCGENE.EQ.6) IDEXBO = 23
      write(29,34) 4,21,IDHKK(IBOSON),0,1,IBOSON+4,0,
     &     PHKK(1,IBOSON),PHKK(2,IBOSON),PHKK(3,IBOSON),
     &     PHKK(4,IBOSON),PHKK(5,IBOSON),VHKK(1,IBOSON),
     &     VHKK(2,IBOSON),VHKK(3,IBOSON)
     &     ,0,0,0
 
      DO I=1,tracknr
c         if (K(I,3).le.nrtrack) then
c...make the mother daughter relation consistent with 2 beam particles
c...and virtual photon added on   
         JM1OUT = 0
         JM2OUT = 0
         JD1OUT = 0
         JD2OUT = 0
         AOUT = 0
         ZOUT = 0
c 2016-12-30 MDB Don't actually change JMOHKK & JDAHKK since we
C                didn't change PHKK etc.
C                Also fix bug where Mother2 wasn't offset
         IF(I.NE.ILEPT) THEN
C               IF(JMOHKK(1,I).GT.0) JMOHKK(1,I)=JMOHKK(1,I)+4
C               IF(JDAHKK(1,I).GT.0) JDAHKK(1,I)=JDAHKK(1,I)+4
C               IF(JDAHKK(2,I).GT.0) JDAHKK(2,I)=JDAHKK(2,I)+4
            IF(JMOHKK(1,I).GT.0) JM1OUT = JMOHKK(1,I)+4
            IF(JMOHKK(2,I).GT.0) JM2OUT = JMOHKK(2,I)+4
            IF(JDAHKK(1,I).GT.0) JD1OUT = JDAHKK(1,I)+4
            IF(JDAHKK(2,I).GT.0) JD2OUT = JDAHKK(2,I)+4
         ELSE
C     Special treatment for scattered lepton
            JM1OUT = 3
         ENDIF
         KSOUT = ISTHKK(I)
         BAMOUT = NOBAM(I)
         IDOUT = IDHKK(I)
         IF (IDHKK(I).EQ.80000) THEN
            ZOUT = IDXRES(I)
            AOUT = IDRES(I)
            IF (IOUTFMT.EQ.2) THEN
               IDOUT = 1000000000 + 10000*ZOUT + 10*AOUT
            ENDIF
         ELSEIF (.NOT. OLDOUT) THEN
            AOUT = NBARY(IDHKK(I))
            IF (MOD(AOUT,3).EQ.0) THEN
               AOUT = AOUT/3
            ELSE
               AOUT = 0
            ENDIF
            ZOUT = PYCHGE(IDHKK(I))
            IF (MOD(ZOUT,3).EQ.0) THEN
               ZOUT = ZOUT/3
            ELSE
               ZOUT = 0
            ENDIF
         ENDIF
         IF (.NOT.OLDOUT) THEN
!!!dump nuclear remnants into final state particles
            IF (ISTHKK(I).EQ.-1) THEN
               KSOUT = 1
               BAMOUT = 3
            ELSEIF (ISTHKK(I).EQ.1001) THEN
               KSOUT = 1
               BAMOUT = 4
            ENDIF
         ENDIF
C         write(29,34) I+4, KSOUT, IDHKK(I), JM2OUT, JM1OUT, 
         write(29,34) I+4, KSOUT, IDOUT, JM2OUT, JM1OUT, 
     &        JD1OUT, JD2OUT, PHKK(1,I), PHKK(2,I), PHKK(3,I),
     &        PHKK(4,I), PHKK(5,I), VHKK(1,I), VHKK(2,I), VHKK(3,I),
     &        AOUT, ZOUT, BAMOUT
         ENDDO
c         if (mcRadCor_EBrems.gt.0.) then
c            write(29,34) nrtrack, 55, 22, 1, 0, 0,
c     &      sngl(dplabg(1)),sngl(dplabg(2)),sngl(-radgamp),
c     &      sngl(radgamE), 0., 0., 0., 0.
c         endif
C 34      format(2(I6,1x,$),I10,1x,$,3(I8,1x,$),5(f15.6,1x,$),
C     &      3(e15.6,1x,$)/)
 34      format(2(I6,1x,$),I10,1x,$,4(I8,1x,$),5(f15.6,1x,$),
     &       3(e15.6,1x,$),3(I8,1x,$)/)

         write(29,*)'=============== Event finished ==============='

         lastgenevent=MYNGEN

c         print*,'output finished'
      RETURN

c...mode 3 is used to print the whole statistics information
3     CONTINUE
      CALL PYSTAT(1)
      CALL PYSTAT(4)

      WRITE(*,*)'The charm mass used is: ', PMAS(4,1)
         
C...Print the Pythia cross section which is needed to get an absolute
C   normalisation the number is in microbarns
      write(*,*)'==================================================='
      if (ITMODE.EQ.0) then
         write(*,*)'FROM PYTHIA (2nd half: ep collisions):'
      elseif (ITMODE.EQ.1) then
         write(*,*)'FROM PYTHIA (en collisions only):'
      elseif (ITMODE.EQ.2) then
         write(*,*)'FROM PYTHIA (ep collisions only):'
      elseif (ITMODE.EQ.3) then
         write(*,*)'FROM PYTHIA (just last few events due to PYINITs):'
      endif
      write(*,*)'Pythia total cross section normalisation:',
     &            pari(1)*1000, ' microbarn'
      write(*,*)'Total Number of generated events pythia', MSTI(5)
      write(*,*)'==================================================='
      write(*,*)'Actual totals from BeAGLE:'
      write(*,*)'Total Number of generated events', NEVENT
      write(*,*)'Total Number of trials', MYNGEN
      write(*,*)'==================================================='
      close(29)

C...Check pdf status       
      call PDFSTA
      
      RETURN

c...mode 4 is used to output the event list to screen in current frame
c...without a lot of reformatting and rearranging (Mark 08/17/2016)
c...using the old event header (Mark 01/02/2017)
4     CONTINUE

      genevent=MYNGEN-lastgenevent
      tracknr = NHKK
      if (mcRadCor_EBrems.gt.0.) then
         nrtrack=tracknr+1
      else
         nrtrack=tracknr
      endif

c...print a title for the event file - use formats from case 2
      write(*,*)' DUMP of /DTEVT1/'
      write(*,*)'============================================'
      write(*,31)
      write(*,*)'============================================'
      write(*,*)' NPOINT(1-4):'
      write(*,*) NPOINT(1),' ',NPOINT(2),' ',NPOINT(3),' ',NPOINT(4)
      write(*,*)'============================================'

***************standard output for event info***************************
      write(*,33) 0, ievent, genevent, ltype, it, itz, pzlep, 
     &     pztarg, pznucl, crang, crori, msti(1), msti(12),
     &     msti(16), pari(34), msti(15), pari(33), pari(53), 
     &     YYOUT, Q2OUT, XBJOUT, W2OUT, NUOUT,
     &     VINT(313), pari(14), pari(15), pari(16), 
     &     pari(18),  pari(22), sngl(py6f2), sngl(py6f1), 
     &     py6r, mcRadCor_Sigrad, mcRadCor_sigcor, radgamEnucl,
     &     VINT(319), BBEA, PHIB, THKB, THKSCL, NCOLLT, NCOLLI,
     &     NTW, NTCW,
     &     NNEVAP, NPEVAP, AREMN, NINC, NINCCH, DFIRST, DAVG,
     &     PXF, PYF, PZF, EEXC(2), RAEVT, USER1, USER2, USER3, 
     &     nrtrack 
      write(*,*)'============================================'

***************standard output for particle info************************
c...use the dpmjet track wide title - EXTENDED!
      write(*,*)'I  ISTHKK(I)  IDHKK(I)  JMOHKK(2,I)  JMOHKK(1,I)
     & JDAHKK(1,I)  JDAHKK(2,I)  PHKK(1,I)  PHKK(2,I)  PHKK(3,I)
     & PHKK(4,I)  PHKK(5,I)  VHKK(1,I) VHKK(2,I) VHKK(3,I) IDRES(I)
     & IDXRES(I)  NOBAM(I), IDBAM(I), IDCH(I)'
      write(*,*)'============================================'
      DO I=1,tracknr
         write(*,35) I,ISTHKK(I),IDHKK(I),JMOHKK(2,I),JMOHKK(1,I),
     &        JDAHKK(1,I),JDAHKK(2,I),PHKK(1,I),PHKK(2,I),PHKK(3,I),
     &        PHKK(4,I),PHKK(5,I),VHKK(1,I),VHKK(2,I),VHKK(3,I),
     &        IDRES(I),IDXRES(I),NOBAM(I),IDBAM(I),IDCH(I)
      ENDDO
      write(*,*)'=============== Event finished ==============='
c     35 is similar to 34, but with 2 extra integers at the end.
c     Shortened some fields 
 35   format(2(I4,1x,$),I15,1x,$,4(I5,1x,$),5(f13.6,1x,$),
     &       3(e14.6,1x,$),5(I4,1x,$)/)
      
      RETURN

      END 


*=====dt_pyf2qpm=========================================================
*used for the calculation of F2 for pythia events sampling in QPM
*formalism
      SUBROUTINE DT_PYF2QPM(Q2,X,F2) 

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
     
c...Local variables
      DIMENSION XPQ(-25:25)
      DOUBLE PRECISION Q2,F2,X

      INTEGER i,j,iter

      call PYPDFU(2212,X,Q2,XPQ)
      F2=4.0*(XPQ(2)+XPQ(-2)+XPQ(4)+XPQ(-4))/9.0+
     &   (XPQ(1)+XPQ(-1)+XPQ(3)+XPQ(-3))/9.0
   
      RETURN

      END    
       
