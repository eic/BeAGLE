*These are the interfaces to GCF input file needed by dpmjet
*
* Written by  Mark Baker  2019-04-27
* Updated                 2020-03-04
*
C...intialize GCF-reading code when using dpmjet
*======================================================================
      subroutine DT_GCFINITQE(INPUT,OUTPUT,IFSEED,NEVTS)
*     input variables:
*           INPUT    Input file with simulated GCF QE physics events.
*           OUTPUT   Output file.
*           IFSEED   Input seed for random #s. 0=use a random seed.
*     input/output variable:
*           NEVTS    # of events to read. If 0, then read all events
*                    If NEVTS># entries or NEVTS=0, set NEVTS=# entries
*
      IMPLICIT NONE
C      include 'pythia.inc'              ! All PYTHIA commons blocks
C      include "mc_set.inc"
C      include "py6strf.inc"
C      include "mcRadCor.inc"
C      include "radgen.inc"
C      include "phiout.inc"
      include "beagle.inc"
C      include "bea_pyqm.inc"

      CHARACTER*60 INPUT,OUTPUT
      INTEGER IFSEED,NEVTS

      DOUBLE PRECISION DT_RNCLUS
      EXTERNAL DT_RNCLUS

* properties of interacting particles

c...target/proj mass, charge and projectile internal ID
      integer IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,ITMMOD,
     &     MODHYP,NHYPER,IDHYP 

      double precision EPN,PPN
      DOUBLE PRECISION ZERO,ONE,TINY10
      INTEGER MAXNCL
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY10=1.0D-10,MAXNCL=260)

C...Added by liang 1/6/12
C...Switches for nuclear correction
      COMMON /PYNUCL/ INUMOD,CHANUM,ORDER,genShd
      SAVE /PYNUCL/
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER,genShd

CC...Pythia event counter (since we keep PYINITing) Mark 2017-01-31
C      COMMON /PYCNTR/ MYNGEN
C      INTEGER MYNGEN
C      SAVE /PYCNTR/

C...Parameters and switch for energy loss
      DOUBLE PRECISION QHAT
      INTEGER QSWITCH
      COMMON /QUENCH/ QHAT, QSWITCH

C... Locals
      integer NEV

      INTEGER NPRT, ievent, genevent, I, tracknr 
      INTEGER  lastgenevent, idum1, idum2, initseed, nrtrack
      REAL trueX, trueW2, trueNu
      DOUBLE PRECISION sqrts, radgamE, radgamp, radgamEnucl
      DOUBLE PRECISION pbeamE, ebeamE, epznucl 
      DOUBLE PRECISION altpbeamE, altpbeam, altsqrts
      CHARACTER PARAM*100
      LOGICAL UseLut, GenLut
      DOUBLE PRECISION SHDFAC
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
C
C      DOUBLE PRECISION Mprot,Mneut,Mnucl,Mlept
C      ! beam type
C      CHARACTER*10 tName
c ---------------------------------------------------------------------
c     Run parameter
c ---------------------------------------------------------------------
      integer*4 today(3), now(3)
c---------------------------------------------------------------------
c     ASCII output file and input file
c ---------------------------------------------------------------------
      integer LINP
      parameter ( LINP=28 )
      CHARACTER*256 inputfilename, outname
      COMMON /OUNAME/ outname

* 14-Dec-2016 MDB Event variable to find datafiles
      CHARACTER*255 ENVDIR
      COMMON /ENVCOM/ ENVDIR
      CHARACTER*255 FILNAM

* Glauber formalism: parameters
      INTEGER NCOMPX,NEB,NQB,KSITEB
      PARAMETER (NCOMPX=20,NEB=8,NQB= 5,KSITEB=50)
      COMMON /DTGLAM/ RASH(NCOMPX),RBSH(NCOMPX),
     &                BMAX(NCOMPX),BSTEP(NCOMPX),
     &                SIGSH,ROSH,GSH,BSITE(0:NEB,NQB,NCOMPX,KSITEB),
     &                NSITEB,NSTATB
      DOUBLE PRECISION RASH, RBSH, BMAX, BSTEP, SIGSH,ROSH,GSH,BSITE
      INTEGER NSITEB,NSTATB

* Locals
      CHARACTER*255 CDUMMY

      ievent=0
      genevent=0
      lastgenevent=0
      tracknr=0
c ---------------------------------------------------------------------
c     Open ascii input file
c ---------------------------------------------------------------------
      inputfilename=INPUT
      outname=OUTPUT
      OPEN(LINP, file=inputfilename,STATUS='UNKNOWN')
      WRITE(*,*) 'the input file is :', inputfilename
      READ(LINP,*) CDUMMY
      WRITE(*,*) 'First line is :',CDUMMY
      READ(LINP,*) NEV, IT, ITZ
      WRITE(*,*) 'NEV, IT, ITZ :',NEV, IT, ITZ
      IF (NEVTS.EQ.0 .OR. NEVTS.GT.NEV) NEVTS=NEV
      READ(LINP,*) CDUMMY
      READ(LINP,*) CDUMMY
      READ(LINP,*) CDUMMY
      READ(LINP,*) CDUMMY
      READ(LINP,*) CDUMMY
      genShd=1
      QSWITCH=0
      QHAT=0.0
      SHDFAC=1.0
      
c...read parameters from dpmjet       
C       INUMOD=IT
C       CHANUM=ITZ
C       mcSet_YMin=real(YMIN)
C       mcSet_YMax=real(YMAX)
C       mcSet_Q2Min=real(Q2MIN)
C       mcSet_Q2Max=real(Q2MAX)

      write(*,*) 'the output file is: ', outname
C       print*,'kinematics cut read by PYTHIA:'
C       print*,YMIN,' < y < ',YMAX,', ',Q2MIN,' < Q2 < ',Q2MAX

C     Getting the date and time of the event generation

      IF (IFSEED.NE.0) THEN
         initseed=IFSEED
      ELSE
         call getseed(initseed,.TRUE.)
      ENDIF
      write(6,*) 'SEED = ', initseed
      call rndmq (idum1,idum2,initseed,' ')
        
C      !default ltype, lName and tName
C      ltype = 11
C      lName = 'gamma/e-'
C      tName = 'p+'

C      !set up nucleon beam
C      Mneut=PYMASS(2112)
C      Mprot=PYMASS(2212)
C      Mlept=PYMASS(ltype)
C      masse=real(Mlept)
C Mark 2015-10-21 This is for the first PYINIT only.
C                 For eA we may have to re-PYINIT event by event,
C                 but start with a neutron since that's the most likely
C Mark 2018-01-24 Use double precision here.

C Note: I may need to setup IJTARG=1 for p or 8 for n

C      MAscl=AZMASS(IT,ITZ,ITMMOD)/INUMOD 

C     For now the input file is defined in the IRF/TRF
C     Nuclear beta=0.
C
      pbeta=0.0D0
      pgamma=1.0D0

C      ebeamE=sqrt(EPN**2+Mlept**2)
C      ebeamEnucl=pgamma*ebeamE-pgamma*pbeta*(-EPN)
C      mcSet_EneBeam=sngl(ebeamEnucl)

C     GenNucDens initializes nuclear spatial density distribution
      if (INRLEV.GT.2) call GenNucDens(ITZ, IT)

C     Initialize RBSH & BMAX parameters
      RASH(1) = 0.0D0
      RBSH(1) = DT_RNCLUS(IT)
      BMAX(1) = 2.0D0*(RASH(1)+RBSH(1))
C      !changed by liang to fit with INC particle status
C      CALL DT_PYDECY(1)

      RETURN
      END
*=====dt_gcfevntqe========================================================
* read in the input file from GCF and dump this event to DTEVT1
* GCF is in the target/ion rest frame, while dpmjet makes many of it's 
* calculations in the naive photon-proton c.m.s frame. 
* Unfortunately, this frame often doesn't exist (e.g. when xB>1).
* Plus, the virtual photon has to be directed to the z+ direction. A 
* rotation for the reference is also necessary.
*
      SUBROUTINE DT_GCFEVNTQE(Q2,YY,MODE,IREJ)
*
*     MDB 2020-02-14 Latest Version
* 
*     input: 
*           MODE  0 = simplified logic w/o INC. Copy GCF to DTEVT1
*           MODE  1 = read in GCF event to internal common & get Q2,Y
*           MODE  2 = Copy GCF event from internal common to DTEVT1
*     output:
*           Q2    Q2 of this current event 
*           YY    Y of this current event  
*           IREJ  Error flag 0=OK
* 
*     Variable OLDOUT from "beagle.inc" controls the input format.

      IMPLICIT NONE
C      include "mc_set.inc"
C      include "mcRadCor.inc"
C      include "radgen.inc"
C      include "phiout.inc"
      include "beagle.inc"
      include "bea_pyqm.inc"

      EXTERNAL IDT_ICIHAD, DCALC, SCLFAC, DT_RNDM
      INTEGER IDT_ICIHAD
      DOUBLE PRECISION DCALC, SCLFAC, DT_RNDM

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 
      INTEGER NEVENT, ICASCA

* event history
      INTEGER NMXHKK
      DOUBLE PRECISION FM2MM, MNGCF
      PARAMETER (NMXHKK=200000)
      PARAMETER (FM2MM=1.0D-12)
      PARAMETER (MNGCF=0.93892D0)

* properties of interacting particles
      INTEGER IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5) 

      INTEGER NHKK,NEVHKK,ISTHKK,IDHKK,JMOHKK,JDAHKK
      DOUBLE PRECISION PHKK, VHKK, WHKK
      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history - IHIST renamed because it's in pythia.inc
      INTEGER IDRES,IDXRES,NOBAM,IDBAM,IDCH,NPOINT,IHISTDPM
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHISTDPM(2,NMXHKK)

* treatment of residual nuclei: 4-momenta
      LOGICAL LRCLPR,LRCLTA
      DOUBLE PRECISION PINIPR,PINITA,PRCLPR,PRCLTA,TRCLPR,TRCLTA
      COMMON /DTRNU1/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5),
     &                TRCLPR(5),TRCLTA(5),LRCLPR,LRCLTA

* treatment of residual nuclei: properties of residual nuclei
      DOUBLE PRECISION AMRCL0,EEXC,EEXCFI
      INTEGER NTOT,NPRO,NN,NH,NHPOS,NQ,NTOTFI,NPROFI
      COMMON /DTRNU2/ AMRCL0(2),EEXC(2),EEXCFI(2),
     &                NTOT(2),NPRO(2),NN(2),NH(2),NHPOS(2),NQ(2),
     &                NTOTFI(2),NPROFI(2)

      DOUBLE PRECISION ZERO,ONE,TINY10
      INTEGER MAXNCL
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY10=1.0D-10,MAXNCL=260)

* Lorentz-parameters of the current interaction from DPMJET
      DOUBLE PRECISION GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PPCM,EPROJ,PPROJ
      COMMON /DTLTRA/ GACMS(2),BGCMS(2),GALAB,BGLAB,BLAB,
     &                UMO,PPCM,EPROJ,PPROJ

* lorentz transformation parameter 
      DOUBLE PRECISION BGTA,GAMM,eveBETA,GAA,COF,COD,SIF,SID
      COMMON /LTPARA/ BGTA(4), GAMM(5), eveBETA(4), GAA
      COMMON /ROTATE/ COF,COD,SIF,SID
      
* properties of photon/lepton projectiles from DPMJET
      DOUBLE PRECISION VIRT,PGAMM,PLEPT0,PLEPT1,PNUCL
      INTEGER IDIREC
      COMMON /DTGPRO/ VIRT,PGAMM(4),PLEPT0(4),PLEPT1(4),PNUCL(4),IDIREC

* kinematics at lepton-gamma vertex from DPMJET
      DOUBLE PRECISION PPL0,PPL1,PPG,PPA
      COMMON /DTLGVX/ PPL0(4),PPL1(4),PPG(4),PPA(4)

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

      DOUBLE PRECISION P1, P2, P3, P4, PTOT, PP1, PP2, PP3, PP4
      INTEGER QCHG

      DOUBLE PRECISION Q2, YY, XX, XXALT
      INTEGER MODE, IREJ

* added by liang to store the output event variables 1/20/12
      COMMON /EVTOUT/ XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT
      DOUBLE PRECISION XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT

C     Extensions for BEAGCF - variables that we would look up in a 
C     Pythia common blocks for MCGENE=5. Note: RECTYPE will be 
C     reported as process.
      INTEGER RECTYPE, LEADTYPE
      DOUBLE PRECISION LEPTONPHI
      COMMON /BEAGCF/ RECTYPE, LEADTYPE, LEPTONPHI

C     Temporary storage for GCF event until DTEVT1 is ready for it.
      INTEGER MAXGCF,NGCF,IOFF
      PARAMETER (MAXGCF=100)
      INTEGER ISTGCF, IDGCF, JMOGCF, JDAGCF, ARESGCF, ZRESGCF, NOBAMGCF
      INTEGER IDBAMGCF
      DOUBLE PRECISION PGCF
      COMMON /GCFDAT/ ISTGCF(MAXGCF),IDGCF(MAXGCF),JMOGCF(2,MAXGCF),
     &     JDAGCF(2,MAXGCF), PGCF(5,MAXGCF),ARESGCF(MAXGCF),
     &     ZRESGCF(MAXGCF),NOBAMGCF(MAXGCF),IDBAMGCF(MAXGCF)

C  Local

      INTEGER ITRK, NINTS, IIMAIN, IIMAINN, ILEPT, III, I5

      integer LINP
      parameter ( LINP=28 )
      
      INTEGER IDUM, IEVENT, IDIM, NDIM
      INTEGER nrTracks
      parameter ( NDIM=4 )
      
      DOUBLE PRECISION PFERR, PPSUM(NDIM), DISTSRC, DTRY, DTEMP
      DOUBLE PRECISION DPNUCL(3), ANORF, PPT
      LOGICAL FNDNEW, KOSHER
      INTEGER NUMNUC,INUC,NUCHIT,ISWAP,ISWTRY,ITEMP

      IREJ = 0
      IF (MODE.EQ.0) CALL DT_EVTINI

C     Event header common to  mode 0 and 1.
      IF (MODE.LT.2) THEN
         IF (OLDOUT) THEN
            READ(LINP,*) IDUM, IEVENT, ltype, IT, ITZ, PZLEP, RECTYPE, 
     &         LEADTYPE, Q2OUT, XBJOUT, NUOUT, LEPTONPHI, PXF, PYF, PZF,
     &         EEXC(2), RAEVT, nrTracks
            IF (IEVENT.LE.5) WRITE(*,*) 'weight: ', RAEVT
         ELSE
            READ(LINP,*) IDUM, IEVENT, ltype, IT, ITZ, PZLEP, RECTYPE, 
     &         LEADTYPE, Q2OUT, XBJOUT, NUOUT, LEPTONPHI, PXF, PYF, PZF,
     &         EEXC(2), RAEVT, SIGEFF, nrTracks
            IF (IEVENT.LE.5) WRITE(*,*) 'lcweight: ', RAEVT, 
     &           ' oldweight',SIGEFF
         ENDIF
         IF (IDUM.NE.0 .OR. nrTracks.NE.9 .OR. EEXC(2).LT.-TINY10 .OR.
     &        IEVENT.NE.NEVENT) THEN
            WRITE(*,*) "IDUM, nrTracks, EEXC(2), IEVENT, NEVENT: ", 
     &           IDUM, nrTracks, EEXC(2), IEVENT, NEVENT
            STOP
     &        'DT_GCFEVNTQE: Bad event header format in GCF input file.'
         ENDIF
C     Protect against slightly negative roundoff errors for E*=0.
         IF (EEXC(2).LT.0.0D0) EEXC(2)=0.0D0
         YYOUT = NUOUT/ABS(PZLEP)
C     W2 doesn't make physical sense for QE. Just M2 for pF=0.
C     Naively: W2OUT = 2.0D0*MNGCF*NUOUT - Q2OUT + MNGCF*MNGCF
         W2OUT = 0.0D0
      ENDIF

C...  Collector to get the momentum and charge sums 
C...  Note: for now, treat the (A-2)* as stable.
      DO IDIM=1,NDIM
         PPSUM(IDIM)=0.0D0
      ENDDO
      QCHG=0

      IF (MODE.GT.0) GOTO(1,2) MODE
      IF (MODE.NE.0) STOP "FATAL ERROR GCFEVNTQE: Illegal mode called"

C     Mode=0 - Read GCF directly into DTEVT1

      DO ITRK = 1, nrTracks
         READ(LINP,*) IDUM, ISTHKK(ITRK), IDHKK(ITRK), JMOHKK(2,ITRK),
     &        JMOHKK(1,ITRK), JDAHKK(1,ITRK), JDAHKK(2,ITRK), 
     &        PHKK(1,ITRK), PHKK(2,ITRK), PHKK(3,ITRK), PHKK(4,ITRK),
     &        PHKK(5,ITRK), IDRES(ITRK), IDXRES(ITRK),  NOBAM(ITRK)
         IF (IDUM.NE.ITRK) STOP
     &        'DT_GCFEVNTQE: Bad track format in GCF input file.'
         IF(NEVENT.LE.5) THEN
            WRITE(*,*) 'Track:',
     &           ITRK,ISTHKK(ITRK),IDHKK(ITRK),JMOHKK(1,ITRK),
     &           JMOHKK(2,ITRK),JDAHKK(1,ITRK),PHKK(1,ITRK),PHKK(2,ITRK)
     &           ,PHKK(3,ITRK),PHKK(4,ITRK),PHKK(5,ITRK),IDRES(ITRK),
     &           IDXRES(ITRK),NOBAM(ITRK)
         ENDIF
         DO IDIM = 1,NDIM
            VHKK(IDIM,ITRK)=0.0D0
            WHKK(IDIM,ITRK)=0.0D0
            IF (ISTHKK(ITRK).EQ.1 .OR. ISTHKK(ITRK).EQ.1000)
     &           PPSUM(IDIM) = PPSUM(IDIM)+PHKK(IDIM,ITRK)
         ENDDO
         IF (ISTHKK(ITRK).EQ.1 .OR. ISTHKK(ITRK).EQ.1000) 
     &        QCHG=QCHG+IDXRES(ITRK)
c...set BAM ID for the particles         
         IDBAM(ITRK)=IDT_ICIHAD(IDHKK(ITRK))
      ENDDO
      NHKK = nrTracks

C     MDB 2017-02-22
C     DT_FOZOCA, DT_SCN4BA need NPOINT(1) pointing to the last nucleon
C     DT_FOZOCA, DT_RESNCL, DT_SCN4BA: NPOINT(4) points to the 1st produced particle
C     DT_CHASTA: NPOINT(3) points to the 1st produced particle.
C
      NPOINT(1)=3
      NPOINT(2)=NPOINT(1)
      NPOINT(3)=NPOINT(1)+4
      NPOINT(4)=NPOINT(1)+4

********flag the two interacting nucleons
      NINTS=2
      IMAIN=1
      IINTER(1)=1
      IINTER(2)=2
      IIMAIN = IINTER(IMAIN)

C      NPOS(J)=0
C      PosAlt(J,KK)=0.0D0

      PosNuc(1)=VHKK(1,IIMAIN)
      PosNuc(2)=VHKK(2,IIMAIN)
      PosNuc(3)=VHKK(3,IIMAIN)
      PosNuc(4)=VHKK(4,IIMAIN)
      
      PFERR = ABS(PXF-PHKK(1,IIMAIN))+ABS(PYF-PHKK(2,IIMAIN))+
     &     ABS(PZF-PHKK(3,IIMAIN))
      IF (PFERR.GT.TINY10) 
     &     STOP 'DT_GCFEVNTQE: FATAL ERROR PXF,PYF,PZF MISMATCH'
      EKF = PHKK(4,IIMAIN)-PHKK(5,IIMAIN)

      THKB = 0.0D0
      THKSCL = 0.0D0
      DAVG = 0.0D0
      DFIRST = 0.0D0

c      write(*,*) 'DPMJET position',PosNuc(1),PosNuc(2),PosNuc(3)
      PYQREC(1)=0.0D0
      PYQREC(2)=0.0D0
      PYQREC(3)=0.0D0
      PYQREC(4)=0.0D0

C     Note: Usersets 0-4,6,7,13 are N/A. 5=Zremn only.
C     Usersets 8-11 are identical since lab=IRF. 
C     Userset 12 (EOUT, ZSUM, ASUM) is nautrally applicable.
C
C     In any case, initialize to zero.
      USER1=0.0D0
      USER2=0.0D0
      USER3=0.0D0


      IF (IDHKK(5).NE.22) STOP 'DT_GCFEVNTQE: Bad event format'

      GOTO 9999

 1    CONTINUE

C     Return Q2 and YY
      Q2 = Q2OUT
      YY = YYOUT

      IF (nrTracks.GT.MAXGCF) STOP 
     &     'FATAL ERROR in GCFEVNTQE: Too many tracks in GCF'
      DO ITRK = 1, nrTracks
         READ(LINP,*) IDUM, ISTGCF(ITRK), IDGCF(ITRK), JMOGCF(2,ITRK),
     &        JMOGCF(1,ITRK), JDAGCF(1,ITRK), JDAGCF(2,ITRK), 
     &        PGCF(1,ITRK), PGCF(2,ITRK), PGCF(3,ITRK), PGCF(4,ITRK),
     &        PGCF(5,ITRK), ARESGCF(ITRK), ZRESGCF(ITRK), NOBAMGCF(ITRK)
         IF (IDUM.NE.ITRK) STOP
     &        'DT_GCFEVNTQE: Bad track format in GCF input file.'
         IF(NEVENT.LE.5) THEN
            WRITE(*,*) 'Track:',
     &           ITRK,ISTGCF(ITRK),IDGCF(ITRK),JMOGCF(1,ITRK),
     &           JMOGCF(2,ITRK),JDAGCF(1,ITRK),PGCF(1,ITRK),PGCF(2,ITRK)
     &           ,PGCF(3,ITRK),PGCF(4,ITRK),PGCF(5,ITRK),ARESGCF(ITRK),
     &           ZRESGCF(ITRK),NOBAMGCF(ITRK)
         ENDIF
         IF (ISTGCF(ITRK).EQ.1 .OR. ISTGCF(ITRK).EQ.1000) THEN
            DO IDIM = 1,NDIM
               PPSUM(IDIM) = PPSUM(IDIM)+PGCF(IDIM,ITRK)
            ENDDO
            QCHG=QCHG+ZRESGCF(ITRK)
         ENDIF
c...set BAM ID for the particles         
         IDBAMGCF(ITRK)=IDT_ICIHAD(IDGCF(ITRK))
      ENDDO
      NGCF = nrTracks

C     Find gamma* angle in IRF-B=LAB frame.
C     Note: this is different than dpm_pythia where COF,SIF,COD,SID
C           refer to the gamma* angle in the HCMS-B frame.
      IF (IDGCF(5).NE.22) STOP 'DT_GCFEVNTQE: Bad event format'
      PTOT=SQRT(PGCF(1,5)**2+PGCF(2,5)**2+PGCF(3,5)**2)
      COD = PGCF(3,5)/PTOT
      PPT = SQRT(PGCF(1,5)**2+PGCF(2,5)**2)
      SID = PPT/PTOT
      IF(PGCF(1,5).GT.ZERO) THEN
         COF = ONE
      ELSE
         COF = -ONE
      ENDIF
      SIF = ZERO      
      IF (PTOT*SID.GT.TINY10) THEN
         COF = PGCF(1,5)/(SID*PTOT)
         SIF = PGCF(2,5)/(SID*PTOT)
         ANORF = SQRT(COF*COF+SIF*SIF)
         COF = COF/ANORF
         SIF = SIF/ANORF
      ENDIF

      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5) .AND. MODE.LT.2) then
         print*
         print*,'DT_GCFEVNTQE mode=',MODE
         print*,'4-momentum & charge totals for Status=1 & 1000:'
         print*,'PP1=',PPSUM(1),' PP2=',PPSUM(2),' PP3=',PPSUM(3),
     &          ' PP4=',PPSUM(4),' QCHG=',QCHG
         print*
         print*,'Rotate to IRF w/ z along gamma* and recalculate:'
         print*
      endif

      DO IDIM = 1,NDIM
         PPSUM(IDIM) = 0.0D0
      ENDDO

C     Rotate the GCF event so that it is in IRF-G
      DO ITRK=1,NGCF
         P1 = PGCF(1,ITRK)
         P2 = PGCF(2,ITRK)
         P3 = PGCF(3,ITRK)
         PGCF(1,ITRK)=COD*(P1*COF+P2*SIF)-P3*SID
         PGCF(2,ITRK)=-P1*SIF+P2*COF
         PGCF(3,ITRK)=SID*(P1*COF+P2*SIF)+P3*COD
         IF (ISTGCF(ITRK).EQ.1 .OR. ISTGCF(ITRK).EQ.1000) THEN
            DO IDIM = 1,NDIM
               PPSUM(IDIM) = PPSUM(IDIM)+PGCF(IDIM,ITRK)
            ENDDO
         ENDIF
      ENDDO

C     Now replace the original incoming nucleon pair w/ that from GCF
C     Change the GCF pair status to history 3 rather than 12
C     GCF versions of (A-2)* are event history, but not final
      ISTGCF(1)=3
      ISTGCF(2)=3
      ISTGCF(3)=3
      ISTGCF(9)=3

C... Flag the identity of the struck nucleon
      idNucPY=LEADTYPE
      idNucBAM=IDT_ICIHAD(idNucPY)

      GOTO 9999

 2    CONTINUE

C     MDB 2017-02-22
C     DT_FOZOCA, DT_SCN4BA need NPOINT(1) pointing to the last nucleon
C     DT_FOZOCA, DT_RESNCL, DT_SCN4BA: NPOINT(4) points to the 1st produced particle
C     DT_CHASTA: NPOINT(3) points to the 1st produced particle.
C
      NPOINT(1)=NHKK
      NPOINT(2)=NPOINT(1)
      NPOINT(3)=NPOINT(1)+5
      NPOINT(4)=NPOINT(1)+5

C     If the struck nucleon is the wrong flavor, repick
C     Note: interacting nucleons have ISTHKK=12, others 14.
C     Note: ITRK=1 refers to the gamma*, so we can skip it.
      NINTS=0
      FNDNEW=.FALSE.
      DO ITRK=2,NHKK
         IF (ISTHKK(ITRK).EQ.12) THEN
            IF (NINTS.GT.0) STOP "FATAL ERROR IN DT_GCFEVNTQE. "//
     &           "TOO MANY INTERACTIONS."
            NINTS=NINTS+1
            FNDNEW = (IDHKK(ITRK).NE.idNucPY)
            IF (FNDNEW) THEN
               ISTHKK(ITRK)=14
            ELSE
               IMAIN=1
               IINTER(1)=ITRK
               IIMAIN = ITRK
            ENDIF
         ENDIF
      ENDDO
      IF (NINTS.NE.1) STOP "FATAL ERROR IN DT_GCFEVNTQE. "//
     &           "NO INTERACTION FOUND."
      IF (FNDNEW) THEN
         IF (idNucPY.EQ.2112) THEN
            NUMNUC=ITZ
         ELSE
            NUMNUC=IT-ITZ
         ENDIF
         NUCHIT = INT(NUMNUC*DT_RNDM(IDUM))+1 
         INUC = 0
         NINTS = 0
         DO ITRK=2,NHKK
            IF (IDHKK(ITRK).EQ.idNucPY .AND. ISTHKK(ITRK).EQ.14) THEN
               INUC=INUC+1
               IF (INUC.EQ.NUCHIT) THEN
                  ISTHKK(ITRK)=12
                  NINTS = NINTS+1
                  IMAIN=1
                  IINTER(1)=ITRK
                  IIMAIN = ITRK
               ENDIF
            ENDIF
         ENDDO
         IF (INUC.NE.NUMNUC) STOP "FATAL ERROR IN DT_GCFEVNTQE. "//
     &           "WRONG # OF NUCLEONS."
         IF (NINTS.NE.1) STOP "FATAL ERROR IN DT_GCFEVNTQE. "//
     &           "WRONG # OF INTERACTIONS."
      ENDIF

C Find the closest nucleon to the struck nucleon
C     IIMAIN is the track # of the incoming struck nucleon
C     IIMAINN is the track # of the incoming recoil nucleon
C     DISTSRC is distance-squared in mm^2 (a very small #)
      DISTSRC=10000.0
      IIMAINN=-1
      DO ITRK=2,NHKK
         IF (ITRK.NE.IIMAIN) THEN
            DTRY = (VHKK(1,IIMAIN)-VHKK(1,ITRK))**2+
     &             (VHKK(2,IIMAIN)-VHKK(2,ITRK))**2+
     &             (VHKK(3,IIMAIN)-VHKK(3,ITRK))**2
            IF (DTRY.LT.DISTSRC) THEN
               IIMAINN = ITRK
               DISTSRC=DTRY
            ENDIF
         ENDIF
      ENDDO
C     If the nearest neighbor does not have the right flavor,
C     swap positions with one which does...
      IF (IDHKK(IIMAINN).NE.RECTYPE) THEN
c... Pick a # between 2 and NHKK and then look for a nucleon to swap flavors
         ISWAP=INT((NHKK-1)*DT_RNDM(IDUM))+2
         ISWTRY=0
         DO WHILE (ISWTRY.LT.NHKK-1)            
            IF ((IDHKK(ISWAP).EQ.RECTYPE) .AND.(ISTHKK(ISWAP).EQ.14))
     &           GOTO 25
            ISWTRY=ISWTRY+1
            ISWAP = ISWAP+1
            IF (ISWAP.GT.NHKK) ISWAP=2
         ENDDO
         STOP "DT_GCFEVNTQE FATAL ERROR: Can't find a nucleon to swap!"
 25      CONTINUE

         IDHKK(ISWAP)=IDHKK(IIMAINN)
         IDHKK(IIMAINN)=RECTYPE
         DO I5=1,5 
            DTEMP = PHKK(I5,IIMAINN)
            PHKK(I5,IIMAINN)=PHKK(I5,ISWAP)
            PHKK(I5,ISWAP) = DTEMP
         ENDDO
         ITEMP = NOBAM(IIMAINN)
         NOBAM(IIMAINN)=NOBAM(ISWAP)
         NOBAM(ISWAP)=ITEMP
         ITEMP = IDBAM(IIMAINN)
         IDBAM(IIMAINN)=IDBAM(ISWAP)
         IDBAM(ISWAP)=ITEMP
         ITEMP = IDCH(IIMAINN)
         IDCH(IIMAINN)=IDCH(ISWAP)
         IDCH(ISWAP)=ITEMP
      ENDIF
      ISTHKK(IIMAINN) = 12
      NINTS = NINTS + 1
      IINTER(2)=IIMAINN

C     Initial (A-2) momentum is -Pstruck-Pspec
C     We want it to be P(A-2) from GCF
C     DPNUCL is the 3-momentum to add to each A-2 incoming nucleon
C     Warning: Due to the Light-cone kinematics, the sum of all
C     A incoming nucleons will no longer have pz=0. This is expected.
      IF (NHKK.NE.IT+1 .OR. IT.LT.3)
     &     STOP "FATAL ERROR IN DT_GCFEVNTQE: INVALID STRUCTURE"
      DO IDIM=1,3
         DPNUCL(IDIM) = (PHKK(IDIM,IIMAIN)+PHKK(IDIM,IIMAINN)+
     &        PGCF(IDIM,3))/DBLE(IT-2)
      ENDDO
      
      IF (IDHKK(IIMAIN).NE.IDGCF(1))
     &     STOP ("FATAL ERROR IN DT_GCFEVNTQE: IIMAIN/GCF(1) mismatch")
      JMOGCF(1,1) = IIMAIN - NPOINT(1)
      JDAHKK(1,IIMAIN)=NPOINT(1)+1
      JDAHKK(2,IIMAIN)=0
      DO IDIM=1,5
         PHKK(IDIM,IIMAIN)=PGCF(IDIM,1)
      ENDDO
      IDXRES(IIMAIN)=ARESGCF(1)
      IDRES(IIMAIN)=ZRESGCF(1)
      NOBAM(IIMAIN)=NOBAMGCF(1)
C
      IF (IDHKK(IIMAINN).NE.IDGCF(2))
     &     STOP ("FATAL ERROR IN DT_GCFEVNTQE: IIMAINN/GCF(2) mismatch")
      JMOGCF(1,2) = IIMAINN - NPOINT(1)
      JDAHKK(1,IIMAINN)=NPOINT(1)+2
      JDAHKK(2,IIMAINN)=0
      DO IDIM=1,5
         PHKK(IDIM,IIMAINN)=PGCF(IDIM,2)
      ENDDO
      IDXRES(IIMAINN)=ARESGCF(2)
      IDRES(IIMAIN)=ZRESGCF(2)
      NOBAM(IIMAINN)=NOBAMGCF(2)
C
C     Make sure that (A-2)* will have the right 3-momentum
      DO ITRK=2,NHKK
         IF (ITRK.NE.IIMAIN .AND. ITRK.NE.IIMAINN) THEN
            DO IDIM=1,3
               PHKK(IDIM,ITRK)=PHKK(IDIM,ITRK)+DPNUCL(IDIM)
            ENDDO
            PHKK(4,ITRK)=SQRT(PHKK(1,ITRK)**2+PHKK(2,ITRK)**2+
     &           PHKK(3,ITRK)**2+PHKK(5,ITRK)**2)
         ENDIF
      ENDDO

C      Alternate definitions of pF...
C      PXF = -PHKK(1,IIMAINN)
C      PYF = -PHKK(2,IIMAINN)
C      PZF = -PHKK(3,IIMAINN)
C      EKF = PHKK(4,IIMAINN)-PHKK(5,IIMAINN)
C
      PXF = PHKK(1,IIMAIN)
      PYF = PHKK(2,IIMAIN)
      PZF = PHKK(3,IIMAIN)
      EKF = PHKK(4,IIMAIN)-PHKK(5,IIMAIN)

C      NPOS(J)=0
C      PosAlt(J,KK)=0.0D0

      PosNuc(1)=VHKK(1,IIMAIN)
      PosNuc(2)=VHKK(2,IIMAIN)
      PosNuc(3)=VHKK(3,IIMAIN)
      PosNuc(4)=VHKK(4,IIMAIN)

C Thickness is twice the integral from 0 to infinity
      THKB = 2.0D0*DCALC(0.0D0)
      IDUM = 0
      THKSCL = THKB*SCLFAC(IDUM)
C     Calculate distance travelled in the nucleus
      IF (NINTS.EQ.2) THEN
         ZMIN = 10000.0
         DO III=1,NINTS       
            ZTRY=VHKK(3,IINTER(III))
            IF (ZTRY.LT.ZMIN) THEN
               ZMIN=ZTRY
            ENDIF
         ENDDO
         DFIRST = DCALC(ZMIN/FM2MM)
         DAVG = 0.5D0*( DCALC(VHKK(3,IIMAIN)/FM2MM)+ 
     &                  DCALC(VHKK(3,IIMAINN)/FM2MM) )
      ELSE
         STOP "DT_GCFEVNTQE FATAL ERROR: NINTS.NE.2"
      ENDIF
      PYQREC(1)=0.0D0
      PYQREC(2)=0.0D0
      PYQREC(3)=0.0D0
      PYQREC(4)=0.0D0

C     Note: Usersets 0-4,6,7,13 are N/A. 5=Zremn only.
C     Usersets 8-11 are identical since lab=IRF. 
C     Userset 12 (EOUT, ZSUM, ASUM) is nautrally applicable.
C
C     In any case, initialize to zero.
      USER1=0.0D0
      USER2=0.0D0
      USER3=0.0D0

      IOFF=NPOINT(1)
      NHKK=NHKK+NGCF
      DO ITRK=1,NGCF
         ISTHKK(ITRK+IOFF)=ISTGCF(ITRK)
         IDHKK(ITRK+IOFF)=IDGCF(ITRK)
         JMOHKK(1,ITRK+IOFF)=JMOGCF(1,ITRK)+IOFF
         JMOHKK(2,ITRK+IOFF)=JMOGCF(2,ITRK)+IOFF
         JDAHKK(1,ITRK+IOFF)=JDAGCF(1,ITRK)+IOFF
         JDAHKK(2,ITRK+IOFF)=JDAGCF(2,ITRK)+IOFF
         IF (JMOHKK(1,ITRK+IOFF).EQ.IOFF) JMOHKK(1,ITRK+IOFF)=0
         IF (JMOHKK(2,ITRK+IOFF).EQ.IOFF) JMOHKK(2,ITRK+IOFF)=0
         IF (JDAHKK(1,ITRK+IOFF).EQ.IOFF) JDAHKK(1,ITRK+IOFF)=0
         IF (JDAHKK(2,ITRK+IOFF).EQ.IOFF) JDAHKK(2,ITRK+IOFF)=0
         IDRES(ITRK+IOFF)=ARESGCF(ITRK)
         IDXRES(ITRK+IOFF)=ZRESGCF(ITRK)
C     Requirements for INC routine (DT_FOZOCA): e' should be ISTHKK=99
C     Real particles must be NOBAM=0. Expect incoming e @ 4, e' @6
         IF (ITRK.EQ.6) THEN
            KOSHER = (
     &           (ABS(IDGCF(ITRK)).EQ.11 .OR. ABS(IDGCF(ITRK)).EQ.13)
     &           .AND. JMOGCF(1,ITRK).EQ.4 .AND. IDGCF(6).EQ.IDGCF(4)
     &           .AND. JMOGCF(1,4).EQ.0)
            IF (.NOT. KOSHER) STOP(
     &        "DT_GCFEVNTQE FATAL ERROR: Scattered lepton out of place")
            NOBAM(ITRK+IOFF)=0
            ISTHKK(ITRK+IOFF)=99
         ELSEIF (ISTGCF(ITRK).EQ.1) THEN
            NOBAM(ITRK+IOFF)=0
         ELSE
            NOBAM(ITRK+IOFF)=NOBAMGCF(ITRK)
         ENDIF
         IDBAM(ITRK+IOFF)=IDBAMGCF(ITRK)
         DO IDIM=1,5
            PHKK(IDIM,ITRK+IOFF)=PGCF(IDIM,ITRK)
         ENDDO
      ENDDO

c     WRITE(*,*) 'Calculated D1st, Davg: ',DFIRST,' ',DAVG

C    The lines below only come into play when we want to use the HCMS
C... Note: DPF(mu) = P(mu)_true - P(mu)_naive is a 4-momentum too.   
C    DPF is the name in the HCMS
C    PXF,PYF,PZF,EKF in the TRF
      DPF(1) = PXF
      DPF(2) = PYF
      DPF(3) = PZF
      DPF(4) = EKF
C      CALL DT_LTNUC(PZF,EKF,DPF(3),DPF(4),3)

C  The GCF event is already in IRF/TRF-G and .
C  For GCF, we will use the convention that -z is along the incoming electron
C  rather than the convention we used for Pythia that +z is along gamma*.

C      Boost from TRF z=g* (P3,P4)  to HCMS z=g* (PYQREC):
C      PYQREC(3)=GACMS(2)*P3-BGCMS(2)*P4
C      PYQREC(4)=GACMS(2)*P4-BGCMS(2)*P3
C      CALL DT_LTNUC(P3,P4,PYQREC(3),PYQREC(4),3)


C      IF (IDHKK(5).NE.22) STOP 'DT_GCFEVNTQE: Bad event format'
C      DO IDIM=1,5
C         GAMM(IDIM)=PHEP(IDIM,5)
C      ENDDO

*********transform pythia entries from lab to c.m.s of gamma+p********
C      eveBETA(1)=-(PHEP(1,2)+GAMM(1))/(PHEP(4,2)+GAMM(4))
C      eveBETA(2)=-(PHEP(2,2)+GAMM(2))/(PHEP(4,2)+GAMM(4))

**!!!!!!!!remember to change the sign of proton beam in LT*****      
C      eveBETA(3)=-(PHEP(3,2)+GAMM(3))/(PHEP(4,2)+GAMM(4))

C      eveBETA(4)=eveBETA(1)*eveBETA(1)+eveBETA(2)*eveBETA(2)+
C     & eveBETA(3)*eveBETA(3)

C      GAA=1./SQRT(1-eveBETA(4))
      
C      eveBETA(1)=eveBETA(1)*GAA
C      eveBETA(2)=eveBETA(2)*GAA
C      eveBETA(3)=eveBETA(3)*GAA

C
C      if(IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
C         WRITE(*,*)'Transformation from Lab to HCMS'
C         WRITE(*,*)'GAA,eveBETA(1-3): ',GAA,' ',eveBETA(1),' ',
C     &        eveBETA(2),' ',eveBETA(3)
C         WRITE(*,*)'Lab GAMM(1-4): ',GAMM(1),' ',GAMM(2),' ',GAMM(3),
C     &        GAMM(4)
C         WRITE(*,*)'HCMS P1,P2,P3,P4: ',P1,' ',P2,' ',P3,' ',P4
C         WRITE(*,*)'Photon angles in HCMS'
C         WRITE(*,*)'COD,SID,COF,SIF: ',COD,' ',SID,' ',COF,' ',SIF
C      endif
********get the position of particles in nucleon rest frame***         
c... we simply use the position of involved nucleon for all the
c... particles         
c... Mark - 2017-02-28 treat recoiling extra nucleons differently
C

 9999 CONTINUE

      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5) .AND. MODE.LT.2) then
         print*
         print*,'DT_GCFEVNTQE mode=',MODE
         print*,'4-momentum & charge totals for Status=1 & 1000:'
         print*,'PP1=',PPSUM(1),' PP2=',PPSUM(2),' PP3=',PPSUM(3),
     &          ' PP4=',PPSUM(4),' QCHG=',QCHG
         print*
      endif

      RETURN
      END


*=====dt_pyout=========================================================
*used for the output of pythia event list and statistics information
      SUBROUTINE DT_GCFOUTQE(MODE)     
 
*     input:
*           MODE: 1:reject statistics - not really used
*                 2:event output
*                 3:total statistics print - JUST CLOSES FILE
*                 4:event output to screen (for debugging) 
*
*     Variable OLDOUT from "beagle.inc" controls the output format.
*                        lcweight           oldweight (internal / output)
*     OLDOUT=.FALSE.      RAEVT             SIGEFF / sigma_rad
*     OLDOUT=.TRUE.         --              RAEVT
*
      IMPLICIT NONE
C      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

C      include 'pythia.inc'              ! All PYTHIA commons blocks
C      include "mc_set.inc"
C      include "py6strf.inc"
C      include "mcRadCor.inc"
C      include "radgen.inc"
C      include "phiout.inc"
      include "beagle.inc"

      EXTERNAL PYCHGE, NBARY, AZMASS, PYMASS
      INTEGER  PYCHGE, NBARY 
      DOUBLE PRECISION AZMASS, PYMASS

* event history
      INTEGER NMXHKK
      DOUBLE PRECISION MNGCF
      PARAMETER (NMXHKK=200000)
      PARAMETER (MNGCF=0.93892D0)
   
      LOGICAL ISHADR
      EXTERNAL ISHADR

      INTEGER NHKK,NEVHKK,ISTHKK,IDHKK,JMOHKK,JDAHKK
      DOUBLE PRECISION PHKK,VHKK,WHKK
      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history - IHIST renamed because it's in pythia.inc
      INTEGER IDRES,IDXRES,NOBAM,IDBAM,IDCH,NPOINT,IHISTDPM
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHISTDPM(2,NMXHKK)

C     Extensions for BEAGCF - variables that we would look up in a 
C     Pythia common blocks for MCGENE=5. Note: RECTYPE will be 
C     reported as process.
      INTEGER RECTYPE, LEADTYPE
      DOUBLE PRECISION LEPTONPHI
      COMMON /BEAGCF/ RECTYPE, LEADTYPE, LEPTONPHI

      DOUBLE PRECISION ZERO,ONE,TINY10
      INTEGER MAXNCL
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY10=1.0D-10,MAXNCL=260)

* Glauber formalism: collision properties
      DOUBLE PRECISION RPROJ,RTARG,BIMPAC
      INTEGER NWTSAM,NWASAM,NWBSAM,NWTACC,NWAACC,NWBACC
      COMMON /DTGLCP/ RPROJ,RTARG,BIMPAC,
     &                NWTSAM,NWASAM,NWBSAM,NWTACC,NWAACC,NWBACC

CC...Pythia event counter (since we keep PYINITing) Mark 2017-01-31
C      COMMON /PYCNTR/ MYNGEN
C      INTEGER MYNGEN
C      SAVE /PYCNTR/

* lorentz transformation parameter
C      COMMON /LTPARA/ BGTA(4), GAMM(5), eveBETA(4), GAA
      COMMON /ROTATE/ COF,COD,SIF,SID
      DOUBLE PRECISION COF,COD,SIF,SID

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 
      INTEGER NEVENT,ICASCA
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
      DOUBLE PRECISION PINIPR,PINITA,PRCLPR,PRCLTA,TRCLPR,TRCLTA
      COMMON /DTRNU1/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5),
     &                TRCLPR(5),TRCLTA(5),LRCLPR,LRCLTA

* treatment of residual nuclei: properties of residual nuclei
      DOUBLE PRECISION AMRCL0,EEXC,EEXCFI
      INTEGER NTOT,NPRO,NN,NH,NHPOS,NQ,NTOTFI,NPROFI
      COMMON /DTRNU2/ AMRCL0(2),EEXC(2),EEXCFI(2),
     &                NTOT(2),NPRO(2),NN(2),NH(2),NHPOS(2),NQ(2),
     &                NTOTFI(2),NPROFI(2)

* treatment of residual nuclei: wounded nucleons
      INTEGER NPW,NPW0,NPCW,NTW,NTW0,NTCW,IPW,ITW
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

      DOUBLE PRECISION P1, P2, P3, P4, PTOT, PP1, PP2, PP3, PP4, PP5

      integer NEV, NPRT, ievent, genevent, I, tracknr
      integer lastgenevent, idum1, idum2, initseed, nrtrack
      DOUBLE PRECISION radgamE, radgamp, radgamEnucl
      DOUBLE PRECISION epznucl
      CHARACTER PARAM*100
      LOGICAL UseLut, GenLut
      INTEGER MODE

C     Local
      DOUBLE PRECISION BETA2,P2TEMP,XCALC
      DOUBLE PRECISION P5SUM(5)
      LOGICAL VERBOSE
      DOUBLE PRECISION GAMTMP1,GAMTMP2,GAMTMP3
      DOUBLE PRECISION BGAMTMP1,BGAMTMP2,BGAMTMP3
      DOUBLE PRECISION ESUM1,ESUM2,ESUM3,PZSUM1,PZSUM2,PZSUM3 
      DOUBLE PRECISION AMASS1,AMASS2,AMASS3,EINA1,EINA2,EINA3
      DOUBLE PRECISION EINE,EIN1,EIN2,EIN3,PZIN
      INTEGER IZSUM,IASUM,ILEPT,IBOSON,J,IDIM,NDIM,ILEPIN,INUCIN
      parameter ( NDIM=4 )
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./

C      if (mcRadCor_EBrems.gt.0.) then
C         radgamEnucl=sqrt(dplabg(1)**2+dplabg(2)**2+dplabg(3)**2)
C         radgamE=pgamma*radgamEnucl-pgamma*pbeta*dplabg(3)
C         radgamp=-pgamma*pbeta*radgamEnucl+pgamma*dplabg(3)
CC         write(*,*) radgamEnucl, radgamE, dplabg(3), radgamp
C      else
C        radgamEnucl=0D0
C        radgamE=0D0
C        radgamp=0D0 
C      endif

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
         IF(J.GT.NPOINT(1)) THEN
            IF (INRLEV.GT.2) THEN
C We did things in a different order...
C     Boost from HCMS-G to IRF-G
               P1=PHKK(1,J)
               P2=PHKK(2,J)
               P3=PHKK(3,J)
               P4=PHKK(4,J)
C               CALL DT_LTNUC(PHKK(3,J),PHKK(4,J),P3,P4,-3)
C     Rotate from gamma* +z back to IRF-B=LAB
C********rotate back from the gamma* +z direction***********
               PP1=COF*(COD*P1+SID*P3)-SIF*P2
               PP2=SIF*(COD*P1+SID*P3)+COF*P2
               PP3=COD*P3-SID*P1
               PP4=P4
C               P1=COF*(COD*PHKK(1,J)+SID*PHKK(3,J))-SIF*PHKK(2,J)
C               P2=SIF*(COD*PHKK(1,J)+SID*PHKK(3,J))+COF*PHKK(2,J)
C               P3=COD*PHKK(3,J)-SID*PHKK(1,J)
C               P4=PHKK(4,J)
CC**************transform back to the lab frame**************      
C               call DT_DALTRA(GAA,-eveBETA(1),-eveBETA(2),-eveBETA(3),
C     &              P1,P2,P3,P4,PTOT,PP1,PP2,PP3,PP4)
CC     WRITE(89,996) J,ISTHKK(J),IDHKK(J),JMOHKK(1,J),
CC     &          JMOHKK(2,J),JDAHKK(1,J),JDAHKK(2,J),
CC     &          PP1,PP2,-PP3,PP4, !remember to change the sign of z back
CC     &              PHKK(5,J)
CC  996      FORMAT(I5,I5,I8,4I5,5F17.5)
C MDB Since the above matrix has det +1, make this a rotation too.
C MDB Rotation not needed for INRLEV=3. Not sure about 1 or 2.
C               PHKK(1,J)=-PP1
C               PHKK(2,J)=PP2      
C               PHKK(3,J)=-PP3
C               PHKK(4,J)=PP4
               PHKK(1,J)=PP1
               PHKK(2,J)=PP2      
               PHKK(3,J)=PP3
               PHKK(4,J)=PP4
CC     PHKK(3,J)=pgamma*(P3+pbeta*P4)
CC     PHKK(4,J)=pgamma*(P4+pbeta*P3)
            ENDIF
c...find the exchanged boson and out e- to make it fit root tree making rules
c...in the following steps
            IF ((ISTHKK(J).EQ.3).AND.(IDHKK(J).EQ.22.OR.IDHKK(J).EQ.23)
     &           .AND. ((JMOHKK(1,J).EQ.(NPOINT(1)+1)).OR.INRLEV.GT.2) 
     &           ) THEN
               IBOSON=J
            ELSEIF (((JMOHKK(1,J).EQ.4).OR.
     &              (ISTHKK(J).EQ.99.AND.INRLEV.GT.2)).AND.
     &              (ABS(IDHKK(J)).EQ.11.OR.ABS(IDHKK(J)).EQ.13)) THEN
               ISTHKK(J)=1
               ILEPT=J
c...2017-01-02 MDB Fill some new event variables
            ELSE
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

C      genevent=MYNGEN-lastgenevent
      tracknr = NHKK
C      if (mcRadCor_EBrems.gt.0.) then
C         nrtrack=tracknr+1
C      else
         nrtrack=tracknr
C      endif

c...print a title for the event file
      If (FIRST) then
C        write(29,*)' PYTHIA EVENT FILE '
        write(29,*)' BEAGLE EVENT FILE '
        write(29,*)'============================================'
        write(29,31) 
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

C      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
C         XCALC = VINT(307)/VINT(309)/(VINT(302)-VINT(4)**2-VINT(303)**2)
C         write(*,*) "MN VINT(4): ",VINT(4),"massp: ",massp
C         write(*,*) "MN VINT(304): ",VINT(304),"MN2 VINT(308)",VINT(308)
C         write(*,*) "X(VINT-calc): ",XCALC," XBJOUT: ",XBJOUT
C         write(*,*) "beamx VINT(305): ",VINT(305),
C     &        "targx VINT(306): ",VINT(306)
C         write(*,*) "Y VINT(309): ",VINT(309)," YYOUT:  ",YYOUT
C         write(*,*) "W2 VINT(2): ",VINT(2)," W2OUT: ",W2OUT
C         write(*,*) "Nu(from X): ",VINT(307)/(2.*VINT(4)*XCALC),
C     &        " NUOUT: ",NUOUT
C         write(*,*) "Q2 VINT(307): ",VINT(307)," Q2OUT: ",Q2OUT
C         write(*,*) "Recalc. W2::",VINT(4)**2+VINT(307)*(1./XCALC-1.)
C         write(*,*) "Recalc. nu from W2:",
C     &              (VINT(2)+VINT(307)-VINT(4)**2)/(2.*VINT(4))
C         write(*,*) ""
C      endif

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
C      IF (tracknr.LT.IT+2) STOP "DT_GCFOUTQE FATAL: Too few tracks."
C      DO I=2,IT+1
      
      IF (INRLEV.LE.2) THEN
         DO I=1,3
            IF(IDHKK(I).NE.2112 .AND. IDHKK(I).NE.2212 .AND. 
     &           IDHKK(I).NE.80000) THEN
               WRITE(*,*) 'DT_GCFOUTQE ERROR: IN-NUCLEON IDHKK=',
     &              IDHKK(I)
            ELSEIF (ISTHKK(I).EQ.18 .OR. ISTHKK(I).EQ.12) THEN
               NWND = NWND + 1
               IF (IDHKK(I).EQ.2212) NWDCH = NWDCH + 1
            ELSEIF (ISTHKK(I).NE.14) THEN
               WRITE(*,*) 'DT_GCEOUTQE ERROR: IN-NUCLEON ISTHKK=',
     &              ISTHKK(I)
            ENDIF
         ENDDO
      ELSE
         DO I=2,IT+1
            IF (ISTHKK(I).EQ.18 .OR. ISTHKK(I).EQ.12) THEN
               NWND = NWND + 1
               IF (IDHKK(I).EQ.2212) NWDCH = NWDCH + 1
            ENDIF
         ENDDO
      ENDIF
      IF (USERSET.GE.8) THEN         
         VERBOSE = (IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5))
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
         IF (INRLEV.GT.1) THEN
            CALL EVTSUM(P5SUM,IZSUM,IASUM,VERBOSE,.FALSE.)
         ELSE
C     For pass-through mode, sum with 1000 treated as stable 
            DO IDIM=1,NDIM
               P5SUM(IDIM)=0.0D0
            ENDDO
            IZSUM=0
            IASUM=0
            DO I=1,NHKK
               IF (ISTHKK(I).EQ.1 .OR. ISTHKK(I).EQ.1000) THEN
                  IZSUM=IZSUM+IDXRES(I)
                  IASUM=IASUM+IDRES(I)
                  DO IDIM=1,NDIM
                     P5SUM(IDIM) = P5SUM(IDIM)+PHKK(IDIM,I)
                  ENDDO
               ENDIF
               P5SUM(5) = SQRT(P5SUM(4)*P5SUM(4)-P5SUM(1)*P5SUM(1)
     &              -P5SUM(2)*P5SUM(2)-P5SUM(1)*P5SUM(1))
            ENDDO
         ENDIF
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
     &  WRITE(*,*)'GCFOUTEP ERROR: NCOLLT .NE. NTW0: ',NCOLLT,' ',NTW0
C      WRITE(*,*)
C      WRITE(*,*) 'Output XBJ= ',XBJOUT
C      WRITE(*,*) 'Output RAVAL= ',RAVAL
C      WRITE(*,*) 'Output RAEVT= ',RAEVT
C      WRITE(*,*)
C
C      Not sure how "xbeamparton = pari(33)" is defined. So 0.0.
C      Leave 0 the "hat" variables and the "radcorr" variables
      if (OLDOUT) SIGEFF=0.0D0
      write(29,33) 0, ievent, genevent, ltype, it, itz, pzlep, 
     &        pztarg, pznucl, crang, crori, RECTYPE, LEADTYPE,
     &        LEADTYPE, 1.0D0, 22, 0.0D0, 0.0D0,
     &        YYOUT, Q2OUT, XBJOUT, W2OUT, NUOUT,
     &        LEPTONPHI, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 
     &        0.0D0, 0.0D0, SIGEFF, 0.0D0, 0.0D0,
     &        0.0D0, BBEA, PHIB, THKB, THKSCL, NCOLLT, NCOLLI,
     &        NWND, NWDCH,
     &        NNEVAP, NPEVAP, AREMN, NINC, NINCCH, DFIRST, DAVG,
     &        PXF, PYF, PZF, EEXC(2), RAEVT, USER1, USER2, USER3, 
     &        nrtrack+4 
C Note: Use E rather than F format for GCF RAEVT (weight)
 33   format((I4,1x,$),(I10,1x,$),4(I4,1x,$),4(f12.6,1x,$),3(I4,1x,$),
     &     I6,1x,$,f9.6,1x,$,I6,1x,$,2(f12.6,1x,$),7(f18.11,3x,$),
     &     7(f19.9,3x,$),(e17.8,1x,$),3(f19.9,3x,$),4(f10.6,1x,$),
     &     9(I5,1x,$),2(f10.6,1x,$),
     &     3(f15.6,1x,$),f12.6,1x,$,4(e17.8,1x,$),I6,/)
      write(29,*)'============================================'

***************standard output for particle info************************
c...add 2 beam information at first to fit into root tree making rule      
c... MDB Change these lines to use the correct pid for lepton+nucleon!
C      I=NPOINT(1)+1   
      IF (INRLEV.GT.2) THEN
         ILEPIN = NPOINT(1)+4
         INUCIN = NPOINT(1)+1
      ELSE
         ILEPIN = NPOINT(1)+1
         INUCIN = 1
      ENDIF
C     Incoming lepton
      write(29,34) 1,21,IDHKK(ILEPIN),0,0,0,0,
     &     PHKK(1,ILEPIN),PHKK(2,ILEPIN),PHKK(3,ILEPIN),PHKK(4,ILEPIN),
     &     PHKK(5,ILEPIN),VHKK(1,ILEPIN),VHKK(2,ILEPIN),VHKK(3,ILEPIN)
     &     ,0,0,0
c...  Incoming nucleon for eicsmear. Use MNGCF for consistent kinematics.
      write(29,34) 2,21,IDHKK(INUCIN),0,0,0,0,0.0D0,0.0D0,0.0D0,
     &     MNGCF,MNGCF,
     &     VHKK(1,INUCIN),VHKK(2,INUCIN),VHKK(3,INUCIN),1,IDXRES(1),0
c...add the scattered lepton 
      write(29,34) 3,21,IDHKK(ILEPT),0,1,ILEPT+4,0,
     &     PHKK(1,ILEPT),PHKK(2,ILEPT),PHKK(3,ILEPT),
     &     PHKK(4,ILEPT),PHKK(5,ILEPT),VHKK(1,ILEPT),
     &     VHKK(2,ILEPT),VHKK(3,ILEPT)
     &     ,0,0,0
c...add the exchanged boson  
      write(29,34) 4,21,IDHKK(IBOSON),0,1,IBOSON+4,0,
     &     PHKK(1,IBOSON),PHKK(2,IBOSON),PHKK(3,IBOSON),
     &     PHKK(4,IBOSON),PHKK(5,IBOSON),VHKK(1,IBOSON),
     &     VHKK(2,IBOSON),VHKK(3,IBOSON)
     &     ,0,0,0
 
      DO I=1,tracknr
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
         IF (IDHKK(I).EQ.80000) THEN
            ZOUT = IDXRES(I)
            AOUT = IDRES(I)
         ELSE
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
!!!dump nuclear remnants into final state particles
         IF (ISTHKK(I).EQ.-1) THEN
            KSOUT = 1
            BAMOUT = 3
         ELSEIF (ISTHKK(I).EQ.1001) THEN
            KSOUT = 1
            BAMOUT = 4
         ENDIF
         write(29,34) I+4, KSOUT, IDHKK(I), JM2OUT, JM1OUT, 
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

C         lastgenevent=MYNGEN

c         print*,'output finished'
      RETURN

c...mode 3 is used to print the whole statistics information
3     CONTINUE
      close(29)

C...Check pdf status       
C      call PDFSTA
      
      RETURN

c...mode 4 is used to output the event list to screen in current frame
c...without a lot of reformatting and rearranging (Mark 08/17/2016)
c...using the old event header (Mark 01/02/2017)
4     CONTINUE

C      genevent=MYNGEN-lastgenevent
      tracknr = NHKK
C      if (mcRadCor_EBrems.gt.0.) then
C         nrtrack=tracknr+1
C      else
         nrtrack=tracknr
C      endif

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
     &     pztarg, pznucl, crang, crori, RECTYPE, LEADTYPE,
     &     LEADTYPE, 1.0D0, 22, 0.0D0, 0.0D0,
     &     YYOUT, Q2OUT, XBJOUT, W2OUT, NUOUT,
     &     LEPTONPHI, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 
     &     0.0D0, 0.0D0, SIGEFF, 0.0D0, 0.0D0,
     &     0.0D0, BBEA, PHIB, THKB, THKSCL, NCOLLT, NCOLLI,
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
*$ CREATE DT_FLUKAIT.FOR
*COPY DT_FLUKAIT
*
*===ficonf=============================================================*
*
      SUBROUTINE DT_FLUKAIT(IT,ITZ,IREJ)

************************************************************************
* This is the subset of the DT_FICONF routine (from S. Roesler) that   *
* interfaces to FLUKA in order handle the evaporation, fission and     *
* Fermi-break-up of residual nucleus.                                  *
*                                                                      *
* In DT_FICONF, the routine first assembles the residual nucleus.      *
* Here we assume that has been done already (e.g. by GCF).             *
*                                                                      *
* Last change 03 May 2019 by M.D. Baker                                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      PARAMETER ( LOUT = 6 )

      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY3=1.0D-3,TINY10=1.0D-10)
      PARAMETER (ANGLGB=5.0D-16)

      EXTERNAL PYCHGE, NBARY
      INTEGER PYCHGE

* event history

      PARAMETER (NMXHKK=200000)

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHIST(2,NMXHKK)

* rejection counter
      COMMON /DTREJC/ IRPT,IRHHA,IRRES(2),LOMRES,LOBRES,
     &                IRCHKI(2),IRFRAG,IRCRON(3),IREVT,
     &                IREXCI(3),IRDIFF(2),IRINC

* central particle production, impact parameter biasing
      COMMON /DTIMPA/ BIMIN,BIMAX,XSFRAC,ICENTR

* particle properties (BAMJET index convention)
      CHARACTER*8  ANAME
      COMMON /DTPART/ ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

* treatment of residual nuclei: 4-momenta
      LOGICAL LRCLPR,LRCLTA
      COMMON /DTRNU1/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5),
     &                TRCLPR(5),TRCLTA(5),LRCLPR,LRCLTA

* treatment of residual nuclei: properties of residual nuclei
      COMMON /DTRNU2/ AMRCL0(2),EEXC(2),EEXCFI(2),
     &                NTOT(2),NPRO(2),NN(2),NH(2),NHPOS(2),NQ(2),
     &                NTOTFI(2),NPROFI(2)

* statistics: residual nuclei
      COMMON /DTSTA2/ EXCDPM(4),EXCEVA(2),
     &                NINCGE,NINCCO(2,3),NINCHR(2,2),NINCWO(2),
     &                NINCST(2,4),NINCEV(2),
     &                NRESTO(2),NRESPR(2),NRESNU(2),NRESBA(2),
     &                NRESPB(2),NRESCH(2),NRESEV(4),
     &                NEVA(2,6),NEVAGA(2),NEVAHT(2),NEVAHY(2,2,240),
     &                NEVAFI(2,2)

      INCLUDE 'beagle.inc'
      INCLUDE '(DIMPAR)'
      INCLUDE '(GENSTK)'
      INCLUDE '(RESNUC)'
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( COUGFM = ELCCGS * ELCCGS / ELCMKS * 1.D-07 * 1.D+13
     &                   * 1.D-09 )
      PARAMETER ( HLFHLF = 0.5D+00 )
      PARAMETER ( FERTHO = 14.33       D-09 )
      PARAMETER ( BEXC12 = FERTHO * 72.40715579499394D+00 )
      PARAMETER ( AMUNMU = HLFHLF * AMELCT - BEXC12 / 12.D+00 )
      PARAMETER ( AMUC12 = AMUGEV - AMUNMU )
      INCLUDE '(NUCDAT)'
      INCLUDE '(PAREVT)'
      INCLUDE '(FHEAVY)'

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA

      DIMENSION INUC(2),IDXPAR(2),IDPAR(2),AIF(2),AIZF(2),AMRCL(2),
     &          PRCL(2,4),MO1(2),MO2(2),VRCL(2,4),WRCL(2,4),
     &          P1IN(4),P2IN(4),P1OUT(4),P2OUT(4)

      DIMENSION EXPNUC(2),EXC(2,260),NEXC(2,260)
      DATA EXC,NEXC /520*ZERO,520*0/
      DATA EXPNUC /4.0D-3,4.0D-3/

C Local variable for High E* protection
      DOUBLE PRECISION TMPMAS,TMPENE
      DOUBLE PRECISION PSUMTOT(4),PSUMRES(4)
      INTEGER ZSUMTOT, ZSUMRES

      IREJ   = 0
C      LRCLPR = .FALSE.
C      LRCLTA = .FALSE.

* skip residual nucleus treatment if not requested or in case
* of central collisions
C      IF ((.NOT.LEVPRT).OR.(ICENTR.GT.0).OR.(ICENTR.EQ.-1)) RETURN

C      DO 1 K=1,2
C         IDPAR(K) = 0
C         IDXPAR(K)= 0
C         NTOT(K)  = 0
C         NTOTFI(K)= 0
C         NPRO(K)  = 0
C         NPROFI(K)= 0
C         NN(K)    = 0
C         NH(K)    = 0
C         NHPOS(K) = 0
C         NQ(K)    = 0
C         EEXC(K)  = ZERO
C         MO1(K)   = 0
C         MO2(K)   = 0
C         DO 2 I=1,4
C            VRCL(K,I) = ZERO
C            WRCL(K,I) = ZERO
C    2    CONTINUE
C    1 CONTINUE
C      NFSP = 0
C      INUC(1) = IP
C      INUC(2) = IT

C Shortcut version for GCF
      KF=2
      I=KF
      IEXREM=9
      IF (ISTHKK(IEXREM).NE.1000) STOP 
     &     'DT_FLUKAIT FATAL ERROR: Excited A-2 remnant not found.' 
      DO 4 K=1,4
         VRCL(KF,K) = VHKK(K,IEXREM)
         WRCL(KF,K) = WHKK(K,IEXREM)
         PRCL(KF,K) = PHKK(K,IEXREM)
 4    CONTINUE
      NQ(KF)=IDXRES(IEXREM)
      AIF(KF)=DBLE(IDRES(IEXREM))
      AIZF(KF)=DBLE(IDXRES(IEXREM))
      INUC(KF)=IT
      NTOT(KF)=IDRES(IEXREM)
      NPRO(KF)=IDXRES(IEXREM)
      NH(KF)=0
      NHPOS(KF)=0
      NN(KF)=NTOT(KF)-NPRO(KF)
      NHPOS(KF)=0
      AMRCL0(KF) = PHKK(5,IEXREM)-EEXC(2)
      AMRCL(KF)=PHKK(5,IEXREM)
      PTORCL = SQRT(PRCL(KF,1)**2+PRCL(KF,2)**2+PRCL(KF,3)**2)
      IF (NEVENT.LE.5) THEN 
         WRITE(*,*) 'A, Z, Mass from GCF, Mass from AZMASS 0,1,3'
         WRITE(*,*) A, Z, PHKK(5,IEXREM)-EEXC(2), 
     &        AZMASS(IDRES(IEXREM),IDXRES(IEXREM),0),
     &        AZMASS(IDRES(IEXREM),IDXRES(IEXREM),1),
     &        AZMASS(IDRES(IEXREM),IDXRES(IEXREM),3)
         WRITE(*,*) 'EEXC, AMRCL, AMRCL0, PRCL(1:4),M_PRCL,PTORCL'
         WRITE(*,*) EEXC(KF), AMRCL(KF), AMRCL0(KF), PRCL(KF,1), 
     &        PRCL(KF,2), PRCL(KF,3), PRCL(KF,4), 
     &        SQRT(PRCL(KF,4)**2-PTORCL**2),PTORCL
      ENDIF

      ICOR   = 0
      INORCL = 0

*
C               LLCPOT = .TRUE.
C               IF ( LLCPOT ) THEN
C                  NNCHIT = MAX ( INUC (I) - NTOT (I), 0 )
C                  DLKPRH = ZERO
C                  RDCORE = 1.14D+00 * DBLE(INUC(I))**(ONE/3.D+00)
C*  Take out roughly one/half of the skin:
C                  RDCORE = RDCORE - 0.5D+00
C                  FRCFLL = RDCORE**3
C                  PRSKIN = (RDCORE+2.4D+00)**3 - FRCFLL
C                  PRSKIN = 0.5D+00 * PRSKIN / ( PRSKIN + FRCFLL )
C                  FRCFLL = ONE - PRSKIN
C                  FRMRDC = FRCFLL + 0.5D+00 * PRSKIN
C                  REDORI = ONE / ( FRMRDC )**(2.D+00/3.D+00)
C                  IF ( NNCHIT .GT. 0 ) THEN
C                     REDCTN = ZERO
C                     DO 1230 NCH = 1, NNCHIT
C                        IF (DT_RNDM(PRFRMI) .LT. PRSKIN) THEN
C                           PRFRMI = (( ONE - 2.D+00 * DLKPRH )
C     &                            * DT_RNDM(PRFRMI))**0.333333333333D+00
C                        ELSE
C                           PRFRMI = ( ONE - 2.D+00 * DLKPRH
C     &                            * DT_RNDM(PRFRMI))**0.333333333333D+00
C                        END IF
C                        REDCTN = REDCTN + PRFRMI**2
C 1230                CONTINUE
C                     REDCTN = REDCTN / DBLE (NNCHIT)
C                  ELSE
C                     REDCTN = 0.5D+00
C                  END IF
C                  EEXC  (I) = EEXC   (I) * REDCTN / REDORI
C                  AMRCL (I) = AMRCL0 (I) + EEXC (I)
C                  PRCL(I,4) = SQRT ( PTORCL**2 + AMRCL(I)**2 )
C               END IF
C**
C               IF (ICASCA.EQ.0) THEN
C                  EXPNUC(I) = EEXC(I)/MAX(1,INUC(I)-NTOT(I))
C                  M = MIN(NTOT(I),260)
C                  EXC(I,M)  = EXC(I,M)+EEXC(I)
C                  NEXC(I,M) = NEXC(I,M)+1
C               ENDIF
C            ENDIF
C         ELSEIF (NTOT(I).EQ.1) THEN
C            WRITE(LOUT,1003) I
C 1003       FORMAT(1X,'FICONF:   warning! NTOT(I)=1? (I=',I3,')')
C            GOTO 9999
C         ELSE
C            AMRCL0(I) = ZERO
C            AMRCL(I)  = ZERO
C            EEXC(I)   = ZERO
C            INORCL    = INORCL+I
C         ENDIF
C    7 CONTINUE

C      PRCLPR(5) = AMRCL(1)
C      PRCLTA(5) = AMRCL(2)

C      IF (ICOR.GT.0) THEN
C         IF (INORCL.EQ.0) THEN
C* one or both residual nuclei consist of one nucleon only, transform
C* this nucleon on mass shell
C            DO 9 K=1,4
C               P1IN(K) = PRCL(1,K)
C               P2IN(K) = PRCL(2,K)
C    9       CONTINUE
C            XM1 = AMRCL(1)
C            XM2 = AMRCL(2)
C            CALL DT_MASHEL(P1IN,P2IN,XM1,XM2,P1OUT,P2OUT,IREJ1)
C            IF (IREJ1.GT.0) THEN
C               WRITE(LOUT,*) 'ficonf-mashel rejection'
C               GOTO 9999
C            ENDIF
C            DO 10 K=1,4
C               PRCL(1,K) = P1OUT(K)
C               PRCL(2,K) = P2OUT(K)
C               PRCLPR(K) = P1OUT(K)
C               PRCLTA(K) = P2OUT(K)
C   10       CONTINUE
C            PRCLPR(5) = AMRCL(1)
C            PRCLTA(5) = AMRCL(2)
C         ELSE
C            IF (IOULEV(3).GT.0)
C     &      WRITE(LOUT,1001) NEVHKK,INT(AIF(1)),INT(AIZF(1)),
C     &                       INT(AIF(2)),INT(AIZF(2)),AMRCL0(1),
C     &                       AMRCL(1),AMRCL(1)-AMRCL0(1),AMRCL0(2),
C     &                       AMRCL(2),AMRCL(2)-AMRCL0(2)
C 1001       FORMAT(1X,'FICONF:   warning! no residual nucleus for',
C     &             ' correction',/,11X,'at event',I8,
C     &             ',  nucleon config. 1:',2I4,' 2:',2I4,
C     &             2(/,11X,3E12.3))
C            IF (NLOOP.LE.500) THEN
C               GOTO 9998
C            ELSE
C               IREXCI(1) = IREXCI(1)+1
C            ENDIF
C         ENDIF
C      ENDIF

* update counter
C     IF (NRESEV(1).NE.NEVHKK) THEN
C        NRESEV(1) = NEVHKK
C        NRESEV(2) = NRESEV(2)+1
C     ENDIF
C      NRESEV(2) = NRESEV(2)+1
C      DO 15 I=1,2
C         EXCDPM(I)   = EXCDPM(I)+EEXC(I)
C         EXCDPM(I+2) = EXCDPM(I+2)+(EEXC(I)/MAX(NTOT(I),1))
C         NRESTO(I) = NRESTO(I)+NTOT(I)
C         NRESPR(I) = NRESPR(I)+NPRO(I)
C         NRESNU(I) = NRESNU(I)+NN(I)
C         NRESBA(I) = NRESBA(I)+NH(I)
C         NRESPB(I) = NRESPB(I)+NHPOS(I)
C         NRESCH(I) = NRESCH(I)+NQ(I)
C   15 CONTINUE

* evaporation
      IF (LEVPRT) THEN
C         DO 13 I=1,2      ! NOTE: I=2 was set above.
* initialize evaporation counter
            EEXCFI(I) = ZERO
            IF ((INUC(I).GT.1).AND.(AIF(I).GT.ONE).AND.
     &          (EEXC(I).GT.ZERO)) THEN
* put residual nuclei into DTEVT1 - NOT NEEDED
C               IDRCL = 80000
C               JMASS = INT( AIF(I))
C               JCHAR = INT(AIZF(I))
*  the following patch is required to transmit the correct excitation
*   energy to Eventd
               IF (ITRSPT.EQ.1) THEN
                  STOP 'DT_FLUKAIT cannot run with ITRSPT=1'
C               CALL DT_EVTPUT(1000,IDRCL,MO1(I),MO2(I),PRCL(I,1),
C     &              PRCL(I,2),PRCL(I,3),PRCL(I,4),JMASS,JCHAR,0)
**sr 22.6.97
C               NOBAM(NHKK) = I
**
C               DO 14 J=1,4
C                  VHKK(J,NHKK) = VRCL(I,J)
C                  WHKK(J,NHKK) = WRCL(I,J)
C   14          CONTINUE
*  interface to evaporation module - fill final residual nucleus into
*  common FKRESN
*   fill resnuc only if code is not used as event generator in Fluka
               ELSE 
                  PXRES  = PRCL(I,1)
                  PYRES  = PRCL(I,2)
                  PZRES  = PRCL(I,3)
                  IBRES  = NPRO(I)+NN(I)+NH(I)
                  ICRES  = NQ(I)
                  ANOW   = DBLE(IBRES)
                  ZNOW   = DBLE(ICRES)
                  PTRES  = SQRT(PXRES**2+PYRES**2+PZRES**2)
*   ground state mass of the residual nucleus (should be equal to AM0T)

                  AMNRES = AMRCL0(I)
                  AMMRES = AMNAMA ( AMNRES, IBRES, ICRES )
                  IF (NEVENT.LE.5) THEN
                     WRITE(*,*)'AMNRES,AMMRES: ',AMNRES,AMMRES
                  ENDIF
*  common FKFINU
                  TV = ZERO
*   kinetic energy of residual nucleus
                  TVRECL = PRCL(I,4)-AMRCL(I)
*   excitation energy of residual nucleus
                  TVCMS  = EEXC(I)
*   Disallow very large E* which leads to infinite loops - MDB
                  IF (EEXCMAX.GT.TINY3 .AND. EEXC(I).GT.EEXCMAX) THEN
                     WRITE(*,*) 'FICONF: ERROR! EEXC>max',EEXC(I),
     &                    'in event',NEVHKK,'.'
                     WRITE(*,*) '  Setting value to',EEXCMAX,' for',
     &                    ' EVEVAP(). Energy not conserved!'
                     WRITE(*,*)
                     TVCMS=EEXCMAX
                     TVRECL=PRCL(I,4)-AMNRES-TVCMS
                  ENDIF

                  PTOLD  = PTRES
                  PTRES  = SQRT(ABS(TVRECL*(TVRECL+
     &                          2.0D0*(AMMRES+TVCMS))))
                  IF (PTOLD.LT.ANGLGB) THEN
                     CALL DT_RACO(PXRES,PYRES,PZRES)
                     PTOLD = ONE
                  ENDIF
                  PXRES = PXRES*PTRES/PTOLD
                  PYRES = PYRES*PTRES/PTOLD
                  PZRES = PZRES*PTRES/PTOLD
                  IF (NEVENT.LE.5) THEN
                     WRITE(*,*) 
     &       'DPMJET PRCL(1-4), AMRCL0, AMRCL, PTOLD, TVRECL, TVCMS'
                     WRITE(*,*) PRCL(I,1),PRCL(I,2),PRCL(I,3),PRCL(I,4),
     &                    AMRCL0(I),AMRCL(I),PTOLD,TVRECL,TVCMS
                     WRITE(*,*) 'FLUKA P(X,Y,Z)RES, PTRES'
                     WRITE(*,*) PXRES,PYRES,PZRES,PTRES
                  ENDIF
* zero counter of secondaries from evaporation
                  NP = 0
* evaporation
                  WE = ONE

                  NPHEAV = 0
                  LRNFSS = .FALSE.
                  LFRAGM = .FALSE.
                  CALL EVEVAP(WE)

* put evaporated particles and residual nuclei to DTEVT1
                  MO = NHKK
                  CALL DT_EVA2HE(MO,EXCITF,I,IREJ1)
               ENDIF
               EEXCFI(I) = EXCITF
               EXCEVA(I) = EXCEVA(I)+EXCITF
            ENDIF
C   13    CONTINUE
      ENDIF

C     Assign, charge and baryon # to new particles.
      DO ILST=10,NHKK
         IF (IDHKK(ILST).NE.80000) THEN
            IDXRES(ILST)=PYCHGE(IDHKK(ILST))/3
            IDRES(ILST)=NBARY(IDHKK(I))/3
         ENDIF
      ENDDO

      IF (NEVENT.LE.5) THEN
         ZSUMTOT=0
         ZSUMRES=0
         DO IDIM=1,4
            PSUMTOT(IDIM)=0.0D0
            PSUMRES(IDIM)=0.0D0
         ENDDO
         WRITE(*,*) 'Post-evaporation E*:',EXCITF
         WRITE(*,*) 'Post-evaporation residual event'
         WRITE(*,*) 
         DO ILST=1,NHKK
            IF (ISTHKK(ILST).EQ.1) THEN
               ZSUMTOT=ZSUMTOT+IDXRES(ILST)
               DO IDIM=1,4
                  PSUMTOT(IDIM)=PSUMTOT(IDIM)+PHKK(IDIM,ILST)
               ENDDO
            ENDIF
            IF (ILST.GE.9) THEN
               WRITE(*,*) ILST,ISTHKK(ILST),IDHKK(ILST),PHKK(1,ILST),
     &           PHKK(2,ILST),PHKK(3,ILST),PHKK(4,ILST),PHKK(5,ILST),
     &           IDRES(ILST),IDXRES(ILST)
               IF (ISTHKK(ILST).EQ.1) THEN
                  ZSUMRES = ZSUMRES + IDXRES(ILST)
                  DO IDIM=1,4
                     PSUMRES(IDIM)=PSUMRES(IDIM)+PHKK(IDIM,ILST)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         WRITE(*,*) 'Q,p4 TOTAL Sum:',
     &        ZSUMTOT,PSUMTOT(1),PSUMTOT(2),PSUMTOT(3),PSUMTOT(4)
         WRITE(*,*) 'Q,p4 Res.  Sum:',
     &        ZSUMRES,PSUMRES(1),PSUMRES(2),PSUMRES(3),PSUMRES(4)
      ENDIF

      RETURN

CC9998 IREXCI(1) = IREXCI(1)+1
C 9998 IREJ   = IREJ+1
C 9999 CONTINUE
C      LRCLPR = .TRUE.
C      LRCLTA = .TRUE.
C      IREJ   = IREJ+1
C      RETURN
      END
*=====dt_gcfkinemqe========================================================
* Get e' and gamma* from GCF. Note: We have rotated the GCF event so that
* e is along -z while dpmjet wants it along +z (ironically how it was in
* the GCF root file to start with.
*
      SUBROUTINE DT_GCFKINEMQE(PPL1,PLTOT,PPG,PGTOT)
*
*     MDB 2020-03-04 Latest Version
* 
*     output: 
*           PPL1 = scattered lepton 4-momentum
*           PLTOT = 3-momentum magnitude for PPL1
*           PPG = gamma* 4-momentum
*           PGTOT = 3-momentum magnitude for PPL1

      IMPLICIT NONE
C      include "beagle.inc"

C     Calling sequence:
      DOUBLE PRECISION PPL1(4),PPG(4)
      DOUBLE PRECISION PLTOT,PGTOT


C     Temporary storage for GCF event until DTEVT1 is ready for it.
      INTEGER MAXGCF,NGCF,IOFF
      PARAMETER (MAXGCF=100)
      INTEGER ISTGCF, IDGCF, JMOGCF, JDAGCF, ARESGCF, ZRESGCF, NOBAMGCF
      INTEGER IDBAMGCF
      DOUBLE PRECISION PGCF
      COMMON /GCFDAT/ ISTGCF(MAXGCF),IDGCF(MAXGCF),JMOGCF(2,MAXGCF),
     &     JDAGCF(2,MAXGCF), PGCF(5,MAXGCF),ARESGCF(MAXGCF),
     &     ZRESGCF(MAXGCF),NOBAMGCF(MAXGCF),IDBAMGCF(MAXGCF)

      INTEGER ISCAT, IGAMM
      PARAMETER (ISCAT=6)
      PARAMETER (IGAMM=5)
C  Local

      INTEGER IDIM
      DOUBLE PRECISION PSIGN

C     Rotate event by 180 degrees about the y-axis
      DO IDIM=1,4
         PPL1(IDIM) = PGCF(IDIM,ISCAT)
         PPG(IDIM) = PGCF(IDIM,IGAMM)
      ENDDO
      PLTOT = SQRT(PPL1(1)*PPL1(1)+PPL1(2)*PPL1(2)+PPL1(3)*PPL1(3))
      PGTOT = SQRT(PPG(1)*PPG(1)+PPG(2)*PPG(2)+PPG(3)*PPG(3))

      RETURN
      END
* Stripped down version of DT_GLAUBE for GCF.
*
*===glaube=============================================================*
*
      SUBROUTINE DT_GCFGLAUBE(NA,NB,IJPROJ,B,INTT,INTA,INTB,JS,JT,NIDX)

************************************************************************
* Calculation of configuartion of interacting nucleons for one event.  *
*    NB / NB    mass numbers of proj./target nuclei           (input)  *
*    B          impact parameter                              (output) *
*    INTT       total number of wounded nucleons                 "     *
*    INTA / INTB number of wounded nucleons in proj. / target    "     *
*    JS / JT(i) number of collisions proj. / target nucleon i is       *
*                                                   involved  (output) *
*    NIDX       index of projectile/target material            (input) *
*               = -2 call within FLUKA transport calculation           *
* This is an update of the original routine SHMAKO by J.Ranft/HJM      *
* This version dated 22.03.96 is revised by S. Roesler                 *
*                                                                      *
* Last change 27.12.2006 by S. Roesler.                                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      PARAMETER ( LINP = 5 ,
     &            LOUT = 6 ,
     &            LDAT = 9 )

      PARAMETER (TINY10=1.0D-10,TINY14=1.0D-14,
     &           ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)

      PARAMETER (NCOMPX=20,NEB=8,NQB= 5,KSITEB=50)

      PARAMETER ( MAXNCL = 260,

     &            MAXVQU = MAXNCL,
     &            MAXSQU = 20*MAXVQU,
     &            MAXINT = MAXVQU+MAXSQU)

* Glauber formalism: parameters
      COMMON /DTGLAM/ RASH(NCOMPX),RBSH(NCOMPX),
     &                BMAX(NCOMPX),BSTEP(NCOMPX),
     &                SIGSH,ROSH,GSH,BSITE(0:NEB,NQB,NCOMPX,KSITEB),
     &                NSITEB,NSTATB

* Glauber formalism: cross sections
C      COMMON /DTGLXS/ ECMNN(NEB),Q2G(NQB),ECMNOW,Q2,
C     &                XSTOT(NEB,NQB,NCOMPX),XSELA(NEB,NQB,NCOMPX),
C     &                XSQEP(NEB,NQB,NCOMPX),XSQET(NEB,NQB,NCOMPX),
C     &                XSQE2(NEB,NQB,NCOMPX),XSPRO(NEB,NQB,NCOMPX),
C     &                XSDEL(NEB,NQB,NCOMPX),XSDQE(NEB,NQB,NCOMPX),
C     &                XETOT(NEB,NQB,NCOMPX),XEELA(NEB,NQB,NCOMPX),
C     &                XEQEP(NEB,NQB,NCOMPX),XEQET(NEB,NQB,NCOMPX),
C     &                XEQE2(NEB,NQB,NCOMPX),XEPRO(NEB,NQB,NCOMPX),
C     &                XEDEL(NEB,NQB,NCOMPX),XEDQE(NEB,NQB,NCOMPX),
C     &                BSLOPE,NEBINI,NQBINI

* Lorentz-parameters of the current interaction
      COMMON /DTLTRA/ GACMS(2),BGCMS(2),GALAB,BGLAB,BLAB,
     &                UMO,PPCM,EPROJ,PPROJ

* properties of photon/lepton projectiles
      COMMON /DTGPRO/ VIRT,PGAMM(4),PLEPT0(4),PLEPT1(4),PNUCL(4),IDIREC

* Glauber formalism: collision properties
      COMMON /DTGLCP/ RPROJ,RTARG,BIMPAC,
     &                NWTSAM,NWASAM,NWBSAM,NWTACC,NWAACC,NWBACC

* Glauber formalism: flags and parameters for statistics
C      LOGICAL LPROD
C      CHARACTER*8 CGLB
C      COMMON /DTGLGP/ JSTATB,JBINSB,CGLB,IOGLB,LPROD

      DIMENSION JS(MAXNCL),JT(MAXNCL)

      NTARG = ABS(NIDX)

      CALL DT_GCFDIAGR(NA,NB,IJPROJ,B,JS,JT,INTT,INTA,INTB,IDIREC,NIDX)
      IF (NIDX.LE.-1) THEN
         RPROJ = RASH(1)
         RTARG = RBSH(NTARG)
      ELSE
         RPROJ = RASH(NTARG)
         RTARG = RBSH(1)
      ENDIF

      RETURN
      END

* Stripped down version of DT_DIAGR for GCF
*
      SUBROUTINE DT_GCFDIAGR(NA,NB,IJPROJ,B,JS,JT,JNT,INTA,INTB,IDIREC,
     &                                                         NIDX)

************************************************************************
* Based on the original version by Shmakov et al.                      *
* This version dated 21.04.95 is revised by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      INCLUDE 'beagle.inc'

      PARAMETER ( LINP = 5 ,
     &            LOUT = 6 ,
     &            LDAT = 9 )

      PARAMETER (ZERO=0.0D0,TINY10=1.0D-10,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (TWOPI  = 6.283185307179586454D+00,
     &           PI     = TWOPI/TWO)

      PARAMETER ( MAXNCL = 260,
     &            MAXVQU = MAXNCL,
     &            MAXSQU = 20*MAXVQU,
     &            MAXINT = MAXVQU+MAXSQU)

* particle properties (BAMJET index convention)
      CHARACTER*8  ANAME
      COMMON /DTPART/ ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      PARAMETER (NCOMPX=20,NEB=8,NQB= 5,KSITEB=50)

* Glauber formalism: parameters
      COMMON /DTGLAM/ RASH(NCOMPX),RBSH(NCOMPX),
     &                BMAX(NCOMPX),BSTEP(NCOMPX),
     &                SIGSH,ROSH,GSH,BSITE(0:NEB,NQB,NCOMPX,KSITEB),
     &                NSITEB,NSTATB

* position of interacted nucleon
      DOUBLE PRECISION PosNuc,BX,BY
      COMMON /NPARINT/ PosNuc(4)

* VDM parameter for photon-nucleus interactions
      COMMON /DTVDMP/ RL2,EPSPOL,INTRGE(2),IDPDF,MODEGA,ISHAD(3)

C  obsolete cut-off information
      DOUBLE PRECISION PTCUT,PTANO,FPS,FPH,PSOMIN,XSOMIN
      COMMON /POCUT1/ PTCUT(4),PTANO(4),FPS(4),FPH(4),PSOMIN,XSOMIN
**

* coordinates of nucleons
      COMMON /DTNUCO/ PKOO(3,MAXNCL),TKOO(3,MAXNCL)

* event history
      PARAMETER (NMXHKK=200000)
      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* interface between Glauber formalism and DPM
      COMMON /DTGLIF/ JSSH(MAXNCL),JTSH(MAXNCL),
     &                INTER1(MAXINT),INTER2(MAXINT)

* statistics: Glauber-formalism
      COMMON /DTSTA3/ ICWP,ICWT,NCSY,ICWPG,ICWTG,ICIG,IPGLB,ITGLB,NGLB

* n-n cross section fluctuations
      PARAMETER (NBINS = 1000)
      COMMON /DTXSFL/ FLUIXX(NBINS),IFLUCT

      DIMENSION JS(MAXNCL),JT(MAXNCL),
     &          JS0(MAXNCL),JT0(MAXNCL,MAXNCL),
     &          JI1(MAXNCL,MAXNCL),JI2(MAXNCL,MAXNCL),JNT0(MAXNCL)
      DIMENSION NWA(0:210),NWB(0:210)

C...Added by liang 1/6/12
C...Switches for nuclear correction
      COMMON /PYNUCL/ INUMOD,CHANUM,ORDER,genShd
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER,genShd

C...Pythia eA shadowing common block from Mark 2017-06-30
      COMMON /PYSHAD/ NKNOTS,RDUMMY,RAVAL,XKNOT(100),YKNOT(100),
     &BKNOT(100),CKNOT(100),DKNOT(100),EKNOT(100),FKNOT(100),SHDFAC
      SAVE /PYSHAD/
      DOUBLE PRECISION XKNOT,YKNOT,BKNOT,CKNOT,DKNOT,EKNOT,FKNOT,RAVAL,
     &SHDFAC
      INTEGER NKNOTS,RDUMMY

C...External function Mark 2016-08-18
      DOUBLE PRECISION SHDMAP

      LOGICAL LFIRST
      DATA LFIRST /.TRUE./

      DATA NTARGO,ICNT /0,0/

C...Mark (08/13/2016) Pythia nuclear shadowing reroll b more often
      NTRMAX=500
      IF (GenShd.GE.2) NTRMAX=5
      NTARG = ABS(NIDX)

      IF (LFIRST) THEN
         LFIRST = .FALSE.
         IF (NCOMPO.EQ.0) THEN
            NCALL  = 0
            NWAMAX = NA
            NWBMAX = NB
            DO 17 I=0,210
               NWA(I) = 0
               NWB(I) = 0
   17       CONTINUE
         ENDIF
      ENDIF
      IF (NTARG.EQ.-1) THEN
         IF (NCOMPO.EQ.0) THEN
            WRITE(LOUT,*) ' DIAGR: distribution of wounded nucleons'
            WRITE(LOUT,'(8X,A,3I7)') 'NCALL,NWAMAX,NWBMAX = ',
     &                                NCALL,NWAMAX,NWBMAX
            DO 18 I=1,MAX(NWAMAX,NWBMAX)
               WRITE(LOUT,'(8X,2I7,E12.4,I7,E12.4)')
     &                          I,NWA(I),DBLE(NWA(I))/DBLE(NCALL),
     &                            NWB(I),DBLE(NWB(I))/DBLE(NCALL)
   18       CONTINUE
         ENDIF
         RETURN
      ENDIF

   16 CONTINUE
      PHIB=0.0D0

      NTRY = 0
    3 CONTINUE
      NTRY = NTRY+1
* initializations
      JNT  = 0
      DO 1 I=1,NA
         JS(I) = 0
    1 CONTINUE
      DO 2 I=1,NB
         JT(I) = 0
    2 CONTINUE
      IF (IJPROJ.EQ.7) THEN
         DO 8 I=1,MAXNCL
            JS0(I) = 0
            JNT0(I)= 0
            DO 9 J=1,NB
               JT0(I,J) = 0
    9       CONTINUE
    8    CONTINUE
      ENDIF

* nucleon configuration
C     IF ((NTARG.NE.NTARGO).OR.(MOD(ICNT,5).EQ.0)) THEN
      IF ((NTARG.NE.NTARGO).OR.(MOD(ICNT,1).EQ.0)) THEN
C        CALL DT_CONUCL(PKOO,NA,RASH,2)
C        CALL DT_CONUCL(TKOO,NB,RBSH(NTARG),1)
         IF (NIDX.LE.-1) THEN
            CALL DT_CONUCL(PKOO,NA,RASH(1),0)
            CALL DT_CONUCL(TKOO,NB,RBSH(NTARG),0)
         ELSE
            CALL DT_CONUCL(PKOO,NA,RASH(NTARG),0)
            CALL DT_CONUCL(TKOO,NB,RBSH(1),0)
         ENDIF
         NTARGO = NTARG
      ENDIF
      ICNT = ICNT+1

* LEPTO: pick out one struck nucleon
C...modified by liang to incorporate the PYTHIA model
C...and Mark (08/13/2016) Use ball Glauber for nuclear shadowing
C...and Mark (09/26/2016) to set impact parameter of photon
C...and Mark (06/30/2017) Reorganize and add partial shadowing
C... NOTE: SIGEFF =effective cross-section in mb
C      RAEVT = RAVAL
C      SIGEFF = SHDFAC*SHDMAP(RAVAL)
C      D2MAX = 0.1D0*SIGEFF/PI
C      IF (GenShd.LT.2 .OR. SIGEFF.LT.0.05) THEN
      JNT     = 1
      JS(1)   = 1
      IDX     = INT(DT_RNDM(X)*NB)+1
      JT(IDX) = 1
      BX      = TKOO(1,IDX)
      BY      = TKOO(2,IDX)
      BBEA    = DSQRT(BX*BX+BY*BY)
      PHIB    = PYANGL(BX,BY)
C     Count both nucleons as wounded.
      NCOLLT  = 2
      NCOLLI  = 2
C      RAEVT = 1.0
C      SIGEFF = 0.0
C
C      GOTO 19
C      ENDIF
C
C   19 CONTINUE
      INTA = 0
      INTB = 0
      DO 6 I=1,NA
        IF (JS(I).NE.0) INTA=INTA+1
    6 CONTINUE
      DO 7 I=1,NB
        IF (JT(I).NE.0) INTB=INTB+1
    7 CONTINUE
      ICWPG = INTA
      ICWTG = INTB
      ICIG  = JNT
      IPGLB = IPGLB+INTA
      ITGLB = ITGLB+INTB
      NGLB = NGLB+1

      IF (NCOMPO.EQ.0) THEN
         NCALL = NCALL+1
         NWA(INTA) = NWA(INTA)+1
         NWB(INTB) = NWB(INTB)+1
      ENDIF

      RETURN
      END
