*These are the interfaces to GCF input file needed by dpmjet
*
* Written by  Mark Baker  2019-04-27
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

      IF (IFSEED.EQ.0) THEN 
         call idate(today)      ! today(1)=day, (2)=month, (3)=year
         call itime(now)        ! now(1)=hour, (2)=minute, (3)=second
         initseed = today(1)+10*today(2)+today(3) + now(1)+5*now(3)
      ELSE
         initseed=IFSEED
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

C     GenNucDens is used even without quenching.
C      call GenNucDens(ITZ, IT)

C      !changed by liang to fit with INC particle status
C      CALL DT_PYDECY(1)

      RETURN
      END
*=====dt_gcfevntqe========================================================
* read in the input file from GCF and dump this event to DTEVT1
* GCF is in the target/ion rest frame, while dpmjet makes many of it's 
* calculations in the photon-proton c.m.s frame. So when we are copying 
* the data common to dpmjet, we must make the right Lorentz transformation. 
* Plus, the virtual photon has to be directed to the z+ direction. A 
* rotation for the reference is also necessary.
*
      SUBROUTINE DT_GCFEVNTQE(Q2,YY)
 
*     output:
*           Q2    Q2 of this current event 
*           YY    Y of this current event  

      IMPLICIT NONE
C      include "mc_set.inc"
C      include "mcRadCor.inc"
C      include "radgen.inc"
C      include "phiout.inc"
      include "beagle.inc"
      include "bea_pyqm.inc"

      EXTERNAL IDT_ICIHAD
      INTEGER IDT_ICIHAD

C* flags for input different options
C      LOGICAL LEMCCK,LHADRO,LSEADI,LEVAPO
C      COMMON /DTFLG1/ IFRAG(2),IRESCO,IMSHL,IRESRJ,IOULEV(6),
C     &                LEMCCK,LHADRO(0:9),LSEADI,LEVAPO,IFRAME,ITRSPT

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
      REAL ZMIN, ZTRY
      DOUBLE PRECISION EEFIX,PZEFIX,BETAFX,GAMMFX

Cc...added by liang & Mark to include pythia energy loss datas
      double precision PAUX, DPF
      COMMON /PFAUX/ PAUX(4), DPF(4)

      DOUBLE PRECISION P1, P2, P3, P4, PTOT, PP1, PP2, PP3, PP4
      INTEGER QCHG

      DOUBLE PRECISION Q2, YY, XX, XXALT

* added by liang to store the output event variables 1/20/12
      COMMON /EVTOUT/ XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT
      DOUBLE PRECISION XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT

C     Extensions for BEAGCF - variables that we would look up in a 
C     Pythia common blocks for MCGENE=5. Note: RECTYPE will be 
C     reported as process.
      INTEGER RECTYPE, LEADTYPE
      DOUBLE PRECISION LEPTONPHI
      COMMON /BEAGCF/ RECTYPE, LEADTYPE, LEPTONPHI
C  Local

      INTEGER ITRK, NINTS, IIMAIN, ILEPT

      integer LINP
      parameter ( LINP=28 )
      
      INTEGER IDUM, IEVENT, IDIM, NDIM
      INTEGER nrTracks
      parameter ( NDIM=4 )
      
      DOUBLE PRECISION PFERR, PPSUM(NDIM)

      CALL DT_EVTINI
      READ(LINP,*) IDUM, IEVENT, ltype, IT, ITZ, PZLEP, RECTYPE, 
     &     LEADTYPE, Q2OUT, XBJOUT, NUOUT, LEPTONPHI, PXF, PYF, PZF,
     &     EEXC(2), RAEVT, nrTracks
      IF (IEVENT.LE.5) WRITE(*,*) RAEVT

      IF (IDUM.NE.0 .OR. nrTracks.NE.9 .OR. EEXC(2).LT.-TINY10 .OR.
     &     IEVENT.NE.NEVENT) STOP
     &     'DT_GCFEVNTQE: Bad event header format in GCF input file.'

C     Protect against slightly negative roundoff errors for E*=0.
      IF (EEXC(2).LT.0.0D0) EEXC(2)=0.0D0
      YYOUT = NUOUT/ABS(PZLEP)
C     W2 doesn't make physical sense for QE. Just M2 for pF=0.
C     Naively: W2OUT = 2.0D0*MNGCF*NUOUT - Q2OUT + MNGCF*MNGCF
      W2OUT = 0.0D0

C... Collector to get the momentum and charge sums 
C... Note: for now, treat the (A-2)* as stable.
      DO IDIM=1,NDIM
         PPSUM(IDIM)=0.0D0
      ENDDO
      QCHG=0
      WRITE(*,*)'I,ISTHKK,IDHKK,JMOHKK(1-2),JDAHKK,PHKK(1-5),'//
     & 'IDRES,IDXRES,NOBAM'

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

CC Thickness is twice the integral from 0 to infinity
C      THKB = 2.0D0*DCALC(0.0D0)
C      IDUM = 0
C      THKSCL = THKB*SCLFAC(IDUM)

C     Calculate distance travelled in the nucleus
C      IF (NINTS.EQ.1) THEN
C         DAVG = DCALC(VHKK(3,IIMAIN)/FM2MM)
C         DFIRST = DAVG
C      ELSE
C         ZMIN = 10000.0
C         DO III=1,NINTS       
C            ZTRY=VHKK(3,IINTER(III))
C            IF (ZTRY.LT.ZMIN) THEN
C               ZMIN=ZTRY
C            ENDIF
C         ENDDO
C         DFIRST = DCALC(ZMIN/FM2MM)
C         DAVG = DCALC(VHKK(3,IIMAIN)/FM2MM)
C      ENDIF

c     WRITE(*,*) 'Calculated D1st, Davg: ',DFIRST,' ',DAVG

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

C    The lines below only come into play when we want to use the HCMS
C... Note: DPF(mu) = P(mu)_true - P(mu)_naive is a 4-momentum too.   
C    DPF is the name in the HCMS
C    PXF,PYF,PZF,EKF in the TRF
C      DPF(1) = PXF
C      DPF(2) = PYF
C      CALL DT_LTNUC(PZF,EKF,DPF(3),DPF(4),3)

C  Unlike Pythia, GCF event is already in IRF/TRF, same frame as PHKK.
C  For GCF, we will use the convention that -z is along the incoming electron
C  rather than the convention we used for Pythia that +z is along gamma*.


CC MDB 2017-08-07 Rename PF as PYQREC.
CC      Boost directly from TRF z=g* to HCMS z=g* in one shot...
CC      WRITE(*,*) 'Boost from TRF to HCMS. gamma, betagamma:',
CC     &      GACMS(2),BGCMS(2)
C      P3=PYQREC(3)
C      P4=PYQREC(4)
C      PYQREC(3)=GACMS(2)*P3-BGCMS(2)*P4
C      PYQREC(4)=GACMS(2)*P4-BGCMS(2)*P3
C      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
C         WRITE(*,*) 'Mycalc Recoil PYQREC(1-4) HCMS: ',PYQREC(1),
C     &        PYQREC(2),PYQREC(3),PYQREC(4)
C      endif
C      CALL DT_LTNUC(P3,P4,PYQREC(3),PYQREC(4),3)
C      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
C         WRITE(*,*) 'DT_LTNUC Recoil PYQREC(1-4) HCMS: ',PYQREC(1),
C     &        PYQREC(2),PYQREC(3),PYQREC(4)
C      endif
Cc...Translate PYJETS into HEPEVT event record
C      CALL PYHEPC(1)

c...Output part
C      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
C         WRITE(*,*) 'Listing HEPEVT as we go. Before rotation which',
C     &        'takes x->-x and z->-z.'
C         WRITE(*,*) 'J, ISTHEP(J), IDHEP(J), JMOHEP(1-2,J), ',
C     &        'JDAHEP(1-2,J)','PHEP(1-5,J)'
C      endif
Cc...First loop to find exchanged boson

      IF (IDHKK(5).NE.22) STOP 'DT_GCFEVNTQE: Bad event format'
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
C*  in the n rest frame rotate virtual photon  angles to +z axis
C*...COD=cos(theta) SID=sin(theta) COF=cos(phi) SIF=sin(phi)
C      call DT_DALTRA(GAA,eveBETA(1),eveBETA(2),eveBETA(3),
C     &GAMM(1),GAMM(2),GAMM(3),GAMM(4),PTOT,P1,P2,P3,P4)
C      PTOT=SQRT(P1**2+P2**2+P3**2)
C      COD = P3/PTOT
C      PPT = SQRT(P1**2+P2**2)
C      SID = PPT/PTOT
C      IF(P1.GT.ZERO) THEN
C         COF = ONE
C      ELSE
C         COF = -ONE
C      ENDIF
C      SIF = ZERO      
C      IF (PTOT*SID.GT.TINY10) THEN
C         COF = P1/(SID*PTOT)
C         SIF = P2/(SID*PTOT)
C         ANORF = SQRT(COF*COF+SIF*SIF)
C         COF = COF/ANORF
C         SIF = SIF/ANORF
C      ENDIF
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
      if(IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5)) then
         print*,'4-momentum & charge totals for Status=1 & 1000:'
         print*,'PP1=',PPSUM(1),' PP2=',PPSUM(2),' PP3=',PPSUM(3),
     &          ' PP4=',PPSUM(4),' QCHG=',QCHG
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

      IMPLICIT NONE
C      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

C      include 'pythia.inc'              ! All PYTHIA commons blocks
C      include "mc_set.inc"
C      include "py6strf.inc"
C      include "mcRadCor.inc"
C      include "radgen.inc"
C      include "phiout.inc"
      include "beagle.inc"

C      EXTERNAL PYCHGE, NBARY, PYMASS
C      INTEGER  NBARY

      EXTERNAL AZMASS, PYMASS
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

C* lorentz transformation parameter
C      COMMON /LTPARA/ BGTA(4), GAMM(5), eveBETA(4), GAA
C      COMMON /ROTATE/ COF,COD,SIF,SID

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

C* flags for input different options
C      LOGICAL LEMCCK,LHADRO,LSEADI,LEVAPO
C      COMMON /DTFLG1/ IFRAG(2),IRESCO,IMSHL,IRESRJ,IOULEV(6),
C     &                LEMCCK,LHADRO(0:9),LSEADI,LEVAPO,IFRAME,ITRSPT

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
      INTEGER IZSUM,IASUM,ILEPT,IBOSON,J,IDIM,NDIM
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
C********rotate back from the gamma* +z direction***********
C            P1=COF*(COD*PHKK(1,J)+SID*PHKK(3,J))+SIF*PHKK(2,J)
C            P2=SIF*(COD*PHKK(1,J)+SID*PHKK(3,J))-COF*PHKK(2,J)
C            P3=COD*PHKK(3,J)-SID*PHKK(1,J)
C            P4=PHKK(4,J)
Cc**************transform back to the lab frame**************      
C            call DT_DALTRA(GAA,-eveBETA(1),-eveBETA(2),-eveBETA(3),
C     &           P1,P2,P3,P4,PTOT,PP1,PP2,PP3,PP4)
Cc     WRITE(89,996) J,ISTHKK(J),IDHKK(J),JMOHKK(1,J),
Cc     &          JMOHKK(2,J),JDAHKK(1,J),JDAHKK(2,J),
Cc     &          PP1,PP2,-PP3,PP4, !remember to change the sign of z back
Cc     &              PHKK(5,J)
Cc  996      FORMAT(I5,I5,I8,4I5,5F17.5)
C            PHKK(1,J)=PP1
C            PHKK(2,J)=PP2      
C            PHKK(3,J)=-PP3
C            PHKK(4,J)=PP4
Cc     PHKK(3,J)=pgamma*(P3+pbeta*P4)
Cc     PHKK(4,J)=pgamma*(P4+pbeta*P3)
c...find the exchanged boson and out e- to make it fit root tree making rules
c...in the following steps
            IF ((ISTHKK(J).EQ.3).AND.(IDHKK(J).EQ.22.OR.IDHKK(J).EQ.23)
     &           .AND. (JMOHKK(1,J).EQ.(NPOINT(1)+1))) THEN
               IBOSON=J
            ELSEIF (ISTHKK(J).EQ.1.AND.JMOHKK(1,J).EQ.4.AND.
     &              (ABS(IDHKK(J)).EQ.11.OR.ABS(IDHKK(J)).EQ.13)) THEN
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
                  IF (IDXRES(J).NE.0) NINCCH = NINCCH + 1
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
      DO I=1,3
         IF(IDHKK(I).NE.2112 .AND. IDHKK(I).NE.2212 .AND. 
     &        IDHKK(I).NE.80000) THEN
            WRITE(*,*) 'DT_GCFOUTQE ERROR: IN-NUCLEON IDHKK=',IDHKK(I)
         ELSEIF (ISTHKK(I).EQ.18 .OR. ISTHKK(I).EQ.12) THEN
            NWND = NWND + 1
            IF (IDHKK(I).EQ.2212) NWDCH = NWDCH + 1
         ELSEIF (ISTHKK(I).NE.14) THEN
            WRITE(*,*) 'DT_GCEOUTQE ERROR: IN-NUCLEON ISTHKK=',ISTHKK(I)
         ENDIF
      ENDDO
      IF (USERSET.EQ.5) THEN         
         USER2 = SIGEFF
      ELSEIF (USERSET.GE.8) THEN         
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
            CALL EVTSUM(P5SUM,IZSUM,IASUM,VERBOSE)
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
      write(29,33) 0, ievent, genevent, ltype, it, itz, pzlep, 
     &        pztarg, pznucl, crang, crori, RECTYPE, LEADTYPE,
     &        LEADTYPE, 1.0D0, 22, 0.0D0, 0.0D0,
     &        YYOUT, Q2OUT, XBJOUT, W2OUT, NUOUT,
     &        LEPTONPHI, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0, 0.0,
     &        0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &        0.0D0, BBEA, PHIB, THKB, THKSCL, NCOLLT, NCOLLI,
     &        NWND, NWDCH,
     &        NNEVAP, NPEVAP, AREMN, NINC, NINCCH, DFIRST, DAVG,
     &        PXF, PYF, PZF, EEXC(2), RAEVT, USER1, USER2, USER3, 
     &        nrtrack+4 
C Note: Use E rather than F format for GCF RAEVT (weight)
 33   format((I4,1x,$),(I10,1x,$),4(I4,1x,$),4(f12.6,1x,$),3(I4,1x,$),
     &     I6,1x,$,f9.6,1x,$,I6,1x,$,2(f12.6,1x,$),7(f18.11,3x,$),
     &     11(f19.9,3x,$),4(f10.6,1x,$),9(I5,1x,$),2(f10.6,1x,$),
     &     3(f15.6,1x,$),f12.6,1x,$,4(e17.8,1x,$),I6,/)
      write(29,*)'============================================'

***************standard output for particle info************************
c...add 2 beam information at first to fit into root tree making rule      
c... MDB Change these lines to use the correct pid for lepton+nucleon!
      I=NPOINT(1)+1   
C     Incoming lepton
      write(29,34) 1,21,IDHKK(I),0,0,0,0,
     &     PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),PHKK(5,I),
     &     VHKK(1,I),VHKK(2,I),VHKK(3,I)
     &     ,0,0,0
c...  Incoming nucleon for eicsmear. Use MNGCF for consistent kinematics.
      write(29,34) 2,21,IDHKK(1),0,0,0,0,0.0D0,0.0D0,0.0D0,MNGCF,MNGCF,
     &        VHKK(1,1),VHKK(2,1),VHKK(3,1),1,IDXRES(1),0
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
C         IF (IDHKK(I).EQ.80000) THEN
            ZOUT = IDXRES(I)
            AOUT = IDRES(I)
C         ELSE
C            AOUT = NBARY(IDHKK(I))
C            IF (MOD(AOUT,3).EQ.0) THEN
C               AOUT = AOUT/3
C            ELSE
C               AOUT = 0
C            ENDIF
C            ZOUT = PYCHGE(IDHKK(I))
C            IF (MOD(ZOUT,3).EQ.0) THEN
C               ZOUT = ZOUT/3
C            ELSE
C               ZOUT = 0
C            ENDIF
C         ENDIF
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
     &     LEPTONPHI, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0, 0.0,
     &     0.0D0, 0.0D0, 0.0D0, 0.0D0,
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

C* flags for input different options
C      LOGICAL LEMCCK,LHADRO,LSEADI,LEVAPO
C      COMMON /DTFLG1/ IFRAG(2),IRESCO,IMSHL,IRESRJ,IOULEV(6),
C     &                LEMCCK,LHADRO(0:9),LSEADI,LEVAPO,IFRAME,ITRSPT

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
      LOGICAL LLCPOT
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



