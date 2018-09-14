*$ CREATE DPMJET.FOR
*COPY DPMJET
*
*
*===program dpmjet=====================================================*
*
      PROGRAM DPMJET

************************************************************************
* DPMJET 3.0: Main program.                                            *
* This version dated 10.03.95 is written by S. Roesler.                *
*                                                                      *
* Last change 08.01.2007 by S. Roesler.                                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER ( LINP = 5 ,
     &            LOUT = 6 ,
     &            LDAT = 9 )

* block data in DPMJET library (uncomment these declarations if library
* option is used)
C     EXTERNAL DT_BDEVAP,DT_BDNOPT,DT_BDPREE,DT_HADPRP,DT_BLKD46,
C    &         DT_BLKD47,DT_RUNTT,DT_NONAME,DT_ZK,DT_BLKD43

C     EXTERNAL PYDATA

* HP compiler settings
C     ON DOUBLE PRECISION UNDERFLOW IGNORE
C     ON DOUBLE PRECISION OVERFLOW IGNORE
C     ON DOUBLE PRECISION INEXACT IGNORE
C     ON DOUBLE PRECISION ILLEGAL IGNORE
C     ON DOUBLE PRECISION DIV 0 IGNORE
C     ON REAL UNDERFLOW IGNORE
C     ON REAL OVERFLOW IGNORE
C     ON REAL INEXACT IGNORE
C     ON REAL ILLEGAL IGNORE
C     ON REAL DIV 0 IGNORE
C     ON INTEGER OVERFLOW IGNORE

    1 CONTINUE
      CALL DT_INIT(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IDT,IGLAU)

      IF ((IDP.EQ. 3).OR.(IDP.EQ. 4).OR.(IDP.EQ. 5).OR.(IDP.EQ. 6).OR.
     &    (IDP.EQ.10).OR.(IDP.EQ.11)) THEN
* lepton-nucleus interactions
         CALL DT_LAEVT(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IDT,
     &                                                           IGLAU)
      ELSE
* hadron/nucleus/photon-nucleus interactions
         WRITE(*,*) 'BeAGLE does not support hadronic or real photo',
     &              'production collisions.'
C         CALL DT_AAEVT(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IGLAU)
      ENDIF

      GOTO 1

      END

*$ CREATE DT_USRHIS.FOR
*COPY DT_USRHIS
*
*===usrhis=============================================================*
*
      SUBROUTINE DT_USRHIS(MODE)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*
* COMMON /DTEVT1/ :
*                   NHKK         number of entries in common block
*                   NEVHKK       number of the event
*                   ISTHKK(i)    status code for entry i
*                   IDHKK(i)     identifier for the entry
*                                (for particles: identifier according
*                                 to the PDG numbering scheme)
*                   JMOHKK(1,i)  pointer to the entry of the first mother
*                                of entry i
*                   JMOHKK(2,i)  pointer to the entry of the second mother
*                                of entry i
*                   JDAHKK(1,i)  pointer to the entry of the first daughter
*                                of entry i
*                   JDAHKK(2,i)  pointer to the entry of the second daughter
*                                of entry i
*                   PHKK(1..3,i) 3-momentum
*                   PHKK(4,i)    energy
*                   PHKK(5,i)    mass
*
* event history
      PARAMETER (NMXHKK=200000)
      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
* extended event history
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHIST(2,NMXHKK)

C  general process information from phojet
      COMMON /POPRCS/ IPROCE,IDNODF,IDIFR1,IDIFR2,IDDPOM,IPRON(15,4)
      INTEGER IPROCE,IDNODF,IDIFR1,IDIFR2,IDDPOM,IPRON

      COMMON /POCKIN/ PTWANT,AS,AH,ALNS,ALNH,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,PTfin,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,
     &                QQPD,QQAL,PDF1(-6:6),PDF2(-6:6),WEIGHT,
     &                PHI1(5),PHI2(5),PHO1(5),PHO2(5),
     &                IA,IB,IC,ID,IV1,IV2,MSPR,IREJSC
      INTEGER IA,IB,IC,ID,IV1,IV2,MSPR,IREJSC
      DOUBLE PRECISION PTWANT,AS,AH,ALNS,ALNH,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                 PT,PTfin,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,QQPD,QQAL,
     &                 PDF1,PDF2,WEIGHT,PHI1,PHI2,PHO1,PHO2

C  added for the hard scattering information from phojet
      INTEGER MSCAHD
      PARAMETER ( MSCAHD = 50 )
      INTEGER LSCAHD,LSC1HD,LSIDX,
     &        NINHD,N0INHD,NIVAL,N0IVAL,NOUTHD,NBRAHD,NPROHD
      DOUBLE PRECISION PPH,PTHD,ETAHD,Q2SCA,PDFVA,XHD,VHD,X0HD
      COMMON /POHSLT/ LSCAHD,LSC1HD,LSIDX(MSCAHD),
     &                PPH(8*MSCAHD,2),PTHD(MSCAHD),ETAHD(MSCAHD,2),
     &                Q2SCA(MSCAHD,2),PDFVA(MSCAHD,2),
     &                XHD(MSCAHD,2),VHD(MSCAHD),X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2),N0INHD(MSCAHD,2),
     &                NIVAL(MSCAHD,2),N0IVAL(MSCAHD,2),
     &                NOUTHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)

      COMMON /DIFFTS/ TTtest
      DOUBLE PRECISION TTtest
      DOUBLE PRECISION tHAT

C  added by liang to search IP from phojet( particle combination )
      INTEGER IPSTATUS
      COMMON /PARCMB/ IPSTATUS


* Glauber formalism: collision properties
      COMMON /DTGLCP/ RPROJ,RTARG,BIMPAC,
     &                NWTSAM,NWASAM,NWBSAM,NWTACC,NWAACC,NWBACC

* properties of interacting particles
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE

* specify output file unit added by liang
      COMMON /OUTPT/ FOUT
      INTEGER FOUT

* out e- and gamma* information
      COMMON /DTGPRO/ VIRT,PGAMM(4),PLEPT0(4),PLEPT1(4),PNUCL(4),IDIREC
      COMMON /DTLGVX/ PPL0(4),PPL1(4),PPG(4),PPA(4)

* added by liang to check the photon flux 12/28/11
      COMMON /FLCHK/ PFXCHK
      DOUBLE PRECISION PFXCHK

* added by liang to store the output event variables 1/20/12
      COMMON /EVTOUT/ XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT
      DOUBLE PRECISION XBJOUT,YYOUT,W2OUT,NUOUT,Q2OUT

* Removed by Mark. 11/10/16. Not really used and not really correct for A.
* added by liang to read beam information to transform back to lab frame
* 1/20/12
C      COMMON /BEAMIF/ pBeam, eBeampr
C      DOUBLE PRECISION pBeam, eBeampr

* some dummy values for outgoing electron and radiated photon
      INTEGER EDUMMY1,EDUMMY2,EDUMMY3,EDUMMY4,EDUMMY5,EDUMMY6,EDUMMY7
      INTEGER GDUMMY1,GDUMMY2,GDUMMY3,GDUMMY4,GDUMMY5,GDUMMY6,GDUMMY7
      DOUBLE PRECISION EDUMMY8
      DOUBLE PRECISION GDUMMY8

      DATA EDUMMY1,EDUMMY2,EDUMMY3,EDUMMY4,EDUMMY5,EDUMMY6,
     &     EDUMMY7,EDUMMY8 
     &     /-1,0,11,0,0,0,0,0.000511/
      DATA GDUMMY1,GDUMMY2,GDUMMY3,GDUMMY4,GDUMMY5,GDUMMY6,
     &     GDUMMY7,GDUMMY8 
     &     /0,0,22,0,0,0,0,0.0/

      GOTO (1,2,3) MODE

*------------------------------------------------------------------
*
    1 CONTINUE
*
* initializations
*
*  Called with MODE=1 once at the beginning of the run.
*
      RETURN
*
*------------------------------------------------------------------
*
    2 CONTINUE

      WRITE(*,*) "Error: CALL DT_USRHIS(2) is no longer supported."
      WRITE(*,*) "MCGENE should be set to 5 - Pythia"

C*
C* scoring of the present event
C*
C*  Called with MODE=2 every time one event has been finished.
C*
C*  The final state particles from the actual event (number NEVHKK)
C*  can be found in DTEVT1 and identified by their status:
C*
C*     ISTHKK(i) = 1    final state particle produced in
C*                      photon-/hadron-/nucleon-nucleon collisions or
C*                      in intranuclear cascade processes
C*                -1    nucleons, deuterons, H-3, He-3, He-4 evaporated
C*                      from excited nucleus and
C*                      photons produced in nuclear deexcitation processes
C*                1001  residual nucleus (ground state)
C*
C*  The types of these particles/nuclei are given in IDHKK as follows
C*
C*     all final state part. except nuclei :
C*       IDHKK(i)=particle identifier according to PDG numbering scheme
C*     nuclei (evaporation products, and residual nucleus) :
C*       IDHKK(i)=80000, IDRES(i)=mass number, IDXRES(i)=charge number
C*
C*  The 4-momenta and masses can be found in PHKK (target nucleus rest frame):
C*                   PHKK(1..3,i) 3-momentum (p_x,p_y,p_z)
C*                   PHKK(4,i)    energy
C*                   PHKK(5,i)    mass
C*
C*
C*
C*  Pick out the final state particles from DTEVT1 in each event for
C*  instance by the following loop (NHKK=number of entries in the present
C*  event) and fill your histograms
CC     DO 20 I=1,NHKK
CC        IF (ABS(ISTHKK(I)).EQ.1) THEN
CC        ELSEIF (ABS(ISTHKK(I)).EQ.1001) THEN
CC        ENDIF
CC  20 CONTINUE
C
Cc...deal with some kinematics issue
Cc  get the hard interaction information      
C      NF1=0
C      NF2=0
C      XPART1=0.0D0
C      XPART2=0.0D0      
C      NUCL=0
C      DO I=1,LSCAHD      
C      NF1=NINHD(I,1)
C      NF2=NINHD(I,2)
C      XPART1=XHD(I,1)
C      XPART2=XHD(I,2)
C      NUCL=NBRAHD(I,2)
C      ENDDO
C      IF(IPROCE.EQ.1) THEN
C         tHat = VHD(1)
C      ELSE IF(IPROCE.EQ.3.OR.IPROCE.EQ.5.OR.IPROCE.EQ.6) THEN
C         tHAT = TTtest
C      ENDIF
C* Lorentz transformation for e- boost from the n rest to real lab frame
C      IF(IJPROJ.EQ.7) THEN
C         P = pBeam
C         E = SQRT(P*P+0.938*0.938)
C         GMA = E/0.938
C         BGT = P/E
C         CALL DT_DALTRA(GMA, 0.0D0,0.0D0,-BGT*GMA,PPL1(1),PPL1(2),
C     &   PPL1(3),PPL1(4),PLTOT,PL1,PL2,PL3,PL4)
C******calculate the theta of out e-*******************************      
C         THETA=ACOS(PL3/PLTOT)
C******transform incoming e- into lab frame************************
C         eBeamE=SQRT(eBeampr**2+EDUMMY8**2)
C         CALL DT_DALTRA(GMA, 0.0D0,0.0D0,-BGT*GMA,0.0D0,0.0D0,eBeampr,
C     &   eBeamE,eDummy,eBeamx,eBeamy,eBeamz,eBeam0)      
C      ENDIF
C
C
Cc...event file header 
C      IF(NEVHKK.EQ.1) THEN
C         WRITE(FOUT,*)'DPMJET EVENT FILE'
C         WRITE(FOUT,*)'==========================================='
C
C         IF(IJPROJ.EQ.1) THEN
C            WRITE(FOUT,*)'ievent, process1, process2, IP, N_coll,',
C     &         ' N_par1, N_par2, nTracks, B_impact'
C         ELSEIF(IJPROJ.EQ.7) THEN
C            WRITE(FOUT,*)'I, ievent, process1, process2, IP, W2, nu,
C     &   Q2, x, y, theta_e, photonFlux, targetparton, 
C     &   projectileparton,xtargetparton, xprojectileparton, tHat, 
C     &   bimpact, N_coll, N_par, nucleon, nrTracks'
C         ENDIF
C
C         WRITE(FOUT,*)'==========================================='
C         IF(IJPROJ.EQ.1) THEN
C            WRITE(FOUT,*)'I  ISTHKK(I)  IDHKK(I) IDRES(I) IDXRES(I)
C     &   JMOHKK(1,I)  JDAHKK(1,I)  JDAHKK(2,I)  PHKK(1,I)  PHKK(2,I)  
C     &   PHKK(3,I)  PHKK(4,I)  PHKK(5,I)  VHKK(1,I) VHKK(2,I) VHKK(3,I)'
C         ELSEIF(IJPROJ.EQ.7) THEN
C            WRITE(FOUT,*)'I  ISTHKK(I)  IDHKK(I)  JMOHKK(1,I) 
C     &   JDAHKK(1,I)  JDAHKK(2,I)  PHKK(1,I)  PHKK(2,I)  PHKK(3,I)
C     &   PHKK(4,I)  PHKK(5,I)  VHKK(1,I) VHKK(2,I) VHKK(3,I)'
C         ENDIF
C         WRITE(FOUT,*)'==========================================='
C      ENDIF
C
Cc...event wise variables output      
C      IF(IJPROJ.EQ.1) THEN
C         WRITE(FOUT,'(8I11,F17.5)') NEVHKK,IPROCE,MSPR,IPSTATUS,
C     &   NWTSAM,NWASAM,NWBSAM,NHKK,BIMPAC
C      ELSEIF(IJPROJ.EQ.7) THEN
Cc...0 initial to be recongonized by eic root maker code
C         WRITE(FOUT,'(5I8,7F17.5,2I6,4F17.5,4I6)') 0,NEVHKK,IPROCE,
C     &   MSPR,IPSTATUS,W2OUT,NUOUT,Q2OUT,XBJOUT,YYOUT,THETA,PFXCHK,
C     &   NF2,NF1,XPART2,XPART1,tHat,BIMPAC,NWTSAM,NWBSAM,NUCL,3+NHKK
C      ENDIF
C      WRITE(FOUT,*)'==========================================='
C
C
Cc...particle wise varaibles output
C      DO  I=1,NHKK
C
C         !initialize output index
C         INDX=I
Cc...add beam particle and out e- info for ep mode output structure      
C         IF(IJPROJ.EQ.7) THEN
Cc...modify the mother/daugther relation when more info put in
C            DO J=1,2
C               IF(JMOHKK(J,I).GT.0) THEN
C                  JMOHKK(J,I)=JMOHKK(J,I)+4
C               ENDIF
C               IF(JDAHKK(J,I).GT.0) THEN
C                  JDAHKK(J,I)=JDAHKK(J,I)+4
C               ENDIF
C            ENDDO
C            IF(I.EQ.1) THEN
C               WRITE(FOUT,996) 1,21,11,EDUMMY4,EDUMMY6,EDUMMY7,
C     &         eBeamx,eBeamy,-eBeamz,eBeam0,EDUMMY8,
C     &         0.0D0,0.0D0,0.0D0
C     &         ,0,0
C
C               WRITE(FOUT,996) 2,21,2212,EDUMMY4,
C     &         EDUMMY6,EDUMMY7,0.0D0,0.0D0,P,E,PPA(4),0.0D0,0.0D0,0.0D0
C     &         ,0,0
C
Cc...output for scattered e-
Cc...mother1 set as 3 to be recongonized by eic root maker code
C               WRITE(FOUT,996) 3,1,11,3,
C     &         EDUMMY6,EDUMMY7,PL1,PL2,-PL3,PL4,EDUMMY8,0.0D0,0.0D0,
C     &         0.0D0 ,0,0
Cc...transform gamma to lab frame
C               CALL DT_DALTRA(GMA, 0.0D0,0.0D0,-BGT*GMA,PPG(1),PPG(2),
C     &         PPG(3),PPG(4),PDUMMY,P1,P2,P3,P4)
Cc...output for gamma
C               WRITE(FOUT,996) I+3,21,IDHKK(I),1,
C     &         JDAHKK(1,I),JDAHKK(2,I),P1,P2,-P3,P4,GDUMMY8,0.0D0,0.0D0,
C     &         0.0D0 ,0,0
C            ENDIF
C
Cc...transform final particles from N-rest to Lab frame for ep mode
C            IF((ABS(ISTHKK(I)).EQ.1).OR.(ISTHKK(I).EQ.1000).OR.
C     &                                  (ISTHKK(I).EQ.1001)) THEN
C               CALL DT_DALTRA(GMA, 0.0D0,0.0D0,-BGT*GMA,PHKK(1,I),
C     &         PHKK(2,I),PHKK(3,I),PHKK(4,I),PDUMMY,P1,P2,P3,P4)
C               PHKK(1,I)=P1
C               PHKK(2,I)=P2
C               PHKK(3,I)=-P3
C               PHKK(4,I)=P4
C            ENDIF
C            !reassign output index
C            INDX=I+4
C
C         ENDIF
C
C         !add extended info, IDRES, IDXRES and so on for AA mode
C         IF(IJPROJ.EQ.1) THEN
C            WRITE(FOUT,995) INDX,ISTHKK(I),IDHKK(I),IDRES(I),IDXRES(I),
C     &     JMOHKK(1,I),JDAHKK(1,I),JDAHKK(2,I),
C     &     PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
C     &     PHKK(5,I),VHKK(1,I),VHKK(2,I),VHKK(3,I)
C 995        FORMAT(2I5,6I9,5F17.5,3E15.6)
C         ELSEIF(IJPROJ.EQ.7) THEN
C            WRITE(FOUT,996) INDX,ISTHKK(I),IDHKK(I),
C     &     JMOHKK(1,I),JDAHKK(1,I),JDAHKK(2,I),
C     &     PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
C     &     PHKK(5,I),VHKK(1,I),VHKK(2,I),VHKK(3,I)
C     &     ,IDRES(I), IDXRES(I)
CC 996        FORMAT(2I5,4I9,5F17.5,3E15.6)
C 996        FORMAT(2I5,4I9,5F17.5,3E15.6,2I9)
C         ENDIF
C         !as VHKK in mm, very small, use scientific precision
C         !to be noticed, VHKK is the variable used for coordinates
C         !in target rest frame by incorporating bimpact in x direction
C         !which means VHKK(1,I)=x(in target rest) + bimpact
C     ENDDO
C      WRITE(FOUT,*)'=============== Event finished ==============='
C
C
C*  At any time during the run a list of the actual entries in DTEVT1 and
C*  DTEVT2 can be obtained (output unit 6) by the following statement:
CC     CALL DT_EVTOUT(4)
C
C      RETURN
*
*------------------------------------------------------------------
*
    3 CONTINUE
*
* output/statistics/histograms etc.
*
*  Called with MODE=3 once after all events have been sampled.
*
      RETURN

      END
