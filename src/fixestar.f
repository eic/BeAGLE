      SUBROUTINE FIXESTAR(ESTI,ESTF,AMRCL0,PRCL,VERB)
C*********************************************************************
 
C... FIXESTAR - Mark D. Baker 2023-06-08
C...
C... We should be working in the naive HCMS of the e+N collision
C... with gamma* momentum along +z and incoming ion along -z  
C... 
C... Change the P+ and M of the remnant to fix E*, while leaving P-,Pt alone
C... For FSP (Final State Particles excluding e' & remnant) , change total
C... P+ and s to conserve overall P+ while leaving total P-,Pt alone
C...
C... Remnant pmu is PRCL(2,1-4), mass AMRCL(2), gsmass AMRCL0(2)
C...
C... The scattered electron should still be status 99 at this point.
C...
      IMPLICIT NONE
C     Calling params.: in    in    in         I/O
      DOUBLE PRECISION ESTI, ESTF, AMRCL0(2), PRCL(2,4)
C              in
      LOGICAL VERB
* event history COMMON

      INTEGER NMXHKK,NHKK,NEVHKK,ISTHKK,IDHKK,JMOHKK,JDAHKK
      DOUBLE PRECISION PHKK,VHKK,WHKK
      PARAMETER (NMXHKK=200000)

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

C     Note: MAXTRY, MAXPRTS should match those in SCLSUBEVT
C     Local variables
      INTEGER MAXTRY, MAXPRTS, NPARTS, NSCLTR
      PARAMETER (MAXTRY=10, MAXPRTS=1000)
      INTEGER IDXHKK(MAXPRTS)
      LOGICAL SCLFAIL
      DOUBLE PRECISION PFSP(5),M2FSP,PFSPPL,PFSPMI,MRMN,PRMNPL,PRMNMI
      DOUBLE PRECISION PFIN(5),M2FIN,DELTA,PFINPL,PFINMI
      DOUBLE PRECISION PTMP(MAXPRTS,5),PTMPTOT
      DOUBLE PRECISION EPS, W2OUT, GA, BGX, BGY, BGZ
      PARAMETER (EPS=1.0D-9)
      INTEGER IP, MU, IDX

      DO MU=1,4
         PFSP(MU)=0.0D0
      ENDDO

      WRITE(*,*) 'In FIXESTAR'
      WRITE(*,*) 'IP, ISTHKK(IP), IDHKK(IP), PHKK(1,IP), PHKK(2,IP)',
     &        ' PHKK(3,IP), PHKK(4,IP), PHKK(5,IP)'
      DO IP=1,NHKK
         WRITE(*,*) IP, ISTHKK(IP), IDHKK(IP), PHKK(1,IP), PHKK(2,IP),
     &        PHKK(3,IP), PHKK(4,IP), PHKK(5,IP)
         IF (ISTHKK(IP).EQ.1) THEN
            DO MU=1,4
               PFSP(MU) = PFSP(MU) + PHKK(MU,IP)
            ENDDO
         ENDIF
      ENDDO
      PFSPPL = PFSP(4) + PFSP(3)
      PFSPMI = PFSP(4) - PFSP(3)
      M2FSP = PFSPPL*PFSPMI-PFSP(1)**2-PFSP(2)**2
      PFSP(5) = DSQRT(M2FSP)

      WRITE(*,*) 'PFSP(1-5):',PFSP(1),PFSP(2),PFSP(3),PFSP(4),PFSP(5)
      WRITE(*,*) 'PFSPPL,PFSPMI:',PFSPPL,PFSPMI
      PRMNPL = PRCL(2,4) + PRCL(2,3)
      PRMNMI = PRCL(2,4) - PRCL(2,3)
      MRMN = DSQRT(PRMNPL*PRMNMI-PRCL(2,2)**2-PRCL(2,1)**2)
      WRITE(*,*) 'MRMN, ESTI+AMRCL0(2)',MRMN, ESTI+AMRCL0(2)
      IF (ABS(MRMN-ESTI-AMRCL0(2)).LT.1D-09) STOP ('FIXESTAR FATAL')

C     Required change in P+remn with fixed P-remn,Ptremn for ESTF
      DELTA = 0.5D0*(ESTF+ESTI+2.0D0*AMRCL0(2))*(ESTF-ESTI)/PRMNMI
      PRCL(2,4) = PRCL(2,4)+DELTA
      PRCL(2,3) = PRCL(2,3)+DELTA
      PFIN(4) = PFSP(4)-DELTA
      PFIN(3) = PFSP(3)-DELTA
      PFIN(2) = PFSP(2)
      PFIN(1) = PFSP(1)
      PFINPL = PFSPPL + 2.0D0*DELTA
      PFINMI = PFSPMI
      M2FIN = PFINPL*PFINMI-PFIN(1)**2-PFIN(2)**2
      PFIN(5) = DSQRT(M2FIN)

C      Boost into the CMS of the final state particles
C      Calling sequence for:
C      SUBROUTINE DT_DALTRA(GA,BGX,BGY,BGZ,PCX,PCY,PCZ,EC,P,PX,PY,PZ,E)
C                             gamma,beta      input 4-v     output 4-v
      BGX = -PFSP(1)/PFSP(5)
      BGY = -PFSP(2)/PFSP(5)
      BGZ = -PFSP(3)/PFSP(5)
      GA = PFSP(4)/PFSP(5)

      NPARTS = 0
      DO IP=1,NMXHKK
         IF (ISTHKK(IP).EQ.1) THEN
            NPARTS = NPARTS + 1
            IDXHKK(NPARTS) = IP
            CALL DT_DALTRA(GA,BGX,BGY,BGZ,
     &           PHKK(1,IP),PHKK(2,IP),PHKK(3,IP),PHKK(4,IP),
     &           PTMPTOT,PTMP(NPARTS,1),PTMP(NPARTS,2),
     &           PTMP(NPARTS,3),PTMP(NPARTS,4))
            PTMP(NPARTS,5) = PHKK(5,IP)
         ENDIF
      ENDDO

C     Scale all particle 3-momenta by a factor ALPHA to reach desired s.
C     This is similar to what happens in pfshift.f
C 

      CALL SCLSUBEVT(M2FSP,M2FIN,EPS,VERB,NPARTS,PTMP,W2OUT,
     &     NSCLTR,SCLFAIL)

      IF (SCLFAIL) THEN
         WRITE(*,*)'LC-based rescale failed after ',NSCLTR,' tries.'
         WRITE(*,*)'HCMS. Before rescale attempt:'
         DO IP=1,NPARTS
            WRITE(*,*)IP," ",PHKK(1,IDXHKK(IP))," ",
     &           PHKK(2,IDXHKK(IP))," ",PHKK(3,IDXHKK(IP))," ",
     &           PHKK(4,IDXHKK(IP))," ",PHKK(5,IDXHKK(IP))
         ENDDO
         WRITE(*,*)'HCMS. After LC-based rescale attempt'
         DO IP=1,NPARTS
            WRITE(*,*)IP," ",PTMP(IP,1)," ",PTMP(IP,2)," ",
     &           PTMP(IP,3)," ",PTMP(IP,4)," ",PTMP(IP,5)
         ENDDO
      ELSEIF (VERB) THEN 
         WRITE(*,*)'LC-based rescale succeeded after ',NSCLTR,
     &        ' tries.'
      ENDIF


C     Boost particles to the desired total 4-momentum in the naive HCMS
C

      BGX = PFIN(1)/PFIN(5)
      BGY = PFIN(2)/PFIN(5)
      BGZ = PFIN(3)/PFIN(5)
      GA = PFIN(4)/PFIN(5)

      DO IP=1,NPARTS
         IDX = IDXHKK(IP)
         CALL DT_DALTRA(GA,BGX,BGY,BGZ,
     &        PTMP(IP,1),PTMP(IP,2),PTMP(IP,3),PTMP(IP,4),
     &        PTMPTOT,PHKK(1,IDX),PHKK(2,IDX),PHKK(3,IDX),PHKK(4,IDX))

      ENDDO

      RETURN
      END
