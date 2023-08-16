      SUBROUTINE FIXESTAR(ESTI,ESTF,AMRCL0,PRCL,VERB,IREJ)
C*********************************************************************
 
C... FIXESTAR - Mark D. Baker 2023-08-05
C...
C... Modify 4-momenta of all particles except e' (status 99 here)
C... To change nuclear excitation energy E* from ESTI to ESTF.
C... If it is not possible, ESTF might need to be modified.
C...
C... ESTI: input E*
C... ESTF: input desired final E*.  Returns actual final E*
C... AMRCL0: ground state mass of the nuclear remnant
C... PRCL: 4-momentm of projectile(1) and target(2) remnant 
C... VERB: verbosity flag
C... IREJ: 0 = success. 1=error - skip event. 2=changed ESTF.
C...
C... We should be working in the naive HCMS of the e+N collision
C... with gamma* momentum along +z and incoming ion along -z  
C... 
C... Change the P+ and M of the remnant to fix E*, while leaving P-,Pt alone
C... For FSP (Final State Particles excluding e' & remnant) , change total
C... P+ and s to conserve overall P+ while leaving total P-,Pt alone
C...
C... The change of s of the FSPs assumes at least 2 particles. If there is 
C... only one final state particle, use DT_MASHEL instead. It allows you 
C... to put 2 particles on mass-shell and 2 is all we've got in this case.
C...
C... If E* (ESTF) is too large, it is not possible to leave all particles
C... on mass shell. In this case, modify ESTF to be between 
C... 0 and 0.99*ESTARMAX
C...
      IMPLICIT NONE
C     Calling params.: in    I/O    in         I/O
      DOUBLE PRECISION ESTI, ESTF, AMRCL0(2), PRCL(2,4)
C              in
      LOGICAL VERB
C              out
      INTEGER IREJ
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
      DOUBLE PRECISION PFIN(5),M2FIN,DELTA,PFINPL,PFINMI,PIN(5)
      DOUBLE PRECISION PTMP(MAXPRTS,5),PTMPTOT, MSUM, ESTARMAX
      DOUBLE PRECISION EPS, W2OUT, GA, BGX, BGY, BGZ, DELTAMAX
      PARAMETER (EPS=1.0D-9)
      INTEGER IP, MU, IDX

      EXTERNAL GETESTAR
      DOUBLE PRECISION GETESTAR

      IREJ=0
      MSUM=0.0D0
      DO MU=1,4
         PFSP(MU)=0.0D0
      ENDDO
      
      NPARTS=0
      DO IP=1,NHKK
C         WRITE(*,*) IP, ISTHKK(IP), IDHKK(IP), PHKK(1,IP), PHKK(2,IP),
C     &        PHKK(3,IP), PHKK(4,IP), PHKK(5,IP)
         IF (ISTHKK(IP).EQ.1) THEN
            NPARTS = NPARTS + 1
            IDXHKK(NPARTS) = IP
            DO MU=1,4
               PFSP(MU) = PFSP(MU) + PHKK(MU,IP)
            ENDDO
            MSUM = MSUM + PHKK(5,IP)
         ENDIF
      ENDDO

      IF (NPARTS.EQ.1) THEN

C     Use PIN(1-4) to hold the incoming cluster: PRCL(2,*), 
C     PFIN(1-5) is the final remnamt

         DO MU=1,4
            PIN(MU) = PRCL(2,MU)
         ENDDO
         PFIN(5)=AMRCL0(2)+ESTF
         CALL DT_MASHEL(PIN,PFSP,PFIN(5),PHKK(5,IDXHKK(1)),
     &                  PFIN,PHKK(1,IDXHKK(1)),IREJ)
         IF (VERB) THEN
            WRITE(*,*) "FIXESTAR with 1 fsp"
            WRITE(*,*) "PIN(1-4),PFSP(1-4),PFIN(5),MFSP,PFIN,PFSPout",
     &        PIN(1),PIN(2),PIN(3),PIN(4),PFSP(1),PFSP(2),PFSP(3),
     &        PFSP(4),PFIN(5),PHKK(5,IDXHKK(1)),PFIN(1),PFIN(2),
     &        PFIN(3),PFIN(4),PHKK(1,IDXHKK(1)),PHKK(2,IDXHKK(1)),
     &        PHKK(3,IDXHKK(1)),PHKK(4,IDXHKK(1))
         ENDIF
         IF (IREJ.NE.0) THEN
            WRITE(*,*) ('FIXESTAR ERROR. DT_MASHEL FAILED.')
            IREJ = 1
            GOTO 9999
         ENDIF
         DO MU=1,4
            PRCL(2,MU) = PFIN(MU)
         ENDDO
         
      ELSEIF (NPARTS.GE.2) THEN

         PFSPPL = PFSP(4) + PFSP(3)
         PFSPMI = PFSP(4) - PFSP(3)
         M2FSP = PFSPPL*PFSPMI-PFSP(1)**2-PFSP(2)**2
         IF (M2FSP.LT.0) STOP ('FIXESTAR FATAL. M2FSP<0.')
         PFSP(5) = DSQRT(M2FSP)

C         WRITE(*,*) 'PFSP(1-5):',PFSP(1),PFSP(2),PFSP(3),PFSP(4),PFSP(5)
C         WRITE(*,*) 'PFSPPL,PFSPMI:',PFSPPL,PFSPMI
         PRMNPL = PRCL(2,4) + PRCL(2,3)
         PRMNMI = PRCL(2,4) - PRCL(2,3)
         MRMN = DSQRT(PRMNPL*PRMNMI-PRCL(2,2)**2-PRCL(2,1)**2)
         WRITE(*,*) 'MRMN, ESTI+AMRCL0(2)',MRMN, ESTI+AMRCL0(2)
         IF (ABS(MRMN-ESTI-AMRCL0(2)).GT.1D-07) STOP ('FIXESTAR FATAL')

         DELTAMAX = 0.5D0*(M2FSP-MSUM*MSUM)/PFSPMI
         ESTARMAX = DSQRT(MRMN*MRMN+2.0D0*PRMNMI*DELTAMAX)-AMRCL0(2)

         IF (ESTARMAX.LE.0.0D0) THEN
            IREJ = 1
C This should be rare enough to print every time
            WRITE (*,*) 'FIXESTAR ERROR: ESTARMAX<0.'
            WRITE (*,*) 'Event is below production threshhold'
            GOTO 9999
         ELSEIF (ESTF.GE.0.999*ESTARMAX) THEN
C This should be rare enough to print every time
            WRITE (*,*) 'FIXESTAR WARNING: ESTAR>ESTARMAX.'
            WRITE (*,*) 'ESTF,ESTARMAX,MRMN,M0,PRMNMI,M2FSP,MSUM:'
            WRITE (*,*) ESTF,ESTARMAX,MRMN,AMRCL0,PRMNMI,M2FSP,MSUM
            ESTF = GETESTAR(2,0.001*ESTARMAX,0.999*ESTARMAX)
            WRITE (*,*) 'Rerolling ESTF: ',ESTF
            IREJ = 2
         ENDIF

C     Required change in P+remn with fixed P-remn,Ptremn for ESTF
         DELTA = 0.5D0*(ESTF+ESTI+2.0D0*AMRCL0(2))*(ESTF-ESTI)/PRMNMI
         PRCL(2,4) = PRCL(2,4)+DELTA
         PRCL(2,3) = PRCL(2,3)+DELTA
         PFIN(4) = PFSP(4)-DELTA
         PFIN(3) = PFSP(3)-DELTA
         PFIN(2) = PFSP(2)
         PFIN(1) = PFSP(1)
         PFINPL = PFSPPL - 2.0D0*DELTA
         PFINMI = PFSPMI
         M2FIN = PFINPL*PFINMI-PFIN(1)**2-PFIN(2)**2
         IF (M2FIN.LE.MSUM*MSUM) THEN
            WRITE(*,*) 'FATAL ERROR: FIXESTAR'
            WRITE(*,*) 'Cannot put everything on mass shell'
            WRITE(*,*) 'This should not happen.'
            IREJ=1
            GOTO 9999
         ENDIF
         PFIN(5) = DSQRT(M2FIN)

C      Boost into the CMS of the final state particles
C      Calling sequence for:
C      SUBROUTINE DT_DALTRA(GA,BGX,BGY,BGZ,PCX,PCY,PCZ,EC,P,PX,PY,PZ,E)
C                             gamma,beta      input 4-v     output 4-v
         BGX = -PFSP(1)/PFSP(5)
         BGY = -PFSP(2)/PFSP(5)
         BGZ = -PFSP(3)/PFSP(5)
         GA = PFSP(4)/PFSP(5)
         
         DO IP=1,NPARTS
            IDX = IDXHKK(IP)
            CALL DT_DALTRA(GA,BGX,BGY,BGZ,
     &           PHKK(1,IDX),PHKK(2,IDX),PHKK(3,IDX),PHKK(4,IDX),
     &           PTMPTOT,PTMP(IP,1),PTMP(IP,2),
     &           PTMP(IP,3),PTMP(IP,4))
            PTMP(IP,5) = PHKK(5,IDX)
         ENDDO

C     Scale all particle 3-momenta by a factor ALPHA to reach desired s.
C     This is similar to what happens in pfshift.f
C 

         CALL SCLSUBEVT(M2FSP,M2FIN,EPS,VERB,NPARTS,PTMP,W2OUT,
     &        NSCLTR,SCLFAIL)

         IF (SCLFAIL) THEN
            WRITE(*,*)'LC-based rescale failed after ',NSCLTR,' tries.'
            WRITE(*,*)'HCMS. Before rescale attempt:'
            DO IP=1,NPARTS
               WRITE(*,*)IP," ",PHKK(1,IDXHKK(IP))," ",
     &              PHKK(2,IDXHKK(IP))," ",PHKK(3,IDXHKK(IP))," ",
     &              PHKK(4,IDXHKK(IP))," ",PHKK(5,IDXHKK(IP))
            ENDDO
            WRITE(*,*)'HCMS. After LC-based rescale attempt'
            DO IP=1,NPARTS
               WRITE(*,*)IP," ",PTMP(IP,1)," ",PTMP(IP,2)," ",
     &              PTMP(IP,3)," ",PTMP(IP,4)," ",PTMP(IP,5)
            ENDDO
         ELSEIF (VERB) THEN 
            WRITE(*,*)'LC-based rescale succeeded after ',NSCLTR,
     &           ' tries.'
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
     &           PTMP(IP,1),PTMP(IP,2),PTMP(IP,3),PTMP(IP,4),
     &           PTMPTOT,PHKK(1,IDX),PHKK(2,IDX),PHKK(3,IDX),PHKK(4,IDX)
     &           )
         ENDDO

      ELSE

         STOP ('FATAL ERROR IN FIXESTAR: NO FINAL STATE PARTICLES!')

      ENDIF
         
 9999 RETURN
      END
