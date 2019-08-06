      SUBROUTINE DEUTFIX(NU,Q2,MDEUT)
C
C     2018-08-25 Mark D. Baker - Initial Version
C
C     Explicit Input: NU, Q2 - event kinematics
C                     MDEUT - Deuteron Mass
C     Implicit I/O:   An event in /PYJETS/ in the ion rest frame
C                     INCLUDING the spectator nucleon.
C
C     This subroutine compensates for a flaw in the handling of 
C     deuteron beams in e+D collisions. BeAGLE (DPMJET) treats a 
C     deuteron as a pair of on mass-shell nucleons with momentum 
C     pF and -pF (3-vector) in the ion rest frame. But this 
C     violates energy conservation in this frame (4-momentum in 
C     general) since the n+p already have more energy than the 
C     deuteron (which is bound!) and any relative momentum between 
C     p + n only compounds the problem.
C
C     For larger nuclei, there is a mean field potential which reduces
C     the energy of outgoing nucleons as well as a nuclear remnant to
C     absorb any remaining 4-momentum imbalance. 
C
C     This routine should be called in the ion rest frame, and currently
C     only for deuterons. It works better if called before hadronization
C     since the iterative procedure converges better for a smaller # of
C     particles.
C
C     Considering the "hadronic" (gamma*+d) part of the event, without 
C     the scattered electron, we would expect the four-momentum in the
C     ion rest frame to be: Wfull = {0, 0, sqrt(nu^2+Q^2); (nu+m_d)}.
C     Due to the incorrect input kinematics, we have too much energy
C     (nu+Ep+En), but the correct 3-momentum. 
C
C     The approach is:
C     Step 1 is to boost into the rest frame of the hadronic event
C
C     Step 2 is to scale all momenta by a common factor
C     ASCALE so that the W2full is correct for gamma* + d,
C     namely W2full = 2 M_d nu + M_d^2 - Q^2
C     In this frame we start with total 4-momentum: (0,0,0;W_oops)
C     and scale to (0,0,0;Wfull) (similar to S/R PFSHIFT).
C
C     For 3 particles the exact formula is really complicated.
C     For >3 particles I don't believe there is a closed form solution.
C     Therefore we'll use an iterative procedure which assumes that
C     ASCALE ~ 1 and converges to the correct value after a couple of
C     iterations.
C
C     The idea is to calculate "alpha" (ASCALE), the factor by which all 
C     momenta should be multiplied by in the gamma*+D HCMS in order to go 
C     from W2oops -> W2F
C
C     We use an approximate formula which just keeps the terms to 
C     O(delta) where delta = alpha - 1 is assumed small.
C
C     Recall that Woops = Sum (E_i) = Sum sqrt(p^2_i + m^2_i)
C     We want WF ~ Sum sqrt[(1+2delta)p^2_i + m^2_i]
C                ~ Sum E_i (1 + delta p^2_i/E^2_i)
C     WF ~ Woops + delta Sum (p^2/E)_i 
C 
C Approximate scale factor ASCALE=1+delta where
C delta=(WF-Woops)/SUM(p^2/E)
C
C     Step 3 is to boost the event with Wfull=(0,0,0;Wfull) back to 
C     a frame with 3-momentum = qz = sqrt(nu^2+Q^2). The ion rest frame,
C     but now with the correct kinematics.
C
C NOTE: The event history is not corrected in any way, but represents
C       the original interaction kinematics.
C
      IMPLICIT NONE
      DOUBLE PRECISION NU, Q2, MDEUT

      include 'beagle.inc'
C      include "py6strf.inc"   ! Temporary! Just use for debug output

C      include 'pythia.inc' - conflicts with IMPLICIT NONE
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.0D0)

CCc...added by liang & Mark to include pythia energy loss datas
C      double precision PAUX, DPF
C      COMMON /PFAUX/ PAUX(4), DPF(4)

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 
      integer NEVENT, ICASCA

C Local
      DOUBLE PRECISION EPSPF
      PARAMETER (EPSPF=1.0D-9)
      INTEGER MAXTRY, MAXPRTS,NDIM
      PARAMETER (MAXTRY=10)
      PARAMETER (MAXPRTS=20)
      PARAMETER (NDIM=4)
      DOUBLE PRECISION W2F, W2TRY(0:MAXTRY), PSUM(NDIM), W2RAW, WRAW
      DOUBLE PRECISION ASCALE(MAXTRY), ASCLFL
      DOUBLE PRECISION S2SUM 
      INTEGER NPRTNS,NLSCAT,IDIM,NSCLTR, ITRK
      INTEGER INDXP(MAXPRTS)
      LOGICAL W2FAIL

      INTEGER INDEX
      DOUBLE PRECISION BETAZ

      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         WRITE (*,*) 'DEUTFIX: About to fix e+D kinematics'
         CALL PYLIST(2)
         WRITE (*,*)
      ENDIF

C     Identify the stable particles and assemble W^mu_oops (PSUM)
      NLSCAT = 0
      NPRTNS = 0
      DO IDIM=1,NDIM
         PSUM(IDIM)=ZERO       ! sum p^mu for all stable except e'
      ENDDO
      DO ITRK=1,N
         IF(K(ITRK,1).EQ.1 .OR. K(ITRK,1).EQ.2) THEN
            IF ( (ABS(K(ITRK,2)).EQ.11 .OR. ABS(K(ITRK,2)).EQ.13) .AND.
     &           K(ITRK,3).EQ.3) THEN
               NLSCAT = NLSCAT+1
            ELSE
               NPRTNS=NPRTNS+1
               IF (NPRTNS.GT.MAXPRTS) 
     &              STOP('DEUTFIX: FATAL ERROR. Too many partons')
               INDXP(NPRTNS)=ITRK
               DO IDIM=1,NDIM
                     PSUM(IDIM)=PSUM(IDIM)+P(ITRK,IDIM)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      IF (NLSCAT.NE.1) 
     &     STOP "ERROR! BAD EVENT CONFIG. Scattered leptons .ne. 1"
      IF (NPRTNS.LT.2)
     &     STOP "ERROR! BAD EVENT CONFIG. Fewer than two particles"

      W2RAW = (PSUM(4)-PSUM(3))*(PSUM(4)+PSUM(3))-PSUM(1)**2-PSUM(2)**2
      WRAW = SQRT(W2RAW)
      W2F = 2.0D0*MDEUT*NU + MDEUT*MDEUT - Q2
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         WRITE(*,*) 'W^mu_full (correct): {0, 0, ', SQRT(NU*NU+Q2),';',
     &        NU+MDEUT,'}'
         WRITE(*,*) 'W^mu_oops found: {,',PSUM(1),', ',PSUM(2),', ',
     &        PSUM(3),'; ',PSUM(4),'}'
         WRITE(*,*) 'W2_full, W2_oops: ',W2F,', ', W2RAW
      ENDIF

C     Step 1: Boost into hadronic rest frame and calculate S2SUM,PSUM
      BETAZ = -PSUM(3)/PSUM(4)
      S2SUM = ZERO
      DO IDIM=1,NDIM
         PSUM(IDIM)=ZERO        ! sum p^mu for all stable except e'
      ENDDO
      DO ITRK=1,NPRTNS
         INDEX = INDXP(ITRK)
         CALL PYROBO(INDEX, INDEX, ZERO, ZERO, ZERO, ZERO, BETAZ)
         S2SUM = S2SUM + 
     &        (P(INDEX,4)-P(INDEX,5))*(ONE+P(INDEX,5)/P(INDEX,4))
         DO IDIM=1,NDIM
            PSUM(IDIM)=PSUM(IDIM)+P(INDEX,IDIM)
         ENDDO
      ENDDO

      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         WRITE(*,*) '(Sum 3-momentum)^2 should be zero: ',
     &        (PSUM(1)*PSUM(1)+PSUM(2)*PSUM(2)+PSUM(3)*PSUM(3))
      ENDIF
      IF ((PSUM(1)*PSUM(1)+PSUM(2)*PSUM(2)+PSUM(3)*PSUM(3)).GT.
     &     EPSPF/10.0D0)
     &        WRITE(*,*) 'DEUTFIX ERROR: TOTAL 3-MOMENTUM NOT ZERO!.'

C     Step 2: Iteratively scale the particle 3-momenta until we reach the 
C     correct W value for the gamma*+D reaction products. 
      NSCLTR=0
      W2TRY(NSCLTR) = PSUM(4)*PSUM(4)
      DO WHILE (NSCLTR.LT.MAXTRY .AND.
     &     ABS(W2TRY(NSCLTR)/W2F-ONE).GT.EPSPF)
         NSCLTR=NSCLTR+1
         ASCALE(NSCLTR) = ONE+(SQRT(W2F)-SQRT(W2TRY(NSCLTR-1)))/S2SUM
         IF ( (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) .OR.
     &        (IOULEV(4).GE.1 .AND. NSCLTR.GT.MAXTRY-3) ) then
            WRITE(*,*)'W2 inaccurate. Iteration # ',NSCLTR
            WRITE(*,*)'Next ASCALE factor = ',ASCALE(NSCLTR)
         ENDIF
C     Zero out our sums. 3-momentum sum should be zero now.
         PSUM(4)=ZERO
         S2SUM=ZERO
         DO ITRK=1,NPRTNS
            INDEX = INDXP(ITRK)
            DO IDIM=1,3
               P(INDEX,IDIM)=ASCALE(NSCLTR)*P(INDEX,IDIM)
            ENDDO
            P(INDEX,4)= SQRT( P(INDEX,5)**2+
     &           (P(INDEX,1)**2+P(INDEX,2)**2+P(INDEX,3)**2))
            PSUM(4) = PSUM(4) + P(INDEX,4)
            S2SUM = S2SUM + (P(INDEX,4)-P(INDEX,5))
     &           *(ONE+P(INDEX,5)/P(INDEX,4))
         ENDDO
         W2TRY(NSCLTR) = PSUM(4)*PSUM(4)
      
         IF ( (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) .OR.
     &        (IOULEV(4).GE.1 .AND. NSCLTR.GT.MAXTRY-3) ) then
            WRITE(*,*)"W2F:",W2F,"Iteration ",NSCLTR," Scaled W2:",
     &           W2TRY(NSCLTR),"W2try(",NSCLTR,")/W2F",
     &           W2TRY(NSCLTR)/W2F,"NPRTNS:",NPRTNS
            WRITE(*,*)'HCMS. After W2 rescale iteration ', NSCLTR
            DO ITRK=1,NPRTNS
               INDEX = INDXP(ITRK)
               WRITE(*,*)ITRK," ",P(INDEX,1)," ",P(INDEX,2)," ",
     &              P(INDEX,3)," ",P(INDEX,4)," ",P(INDEX,5)
            ENDDO
            WRITE(*,*)"Leads to ..."
            WRITE(*,*)"PSUM(1-4): 0, 0, 0, ",PSUM(4)
         ENDIF
      ENDDO

      W2FAIL = (ABS(W2TRY(NSCLTR)/W2F-1).GT.EPSPF)

      IF (.NOT. W2FAIL) THEN
         IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
            WRITE(*,*)'Deuteron W2-rescale succeeded after ',NSCLTR,
     &           ' tries.'
         ENDIF
      ELSE
         WRITE(*,*)'Deuteron W2-rescale failed after ',NSCLTR,' tries.'
         WRITE(*,*)'g*+D-HCMS. After Deuteron W2 rescale attempt'
         DO ITRK=1,NPRTNS
            WRITE(*,*)ITRK," ",P(INDXP(ITRK),1)," ",P(INDXP(ITRK),2),
     &           " ",P(INDXP(ITRK),3)," ",P(INDXP(ITRK),4)," ",
     &           P(INDXP(ITRK),5)
         ENDDO
      ENDIF
      
      IF ( IOULEV(4).GE.2 .AND. 
     &     (NEVENT.LE.IOULEV(5) .OR. NSCLTR.GT.MAXTRY-3) ) THEN
         WRITE(*,*)
         WRITE(*,*)'Iteration   W2/W2F'
         WRITE(*,*)'          0   ',W2RAW/W2F
         ASCLFL = ONE
         DO ITRK=1,NSCLTR
            WRITE(*,*) ITRK, "  ", W2TRY(ITRK)/W2F
            ASCLFL = ASCLFL * ASCALE(ITRK)
         ENDDO
         IF (W2FAIL) THEN
            WRITE(*,*)'Failed to converge to level of: ',EPSPF
         ELSE
            WRITE(*,*)'Success. Converged to level of: ',EPSPF
         ENDIF
         WRITE(*,*)'WF/WRAW:     ',SQRT(W2F)/WRAW
         WRITE(*,*)'ASCALE(1):   ',ASCALE(1)
         WRITE(*,*)'ASCALE(full):',ASCLFL
      ENDIF

      IF (USERSET.EQ.6) THEN
         USER1=W2F
         USER2=W2TRY(NSCLTR)/W2F - ONE
C         USER3=DBLE(NPRTNS)
         IF (W2FAIL) THEN
            USER3 = -1.0D0
         ELSE
            USER3=DBLE(NSCLTR)
         ENDIF
      ENDIF

C     Step 3: Boost back into the ion rest frame and calculate PSUM
      BETAZ = ONE/SQRT(ONE + W2TRY(NSCLTR)/(Q2+NU*NU))
      DO IDIM=1,NDIM
         PSUM(IDIM)=ZERO        ! sum p^mu for all stable except e'
      ENDDO
      DO ITRK=1,NPRTNS
         INDEX = INDXP(ITRK)
         CALL PYROBO(INDEX, INDEX, ZERO, ZERO, ZERO, ZERO, BETAZ)
         DO IDIM=1,NDIM
            PSUM(IDIM)=PSUM(IDIM)+P(INDEX,IDIM)
         ENDDO
      ENDDO

      W2RAW = (PSUM(4)-PSUM(3))*(PSUM(4)+PSUM(3))+PSUM(1)**2+PSUM(2)**2
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         WRITE(*,*) 'After Deuteron kinematic fix:'
         WRITE(*,*) 'W^mu_full (correct): {0, 0, ', SQRT(NU*NU+Q2),';',
     &        NU+MDEUT,'}'
         WRITE(*,*) 'W^mu  found: {,',PSUM(1),', ',PSUM(2),', ',
     &        PSUM(3),'; ',PSUM(4),'}'
         WRITE(*,*) 'W2_full, W2_fixed: ',W2F,', ', W2RAW
      ENDIF

      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*)"PYLIST: After deuteron fix"
         CALL PYLIST(2)
      ENDIF
      
      RETURN
      END
