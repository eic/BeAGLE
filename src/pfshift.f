      SUBROUTINE PFSHIFT(WRAW,W2RAW,MNUCL)
C
C     2018-07-13 Mark D. Baker - Initial Version
C
C     This subroutine adds the struck nucleon pF back to the
C     Pythia subevent partonic skeleton after the fact.
C     We start and end in the TRF, but operate in the naive HCMS
C     defined by gamma* + N at rest in nuclear TRF.
C
C     Step 1 is to scale all momenta in the HCMS by a common factor
C     ASCALE so that the W2 is correct for gamma* + moving N rather
C     than gamma* + N at rest in TRF.
C     For 2 particles there is a reasonable exact formula.
C     For 3 particles the exact formula is really complicated.
C     For >3 particles I don't believe there is a closed form solution.
C     Therefore we'll use an iterative procedure which assumes that
C     ASCALE ~ 1. For N=2, the 1st step of the procedure uses the exact
C     formula and should therefore converge immediately. 
C
C     Step 2 is to boost all momenta so the that subevent has the
C     correct momentum in the HCMS (and ultimately TRF)
C 
C     WRAW, W2RAW, MNUCL are input parameters only - from VINT(1),(2),(4)
C
      IMPLICIT NONE
      DOUBLE PRECISION WRAW, W2RAW, MNUCL

      include 'beagle.inc'
C      include "py6strf.inc"   ! Temporary! Just use for debug output

C      include 'pythia.inc' - conflicts with IMPLICIT NONE
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V

* Lorentz-parameters of the current interaction from DPMJET
      COMMON /DTLTRA/ GACMS(2),BGCMS(2),GALAB,BGLAB,BLAB,
     &                UMO,PPCM,EPROJ,PPROJ
      DOUBLE PRECISION GACMS,BGCMS,GALAB,BGLAB,BLAB
      DOUBLE PRECISION UMO, PPCM, EPROJ, PPROJ

Cc...added by liang & Mark to include pythia energy loss datas
      double precision PAUX, DPF
      COMMON /PFAUX/ PAUX(4), DPF(4)

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 
      integer NEVENT, ICASCA

C Local
      DOUBLE PRECISION EPSPF
      PARAMETER (EPSPF=1.0D-9)
      INTEGER MAXTRY, MAXPRTS,NDIM,NDIMM
      PARAMETER (MAXTRY=10)
      PARAMETER (MAXPRTS=20)
      PARAMETER (NDIM=4)
      PARAMETER (NDIMM=5)
      DOUBLE PRECISION W2F, W2TRY(MAXTRY), PSUM(NDIM) ! PSMTRY(NDIM)
C      DOUBLE PRECISION PSMTR1(NDIM),PSMTR2(NDIM), EETEMP
C      DOUBLE PRECISION PSUM1(NDIM),PSUM2(NDIM), SUMM2, SUMM1
      DOUBLE PRECISION ASCALE(MAXTRY), PHIGH2, SQRM1, SQRM2, ASCLFL
      DOUBLE PRECISION ASCTMP1 ! ASCTMP2, ASCTMP3
      DOUBLE PRECISION W2TTMP1 ! W2TTMP2, W2TTMP3 
      DOUBLE PRECISION S2SUM ! DELTAA, EPSIL1, S4SUM, BETA2
      DOUBLE PRECISION FERBX, FERBY, FERBZ, FERGAM
      INTEGER NPRTNS,NLSCAT,IDIM,NSCLTR, ITRK
      DOUBLE PRECISION PTMP(MAXPRTS,NDIMM)
      INTEGER INDXP(MAXPRTS)
      LOGICAL W2FAIL

      IF (WRAW.LT.1.0 .OR. ABS(W2RAW-WRAW*WRAW).GT.0.001 .OR.
     & MNUCL.LT.0.9 .OR. MNUCL.GT.1.0) THEN
         WRITE(*,*)'WRAW,W2RAW,MNUCL: ',WRAW,W2RAW,MNUCL
         STOP 'PFSHIFT: FATAL ERROR. BAD KINEMATICS!'
      ENDIF
C     Boost into naive HCMS  (assumes nucleon at rest in A-TRF)
      CALL PYROBO(0,0,0.0D0,0.0D0,0.0D0,0.0D0,-BGCMS(2)/GACMS(2))

      W2F = W2RAW + 2.*WRAW*DPF(4) - 2.*MNUCL*EKF
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         write(*,*) "W2 ignore pF:", W2RAW
         write(*,*) "W2 corrected: ", W2F
         write(*,*) "gamma*beta, gamma, beta:",BGCMS(2),GACMS(2),
     &        BGCMS(2)/GACMS(2)
      ENDIF
      
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         write(*,*) "DT_PYEVNTEP: HCMS g*=z, e' px>0 py=0"
         CALL PYLIST(2)
      ENDIF

C     Analyze the events and extract the "partonic" skeleton
      NPRTNS=0
      NLSCAT=0
      S2SUM=0.0D0
C      S4SUM=0.0D0
      DO IDIM=1,NDIM
         PSUM(IDIM)=0.0D0       ! sum p^mu for all stable except e'
C         PSUM1(IDIM)=0.0D0      ! sum p^mu for pz<0: P1
C         PSUM2(IDIM)=0.0D0      ! sum p^mu for pz>=0: P2
c         PSMTRY(IDIM)=0.0D0     ! as above, but after momentum scaling
C         PSMTR1(IDIM)=0.0D0
C         PSMTR2(IDIM)=0.0D0
      ENDDO
C      SUMM1=0.0D0               ! sum m^2 for pz<0 - along P1
C      SUMM2=0.0D0               ! sum m^2 for pz>0 - along P2
      DO ITRK=1,N
         IF(K(ITRK,1).EQ.1 .OR. K(ITRK,1).EQ.2) THEN
            IF ( (ABS(K(ITRK,2)).EQ.11 .OR. ABS(K(ITRK,2)).EQ.13) .AND.
     &           K(ITRK,3).EQ.3) THEN
               NLSCAT = NLSCAT+1
            ELSE
               NPRTNS=NPRTNS+1
               IF (NPRTNS.GT.MAXPRTS) 
     &              STOP('PFSHIFT: FATAL ERROR. Too many partons')
               INDXP(NPRTNS)=ITRK
               DO IDIM=1,NDIMM
                  PTMP(NPRTNS,IDIM)=P(ITRK,IDIM)
                  IF (IDIM.LE.NDIM) THEN
                     PSUM(IDIM)=PSUM(IDIM)+PTMP(NPRTNS,IDIM)
                  ENDIF
               ENDDO
               S2SUM = S2SUM + (PTMP(NPRTNS,4)-PTMP(NPRTNS,5))
     &              *(1.0D0+PTMP(NPRTNS,5)/PTMP(NPRTNS,4))
            ENDIF
         ENDIF
      ENDDO
      IF (NLSCAT.NE.1) 
     &     STOP "ERROR! BAD EVENT CONFIG. Scattered leptons .ne. 1"
         
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         W2TRY(1) = PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2
         WRITE(*,*)"W2 from Pythia:",W2RAW,"W2 calc.:",W2TRY(1)
         WRITE(*,*)"PSUM(1-4):",PSUM(1),PSUM(2),PSUM(3),PSUM(4)
      ENDIF
C
C Calculate "alpha", the factor by which all momenta should be 
C multiplied by in the HCMS in order to go from W2RAW -> W2F
C
C We use an exact formula in the case where there are only two
C particles.
C 
C For N>2, we use an approximate formula which just keeps the
C terms to O(delta) where delta = alpha - 1 is assumed small.
C
C Approximate scale factor ASCALE=1+delta where
C delta=(WF-W0)/SUM(p^2/E)
C
      IF (NPRTNS.GT.2) THEN
         ASCALE(1) = 1.0D0 + (SQRT(W2F)-WRAW)/S2SUM
      ELSEIF (NPRTNS.EQ.2) THEN
         PHIGH2 =( PTMP(1,1)**2+PTMP(1,2)**2+PTMP(1,3)**2 +
     &             PTMP(2,1)**2+PTMP(2,2)**2+PTMP(2,3)**2 )/2.0D0
         SQRM1 = PTMP(1,5)*PTMP(1,5)
         SQRM2 = PTMP(2,5)*PTMP(2,5)
         ASCALE(1) = (W2F-2.0D0*(SQRM1+SQRM2)+((SQRM1-SQRM2)**2)/W2F )/
     &        (4.0D0*PHIGH2)
         ASCALE(1)=SQRT(ASCALE(1))
      ELSE
         STOP ('PFSHIFT: FATAL ERROR. NPRTNS<2!')
      ENDIF
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         WRITE(*,*)'HCMS. Before pF-based W2 rescale. p5:'
         DO ITRK=1,NPRTNS
            WRITE(*,*)ITRK," ",PTMP(ITRK,1)," ",PTMP(ITRK,2)," ",
     &           PTMP(ITRK,3)," ",PTMP(ITRK,4)," ",PTMP(ITRK,5)
         ENDDO
         WRITE(*,*)"W2F,PHIGH2,SQRM1,SQRM2,NUMER,DENOM:",
     &        W2F,PHIGH2,SQRM1,SQRM2,
     &        (W2F*W2F-2.0D0*W2F*(SQRM1+SQRM2)+(SQRM1-SQRM2)**2),
     &        (4.0D0*PHIGH2*W2F)
         WRITE(*,*)"ASCALE(1):",ASCALE(1)
      ENDIF

      S2SUM = 0.0D0
      DO IDIM=1,NDIM
         PSUM(IDIM)=0.0D0
C         PSUM1(IDIM)=0.0D0
C         PSUM2(IDIM)=0.0D0
      ENDDO
      DO ITRK=1,NPRTNS
         DO IDIM=1,NDIM
            IF (IDIM.LE.3) THEN
               PTMP(ITRK,IDIM) = ASCALE(1)*PTMP(ITRK,IDIM)
               PSUM(IDIM)=PSUM(IDIM)+PTMP(ITRK,IDIM)
            ELSE
C     Note: P already scaled. Just recalc E.
               PTMP(ITRK,4)= SQRT( PTMP(ITRK,5)**2+
     &              (PTMP(ITRK,1)**2+PTMP(ITRK,2)**2+PTMP(ITRK,3)**2))
               PSUM(4) = PSUM(4) + PTMP(ITRK,4)
            ENDIF
         ENDDO
         S2SUM = S2SUM + (PTMP(ITRK,4)-PTMP(ITRK,5))
     &        *(1.0D0+PTMP(ITRK,5)/PTMP(ITRK,4))
      ENDDO
      W2TRY(1)  = (PSUM(4)-PSUM(3))*(PSUM(4)+PSUM(3))
     &     -PSUM(1)**2-PSUM(2)**2
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         WRITE(*,*)"WF/WRAW, ASCALE(1):", SQRT(W2F/W2RAW),ASCALE(1)
         WRITE(*,*)"W2F:",W2F,"Scaled W2:",W2TRY(1),"W2try/W2F",
     &        W2TRY(1)/W2F,"NPRTNS:",NPRTNS
         WRITE(*,*)'HCMS. After pF-based W2 rescale'
         DO ITRK=1,NPRTNS
            WRITE(*,*)ITRK," ",PTMP(ITRK,1)," ",PTMP(ITRK,2)," ",
     &           PTMP(ITRK,3)," ",PTMP(ITRK,4)," ",PTMP(ITRK,5)
         ENDDO
         WRITE(*,*)"PSUM(1-4):",PSUM(1),PSUM(2),PSUM(3),PSUM(4)
      ENDIF
C     2nd and subsequent iterations
      NSCLTR=1
      DO WHILE (NSCLTR.LT.MAXTRY .AND.
     &     ABS(W2TRY(NSCLTR)/W2F-1.0D0).GT.EPSPF)
         NSCLTR=NSCLTR+1
         PHIGH2 = ASCALE(NSCLTR-1)*ASCALE(NSCLTR-1)*PHIGH2
         ASCALE(NSCLTR) = 1.0D0+(SQRT(W2F)-SQRT(W2TRY(NSCLTR-1)))/S2SUM
         IF ( (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) .OR.
     &        (IOULEV(4).GE.1 .AND. NSCLTR.GT.MAXTRY-3) ) then
            WRITE(*,*)'W2 inaccurate. Iteration # ',NSCLTR
            WRITE(*,*)'Next ASCALE factor = ',ASCALE(NSCLTR)
         ENDIF
         S2SUM=0.0D0
C     3-momenta just scale
         DO IDIM=1,3
            PSUM(IDIM)=ASCALE(NSCLTR)*PSUM(IDIM)
         ENDDO
         PSUM(4)=0.0D0
         DO ITRK=1,NPRTNS
            DO IDIM=1,3
               PTMP(ITRK,IDIM)=ASCALE(NSCLTR)*PTMP(ITRK,IDIM)
            ENDDO
            PTMP(ITRK,4)= SQRT( PTMP(ITRK,5)**2+
     &           (PTMP(ITRK,1)**2+PTMP(ITRK,2)**2+PTMP(ITRK,3)**2))
            PSUM(4) = PSUM(4) + PTMP(ITRK,4)
            S2SUM = S2SUM + (PTMP(ITRK,4)-PTMP(ITRK,5))
     &           *(1.0D0+PTMP(ITRK,5)/PTMP(ITRK,4))
         ENDDO
         W2TRY(NSCLTR) = (PSUM(4)-PSUM(3))*(PSUM(4)+PSUM(3))
     &        -PSUM(1)**2-PSUM(2)**2
      
         IF ( (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) .OR.
     &        (IOULEV(4).GE.1 .AND. NSCLTR.GT.MAXTRY-3) ) then
            WRITE(*,*)"W2F:",W2F,"Iteration ",NSCLTR," Scaled W2:",
     &           W2TRY(NSCLTR),"W2try(",NSCLTR,")/W2F",
     &           W2TRY(NSCLTR)/W2F,"NPRTNS:",NPRTNS
            WRITE(*,*)'HCMS. After W2 rescale iteration ', NSCLTR
            DO ITRK=1,NPRTNS
               WRITE(*,*)ITRK," ",PTMP(ITRK,1)," ",PTMP(ITRK,2)," ",
     &              PTMP(ITRK,3)," ",PTMP(ITRK,4)," ",PTMP(ITRK,5)
            ENDDO
            WRITE(*,*)"Leads to ..."
            WRITE(*,*)"PSUM(1-4):",PSUM(1),PSUM(2),PSUM(3),PSUM(4)
         ENDIF
      ENDDO

      W2FAIL = (ABS(W2TRY(NSCLTR)/W2F-1).GT.EPSPF)

      IF (.NOT. W2FAIL) THEN
         IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
            WRITE(*,*)'pF-based W2-rescale succeeded after ',NSCLTR,
     &           ' tries.'
         ENDIF
      ELSE
         WRITE(*,*)'pF-based W2-rescale failed after ',NSCLTR,' tries.'
         WRITE(*,*)'HCMS. Before rescale attempt:'
         DO ITRK=1,NPRTNS
            WRITE(*,*)ITRK," ",P(INDXP(ITRK),1)," ",P(INDXP(ITRK),2),
     &           " ",P(INDXP(ITRK),3)," ",P(INDXP(ITRK),4)," ",
     &           P(INDXP(ITRK),5)
         ENDDO
         WRITE(*,*)'HCMS. After pF-based W2 rescale attempt'
         DO ITRK=1,NPRTNS
            WRITE(*,*)ITRK," ",PTMP(ITRK,1)," ",PTMP(ITRK,2)," ",
     &           PTMP(ITRK,3)," ",PTMP(ITRK,4)," ",PTMP(ITRK,5)
         ENDDO
      ENDIF
      
      IF ( IOULEV(4).GE.2 .AND. 
     &     (NEVENT.LE.IOULEV(5) .OR. NSCLTR.GT.MAXTRY-3) ) THEN
         WRITE(*,*)
         WRITE(*,*)'Iteration   W2/W2F'
         WRITE(*,*)'          0   ',W2RAW/W2F
         ASCLFL = 1.0D0
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

      IF (USERSET.EQ.3) THEN
         USER1=W2F
         USER2=W2TRY(NSCLTR)/W2F - 1.0D0
C         USER3=DBLE(NPRTNS)
         IF (W2FAIL) THEN
            USER3 = -1.0D0
         ELSE
            USER3=DBLE(NSCLTR)
         ENDIF
      ENDIF
      IF (IFERPY.GT.1) THEN
         FERBX = DPF(1)/SQRT(W2TRY(NSCLTR)) 
         FERBY = DPF(2)/SQRT(W2TRY(NSCLTR)) 
         FERBZ = DPF(3)/SQRT(W2TRY(NSCLTR))
         FERGAM = SQRT(1.0D0+FERBX*FERBX+FERBY*FERBY+FERBZ*FERBZ)
         FERBX = FERBX/FERGAM
         FERBY = FERBY/FERGAM
         FERBZ = FERBZ/FERGAM
         DO ITRK=1,NPRTNS
            DO IDIM=1,NDIM
               P(INDXP(ITRK),IDIM)=PTMP(ITRK,IDIM)
            ENDDO
            IF (ABS(P(INDXP(ITRK),5)-PTMP(ITRK,5)).GT.EPSPF)
     &           STOP 'PFSHIFT: FATAL ERROR'
            CALL PYROBO(INDXP(ITRK),INDXP(ITRK),0.0D0,0.0D0,
     &           FERBX,FERBY,FERBZ)
         ENDDO
         IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
            WRITE(*,*)"PYLIST: After pF boost"
            CALL PYLIST(2)
         ENDIF
      ENDIF
      
C     Boost back into the TRF
      CALL PYROBO(0,0,0.0D0,0.0D0,0.0D0,0.0D0,BGCMS(2)/GACMS(2))
      
      RETURN
      END
