      SUBROUTINE PFSHIFT(WRAW,W2RAW,Q2,NU,MNUCL,PDEUT,PSPEC,IKIN,IREJ)
C
C     2018-07-13 Mark D. Baker - Initial Version
C
C     This subroutine adds the struck nucleon pF back to the
C     Pythia subevent partonic skeleton after the fact.
C     We start and end in the TRF, but operate in the naive HCMS
C     defined by gamma* + N at rest in nuclear TRF.
C
C     NOTE: For now, we are ignoring the possibility of PDEUT(1-3)
C           being non-zero
C
C     Input: WRAW, W2RAW - from Pythia subevent VINT(1),VINT(2)
C            Q2, NU - from Pythia
C            MNUCL - struck nucleon for IKIN=1 (VINT(4)), spectator for IKIN=2
C            PDEUT - Deuteron "5"-momentum, irrelevant for IKIN=1
C            PSPEC - Irrelevant for IKIN=1
C                    Spectator "5"-momentum, for IKIN=2
C                    "5"-momentum of sum of spec1+spec2 for IKIN=3
C            IKIN - model for correct Pythia subsystem Wmu   
C
C     I/O: IREJ = error flag.
C          Input: Expect IREJ=0. Otherwise skip.
C          Output: 0=OK. 1=Kill event: bad kinematics
C
C     Common block input: P(mu)_true-P(mu)_naive OR use k
C                         PXF,PYF,PZF,EKF = k & E(k)-m in the TRF
C                         DPF(4) is in the naive HCMS 
C
C     There are two models for the "correct" Pythia subsystem Wmu.
C     In both cases, q comes from Pythia.
C     IKIN=1: incoming struck nucleon p on mass shell but with additional
C              3-Fermi momentum. Then use W=p+q
C     IKIN=2: Appropriate especially for SRC & deuteron case. Could be used
C              for bulk Fermi momentum too...
C              W = q + P - pspec'
C              where P is an incoming on-shell deuteron 
C              pspec' given by an on-shell spectator with 3-momentum -k
C              For A=3 case, W = q + P - pspec1' - pspec2' and 
C              pspec' = pspec1' + pspec2'
C     IKIN=0: Calculate IKIN=1, for filling in USERSET 3, but don't
C     actually change the event record.
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
      IMPLICIT NONE
      DOUBLE PRECISION WRAW, W2RAW, Q2, NU, MNUCL, PDEUT(5), PSPEC(5)
      INTEGER IKIN, IREJ

      DOUBLE PRECISION W2MIN
      PARAMETER (W2MIN=4.0D0)

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
      PARAMETER (MAXPRTS=50)
      PARAMETER (NDIM=4)
      PARAMETER (NDIMM=5)
      DOUBLE PRECISION W2F, W2TRY(MAXTRY), PSUM(NDIM) 
      DOUBLE PRECISION ASCALE(MAXTRY), PHIGH2, SQRM1, SQRM2, ASCLFL
      DOUBLE PRECISION ASCTMP1, W2TTMP1, S2SUM 
      DOUBLE PRECISION FERBX, FERBY, FERBZ, FERGAM
      INTEGER NPRTNS,NLSCAT,IDIM,NSCLTR, ITRK
      DOUBLE PRECISION PTMP(MAXPRTS,NDIMM)
      INTEGER INDXP(MAXPRTS)
      LOGICAL W2FAIL
      DOUBLE PRECISION WZIRF,W0IRF,WZHCMS,W0HCMS

      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) then
         WRITE(*,*)
     &        "PFSHIFT(WRAW,W2RAW,Q2,NU,MNUCL,PDEUT,PSPEC,IKIN,IREJ)"
         WRITE(*,*)WRAW,W2RAW,Q2,NU,MNUCL
         WRITE(*,*)PDEUT(1),PDEUT(2),PDEUT(3),PDEUT(4),PDEUT(5)
         WRITE(*,*)PSPEC(1),PSPEC(2),PSPEC(3),PSPEC(4),PSPEC(5)
         WRITE(*,*) IKIN,IREJ
      ENDIF

C     Skip event if it is already flagged as bad (shouldn't happen!)
      IF (IREJ.NE.0) THEN
         WRITE(*,*)"PHSHIFT WARNING: IREJ non-zero on input. Skipping."
         GOTO 9999
      ENDIF

      IF (WRAW.LT.1.0 .OR. ABS(W2RAW-WRAW*WRAW).GT.0.001 .OR.
     & MNUCL.LT.0.9 .OR. MNUCL.GT.1.0) THEN
         WRITE(*,*)'WRAW,W2RAW,MNUCL: ',WRAW,W2RAW,MNUCL
         STOP 'PFSHIFT: FATAL ERROR. BAD KINEMATICS!'
      ENDIF
C     Boost into naive HCMS  (assumes nucleon at rest in A-TRF)
      CALL PYROBO(0,0,0.0D0,0.0D0,0.0D0,0.0D0,-BGCMS(2)/GACMS(2))

      IF (IKIN.EQ.0 .OR. IKIN.EQ.1) THEN
         W2F = W2RAW + 2.*WRAW*DPF(4) - 2.*MNUCL*EKF
      ELSEIF (IKIN.EQ.2) THEN
         W2F = 2.*(PDEUT(4)-PSPEC(4))*NU - Q2 + PDEUT(5)*PDEUT(5) + 
     &        PSPEC(5)*PSPEC(5) - 2.*PSPEC(4)*PDEUT(4) 
     &        - 2*(PDEUT(3)-PSPEC(3))*SQRT(Q2+NU*NU) +
     &        2.*(PDEUT(1)*PSPEC(1)+PDEUT(2)*PSPEC(2)+
     &            PDEUT(3)*PSPEC(3))
      ELSE
         WRITE(*,*) 'PFSHIFT ERROR: Illegal IKIN =',IKIN
         W2F = W2RAW
      ENDIF

      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         write(*,*) "W2 ignore pF:", W2RAW
         write(*,*) "W2 corrected: ", W2F
         write(*,*) "gamma*beta, gamma, beta:",BGCMS(2),GACMS(2),
     &        BGCMS(2)/GACMS(2)
      ENDIF

C     Don't allow illegal kinematics (usually spec. kz<<0)
      IF (W2F.LT.W2MIN) THEN
         IF (IOULEV(1).GT.0 .OR. 
     &        (IOULEV(4).GE.1 .AND. NEVENT.LE.IOULEV(5))) THEN
            WRITE(*,*) 'PFSHIFT: Illegal kinematics. W2RAW, W2F= ',
     &           W2RAW, W2F
            WRITE(*,*) 'PSPEC(5)=',
     &           PSPEC(1),PSPEC(2),PSPEC(3),PSPEC(4),PSPEC(5)
         ENDIF
         IREJ=1
         GOTO 9999
      ENDIF
      
      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         write(*,*) "DT_PYEVNTEP: HCMS g*=z, e' px>0 py=0"
         CALL PYLIST(2)
      ENDIF

C     Analyze the events and extract the "partonic" skeleton
      NPRTNS=0
      NLSCAT=0
      S2SUM=0.0D0
      DO IDIM=1,NDIM
         PSUM(IDIM)=0.0D0       ! sum p^mu for all stable except e'
      ENDDO
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
C
      IF (IKIN.GT.0) THEN
         FERBX = DPF(1)/SQRT(W2TRY(NSCLTR)) 
         FERBY = DPF(2)/SQRT(W2TRY(NSCLTR)) 
         IF (IKIN.EQ.1) THEN
            FERBZ = DPF(3)/SQRT(W2TRY(NSCLTR))
         ELSEIF (IKIN.EQ.2) THEN
C     Note: Wmu transverse quantities are PXF,PYF=DPF(1,2). WZ,W0:
            WZIRF = SQRT(Q2+NU*NU)+PZF
            W0IRF = NU + PDEUT(5) - PSPEC(4)
            CALL DT_LTNUC(WZIRF,W0IRF,WZHCMS,W0HCMS,3)
            FERBZ = WZHCMS/SQRT(W2TRY(NSCLTR))
         ELSE
            FERBZ = 0.0D0
            FERBY = 0.0D0
            FERBX = 0.0D0
         ENDIF
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
      
 9999 RETURN
      END
