      SUBROUTINE SCLSUBEVT(W2RAW,W2F,EPS,VERB,NPARTS,PTMP,W2OUT,NSCLTR,
     &     W2FAIL)
C           out  )        ( in    in in   in   in    I/O  out    out  
C
C     2023-06 Mark D. Baker - Initial version
C
C     Scale a set of particle 3-momenta in their own cm
C     in order to change their overall s
C
C     Originally part of PFSHIFT, but also needed in FIXESTAR
C
C W2RAW = Invariant mass^2 s of the particles in PTMP
C W2F = Desired invariant s for the particles
C EPS = Precision required
C VERB = Verbosity 
C NPARTS = # of particles 
C PTMP = array of particle "5-momenta" = 4 momentum & mass
C W2OUT = final s for the particles 
C NSCLTR = # of interations needed to converge
C W2FAIL = .TRUE. if we fail to converge, .FALSE. if we succeed.
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
C S2SUM = sum(p^2/E) for all particles
C 

C     Note: MAXTRY, MAXPRTS should match those in the calling routines

      IMPLICIT NONE

C     Constants and calling parameters
      LOGICAL VERB, W2FAIL
      INTEGER NPARTS, NSCLTR 
      INTEGER MAXTRY, MAXPRTS, NDIM, NDIMM
      PARAMETER (MAXTRY=10, MAXPRTS=1000)
      PARAMETER (NDIM=4, NDIMM=5)
      DOUBLE PRECISION EPS,W2RAW,W2F,W2OUT,PTMP(MAXPRTS,NDIMM)

C     Local variables
      INTEGER ITRK, IDIM, II
      DOUBLE PRECISION ASCALE(MAXTRY), W2TRY(MAXTRY), ASCLFL, S2SUM
      DOUBLE PRECISION WRAW, PHIGH2, SQRM1, SQRM2, PSUM(NDIM)

C     Protect against NaNs
      IF (W2RAW.LT.0.0) STOP ('SCLSUBEVT FATAL: W2RAW<0')
      IF (W2F.LT.0.0) STOP ('SCLSUBEVT FATAL: W2F<0')

      WRAW=DSQRT(W2RAW)
      S2SUM = 0.0D0
      DO ITRK=1,NPARTS
         S2SUM = S2SUM + (PTMP(ITRK,4)-PTMP(ITRK,5))
     &        *(1.0D0+PTMP(ITRK,5)/PTMP(ITRK,4))
      ENDDO

      IF (NPARTS.GT.2) THEN
         ASCALE(1) = 1.0D0 + (SQRT(W2F)-WRAW)/S2SUM
      ELSEIF (NPARTS.EQ.2) THEN
         PHIGH2 =( PTMP(1,1)**2+PTMP(1,2)**2+PTMP(1,3)**2 +
     &             PTMP(2,1)**2+PTMP(2,2)**2+PTMP(2,3)**2 )/2.0D0
         SQRM1 = PTMP(1,5)*PTMP(1,5)
         SQRM2 = PTMP(2,5)*PTMP(2,5)
         ASCALE(1) = (W2F-2.0D0*(SQRM1+SQRM2)+((SQRM1-SQRM2)**2)/W2F )/
     &        (4.0D0*PHIGH2)
         IF (ASCALE(1).LT.0.0) THEN
            WRITE(*,*) 'FATAL ERROR IN SCLSUBEVT. ASCALE(1)**2<0'
            WRITE(*,*) 'PTMP(1,*):',(PTMP(1,II),ii=1,5)
            WRITE(*,*) 'PTMP(2,*):',(PTMP(2,II),ii=1,5)
            WRITE(*,*) 'W2F,SQPM1,SQRM2,PHIGH2',W2F,SQRM1,SQRM2,PHIGH2
            STOP('STOPPING')
         ELSE
            ASCALE(1)=SQRT(ASCALE(1))
         ENDIF
      ELSE
         STOP ('PFSHIFT: FATAL ERROR. NPARTS<2!')
      ENDIF
      IF (VERB) THEN
         WRITE(*,*)'True subevent CMS. Before W2 rescale. p5:'
         DO ITRK=1,NPARTS
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
      DO ITRK=1,NPARTS
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
      IF (VERB) THEN
         WRITE(*,*)"WF/WRAW, ASCALE(1):", SQRT(W2F/W2RAW),ASCALE(1)
         WRITE(*,*)"W2F:",W2F,"Scaled W2:",W2TRY(1),"W2try/W2F",
     &        W2TRY(1)/W2F,"NPARTS:",NPARTS
         WRITE(*,*)'HCMS. After pF-based W2 rescale'
         DO ITRK=1,NPARTS
            WRITE(*,*)ITRK," ",PTMP(ITRK,1)," ",PTMP(ITRK,2)," ",
     &           PTMP(ITRK,3)," ",PTMP(ITRK,4)," ",PTMP(ITRK,5)
         ENDDO
         WRITE(*,*)"PSUM(1-4):",PSUM(1),PSUM(2),PSUM(3),PSUM(4)
      ENDIF
C     2nd and subsequent iterations
      NSCLTR=1
      DO WHILE (NSCLTR.LT.MAXTRY .AND.
     &     ABS(W2TRY(NSCLTR)/W2F-1.0D0).GT.EPS)
         NSCLTR=NSCLTR+1
         PHIGH2 = ASCALE(NSCLTR-1)*ASCALE(NSCLTR-1)*PHIGH2
         IF (W2TRY(NSCLTR-1).LT.0.0) THEN
            WRITE (*,*) 'SCLSUBEVT ERROR: W2TRY<0:, ',W2TRY(NSCLTR-1)
            STOP
         ENDIF
         ASCALE(NSCLTR) = 1.0D0+(SQRT(W2F)-SQRT(W2TRY(NSCLTR-1)))/S2SUM
         IF ( VERB .OR.  NSCLTR.GT.MAXTRY-3 )  THEN
            WRITE(*,*)'W2 inaccurate. Iteration # ',NSCLTR
            WRITE(*,*)'Next ASCALE factor = ',ASCALE(NSCLTR)
         ENDIF
         S2SUM=0.0D0
C     3-momenta just scale
         DO IDIM=1,3
            PSUM(IDIM)=ASCALE(NSCLTR)*PSUM(IDIM)
         ENDDO
         PSUM(4)=0.0D0
         DO ITRK=1,NPARTS
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
      
         IF ( VERB .OR.  NSCLTR.GT.MAXTRY-3 ) THEN
            WRITE(*,*)"W2F:",W2F,"Iteration ",NSCLTR," Scaled W2:",
     &           W2TRY(NSCLTR),"W2try(",NSCLTR,")/W2F",
     &           W2TRY(NSCLTR)/W2F,"NPARTS:",NPARTS
            WRITE(*,*)'HCMS. After W2 rescale iteration ', NSCLTR
            DO ITRK=1,NPARTS
               WRITE(*,*)ITRK," ",PTMP(ITRK,1)," ",PTMP(ITRK,2)," ",
     &              PTMP(ITRK,3)," ",PTMP(ITRK,4)," ",PTMP(ITRK,5)
            ENDDO
            WRITE(*,*)"Leads to ..."
            WRITE(*,*)"PSUM(1-4):",PSUM(1),PSUM(2),PSUM(3),PSUM(4)
         ENDIF
      ENDDO

      W2FAIL = (ABS(W2TRY(NSCLTR)/W2F-1).GT.EPS)

      IF ( VERB .OR.  NSCLTR.GT.MAXTRY-3 ) THEN
         WRITE(*,*)
         WRITE(*,*)'Iteration   W2/W2F'
         WRITE(*,*)'          0   ',W2RAW/W2F
         ASCLFL = 1.0D0
         DO ITRK=1,NSCLTR
            WRITE(*,*) ITRK, "  ", W2TRY(ITRK)/W2F
            ASCLFL = ASCLFL * ASCALE(ITRK)
         ENDDO
         IF (W2FAIL) THEN
            WRITE(*,*)'Failed to converge to level of: ',EPS
         ELSE
            WRITE(*,*)'Success. Converged to level of: ',EPS
         ENDIF
         WRITE(*,*)'WF/WRAW:     ',SQRT(W2F)/WRAW
         WRITE(*,*)'ASCALE(1):   ',ASCALE(1)
         WRITE(*,*)'ASCALE(full):',ASCLFL
      ENDIF

      W2OUT = W2TRY(NSCLTR)
      RETURN
      END
