      SUBROUTINE EVTSUM(P5SUM,IZSUM,IASUM,ESVERB,PRHTOI)
C
C     Output variables: Sums for all stable particles:
C            P5SUM: "5"-vector sum (5=e,5=sqrts) 
C            IZSUM: Total charge
C            IASUM: Total baryon #
C     Control variable: ESVERB: Printout information?
C                       PRHTOI: Print IRF in addition to HCMS? 
C
      IMPLICIT NONE
      include "beagle.inc"
      include "dtevts.inc"

C      Note: IDRES & IDXRES are filled late for many particles.
C            Do not assume they exist unless ID=80000

      INTEGER PYCHGE,NBARY
      EXTERNAL PYCHGE,NBARY

      DOUBLE PRECISION P5SUM(5), EHTEMP, PZHTEMP
      LOGICAL ESVERB, PRHTOI
      INTEGER IZSUM,IASUM

      INTEGER ITRKES,IDIMES

      IZSUM=0
      IASUM=0
      DO IDIMES=1,4
         P5SUM(IDIMES)=0.0D0
      ENDDO
    
C     Note: Internally we use ISTHKK = +/-1 and 1001 for stable particles.
      DO ITRKES=1,NHKK
         IF (ABS(ISTHKK(ITRKES)).EQ.1 .OR. ISTHKK(ITRKES).EQ.1001) THEN
            IF (IDHKK(ITRKES).EQ.80000) THEN
               IZSUM = IZSUM + IDXRES(ITRKES)
               IASUM = IASUM + IDRES(ITRKES)
            ELSE
               IZSUM = IZSUM + PYCHGE(IDHKK(ITRKES))/3
               IASUM = IASUM + NBARY(IDHKK(ITRKES))/3
            ENDIF
            DO IDIMES=1,4
               P5SUM(IDIMES)=P5SUM(IDIMES)+PHKK(IDIMES,ITRKES)
            ENDDO
         ENDIF
      ENDDO
      P5SUM(5)=SQRT( (P5SUM(4)-P5SUM(3))*(P5SUM(4)+P5SUM(3))
     &     -P5SUM(2)*P5SUM(2)-P5SUM(1)*P5SUM(1) )
      
      IF (ESVERB) THEN
         WRITE(*,*) 'P5SUM: ',P5SUM(1),' ',P5SUM(2),' ',P5SUM(3),' ',
     &        P5SUM(4),' ',P5SUM(5)
         IF (PRHTOI) THEN
            CALL DT_LTNUC(P5SUM(3),P5SUM(4),PZHTEMP,EHTEMP,-3)
            WRITE(*,*) 'IRF:   ',P5SUM(1),' ',P5SUM(2),' ',PZHTEMP,' ',
     &        EHTEMP,' ',P5SUM(5)
         ENDIF
         WRITE(*,*) 'ZSUM, ASUM: ',IZSUM,' ',IASUM
      ENDIF

      RETURN
      END
