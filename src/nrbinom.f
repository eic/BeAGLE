      FUNCTION NRBINOM(N,P,EXACT)
C
C     2017-05-30 MDB 
C
C     Roll a # according to the binomal distribution for N tries
C     with probability P. For large values of np(1-p) approximate
C     with a Gaussian unless EXACT=.TRUE.
C
      IMPLICIT NONE
      INTEGER NRBINOM, N, I
      DOUBLE PRECISION P,VAR,MU,PT,PHI,PYR,TWOPI
      LOGICAL EXACT
      EXTERNAL PYR

      PARAMETER (TWOPI = 6.283185307179586d0)

      MU = DBLE(N)*P
      VAR = MU*(1.0-P)
C     
      IF (EXACT .OR. VAR.LT.20.0) THEN
         NRBINOM=0
         DO I=1,N
            IF (PYR(0).LE.P) NRBINOM = NRBINOM+1
         ENDDO
      ELSE
C        Gaussian approximation
         PT = SQRT(VAR)*SQRT(-LOG(PYR(0)))
         PHI = TWOPI*PYR(0)
         NRBINOM = NINT(MU+PT*COS(PHI))
         IF (NRBINOM.LT.0) THEN
            NRBINOM = 0
         ELSEIF (NRBINOM.GT.N) THEN
            NRBINOM = N
         ENDIF
      ENDIF

      RETURN
      END
