      DOUBLE PRECISION FUNCTION ANEAR(A)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION A
C...
C... Mark Baker 03-Oct-2017  Find nearest valid A for EPS09LO
C...

      INTEGER NAVALS
      PARAMETER (NAVALS=17)

      INTEGER ALAST, ANLAST

      INTEGER AVAL(NAVALS)
      REAL EDGE(NAVALS+1)  ! EDGE(i) is least A assigned to AVAL(i)

C     Edge(i) = sqrt(AVAL(i) * AVAL(i-1))  (AVAL(0)=1)
      DATA ALAST/ 0 /
      DATA ANLAST / 0 /
      DATA AVAL/  4,   6,   9,   12,   16,   27,   40,   56,   64, 108, 
     &          115, 117, 184,  195,  197,  208, 238 /
      DATA EDGE/1.9, 4.9, 7.3, 10.4, 13.9, 20.8, 32.9, 47.3, 59.9, 83.1,
     &          111.4, 115.99, 146.7, 189.4, 195.99, 202.4, 222.5, 999./

      IF (A.LE.0 .OR. A.GE.EDGE(NAVALS+1)) THEN
         WRITE(*,*) 'A=',A
         STOP "ERROR in ANEAR. Invalid A"
      ELSEIF (NINT(A).EQ.ALAST) THEN
         ANEAR = DBLE(ANLAST)
      ELSEIF (NINT(A).EQ.1) THEN
         ANEAR = A
      ELSE 
         ANEAR = 0.0D0
         I = 1
         DO WHILE (I.LE.NAVALS .AND. NINT(ANEAR).EQ.0) 
            IF (NINT(A).EQ.AVAL(I)) THEN
               ANEAR = A
            ELSEIF (EDGE(I).LE.A .AND. A.LE.EDGE(I+1)) THEN
               ANEAR = DBLE(AVAL(I))
               WRITE(*,*)'EPS09LO unavailable for A=',NINT(A),
     &              'using',NINT(ANEAR)
            ENDIF
            I = I + 1
         ENDDO
      ENDIF

      ALAST = NINT(A)
      ANLAST = NINT(ANEAR)
      RETURN
      END

