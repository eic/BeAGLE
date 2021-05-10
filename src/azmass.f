      DOUBLE PRECISION FUNCTION AZMASS(A,Z,IMETH)
C...
C... Mark Baker 21-Oct-2016
C...            28-Feb-2017 correct for atomic electrons
C...
C... Input A, Z and output the mass in GeV/c^2 using method IMETH
C...
C... IMETH = -1: Use 1 for backward compatibility. 
C... IMETH = 0: Best default. Currently using IMETH=3
C...
C... IMETH = 1: Lookup table based on ??, subtracting atomic electron masses,
C...            but not atomic binding energy. Note: A=1 (p & n) do NOT refer
C...            to atoms, but to the particles.
C... IMETH = 2: As above, but also adding back atomic binding energy.
C...
C... IMETH = 3: Calculate as in DT_RESNCL in DPMJET-F
C...
      IMPLICIT NONE

      INTEGER A,Z,IMETH

      EXTERNAL EXMSAZ
      DOUBLE PRECISION EXMSAZ, EXRESU

      INTEGER NTAB, IZDUM
      PARAMETER (NTAB=20)

      DOUBLE PRECISION AMU, MELEC
      PARAMETER (AMU=0.931494095)
      PARAMETER (MELEC=0.510998946D-03)

      DOUBLE PRECISION EMVGEV,AMUGEV,AMPRTN,AMNTRN,AMELCT,HLFHLF
      DOUBLE PRECISION FERTHO,BEXC12,AMUNMU,AMUC12
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( HLFHLF = 0.5D+00 )
      PARAMETER ( FERTHO = 14.33       D-09 )
      PARAMETER ( BEXC12 = FERTHO * 72.40715579499394D+00 )
      PARAMETER ( AMUNMU = HLFHLF * AMELCT - BEXC12 / 12.D+00 )
      PARAMETER ( AMUC12 = AMUGEV - AMUNMU )

      INTEGER AVAL(NTAB)
      INTEGER ZVAL(NTAB)
      DOUBLE PRECISION MASS(NTAB)

      DATA AVAL/1,1,2,3,3, 4,7,9,11,12, 27,28,56,63,107,
     &          109,131,179,208,238/
      DATA ZVAL/0,1,1,1,2, 2,3,4,5,6, 13,14,26,29,47,
     &          47,54,79,82,92/
      DATA MASS/1.0086649, 1.0072765, 2.0135532, 3.0160492, 3.0160293, 
     &          4.0026033, 7.0160046, 9.0121822, 11.009305, 12.000000, 
     &          26.981539, 27.976927, 55.934938, 62.929598, 106.90510, 
     &          108.90475, 130.90508, 196.96657, 207.97665, 238.05079/

      LOGICAL FIRSTERR
      DATA FIRSTERR /.TRUE./

      INTEGER I
C     Start with a guess that the mass is A amu
      IF (IMETH.EQ.3 .OR. IMETH.EQ.0) THEN
         EXRESU=EXMSAZ(DBLE(A),DBLE(Z),.TRUE.,IZDUM)
C         WRITE(*,*)'A,AMUC12,EMVGEV,EXMSAZ(A,Z,.TRUE.,IZDUM):',
C     &        A,AMUC12,EMVGEV,EXRESU
         AZMASS = A*AMUC12+EMVGEV*EXRESU
C         WRITE(*,*)'Calc. A*AMUC12+EMVGEV*EXMSAZ:',AZMASS
      ELSE
         AZMASS = DBLE(A)
         DO I=1,NTAB
            IF (Z.EQ.ZVAL(I)) THEN
               IF (A.EQ.AVAL(I)) THEN
                  AZMASS=MASS(I)
               ENDIF
            ENDIF
         ENDDO
         AZMASS = AZMASS*AMU
         IF (A.GT.1) THEN
            AZMASS=AZMASS-DBLE(Z)*MELEC
            IF (IMETH.EQ.2) THEN
               IF (Z.EQ.6 .AND. A.EQ.12) THEN
                  AZMASS = AZMASS + 632.16D-09
               ELSEIF (Z.EQ.79 .AND. A.EQ.197) THEN
                  AZMASS = AZMASS + 322187.9637D-09
               ELSEIF (Z.EQ.82 .AND. A.EQ.208) THEN
                  AZMASS = AZMASS + 355823.5114D-09
               ELSEIF (FIRSTERR) THEN
                  WRITE(*,*) 'AZMASS Warning. Method 2 not allowed'
     &                 //"w/ this A,Z:",A," ",Z,". Using method 3."
                  EXRESU=EXMSAZ(DBLE(A),DBLE(Z),.TRUE.,IZDUM)
                  AZMASS = A*AMUC12+EMVGEV*EXRESU
                  FIRSTERR = .FALSE.
               ENDIF
            ENDIF
         ENDIF
      ENDIF

C      WRITE(*,*) 'Using IMETH, A, Z, AZMASS (GeV)',IMETH,A,Z,AZMASS
      RETURN
      END

