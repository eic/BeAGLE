      INTEGER FUNCTION NBARY(KF)
      INTEGER KF, AKF, TENS
C
C... Mark Baker 06-Oct-2016 
C
C... Return 3 x baryon # given KF code
C...
C... Quarks are 1-8
C... Diquarks are between 1100 and 5999 with 10s digit = 0
C... Baryons are between 1100 and 5999 and non-zero 10s digit
C... Anti- quarks, diquarks, baryons are the negative of this

      NBARY = 0
      AKF = ABS(KF)
      IF (1.LE.AKF .AND. AKF.LE.8) THEN
         NBARY = 1
      ELSEIF (1100.LE.AKF .AND. AKF.LE.5999) THEN
         TENS = MOD(AKF,100) / 10
         IF (TENS.EQ.0) THEN
            NBARY = 2
         ELSE
            NBARY = 3
         ENDIF
      ENDIF
      NBARY = SIGN(NBARY, KF)

      RETURN
      END

