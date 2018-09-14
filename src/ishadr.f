      LOGICAL FUNCTION ISHADR(KF)
      INTEGER KF, AKF, TENS
C
C... Mark Baker 02-JAN-2017 
C
C... Return .TRUE. if a given KF code refers to a conventional hadron.
C... Negative numbers are antiparticles.
C...
C... 1-8    NO quark 
C... 11-18  NO leptons 
C... 21-42  NO force particles
C... 45-46  NO supersymmetric force particles? H_30, A_20
C... 81-100 NO Non-standard (strings, junctions, etc.)
C... 110    NO reggeon
C... 111-553 YES mesons
C... 990    NO pomeron
C... 1100-5999 
C...        10s digit = 0 NO diquarks
C...        10s digit > 0 YES baryons
C... 10000-99999 YES? more mesons
C... 100443,100553, YES? psi', Upsilon'
C... 1xxxxxx NO SUSY states
C... 2xxxxxx NO SUSY states
C... 3xxxxxx NO Technicolor states
C... 4xxxxxx NO Excited fermions
C... 5xxxxxx NO Excited graviton
C... 99xxxxx NO more unusual states

      AKF = ABS(KF)
      TENS = MOD(AKF,100) / 10
      ISHADR = ( (111.LE.AKF .AND. AKF.LE.599) .OR.
     &           (1100.LE.AKF .AND. AKF.LE.5999 .AND. TENS.GT.0) .OR.
     &           (10000.LE.AKF .AND. AKF.LE.999999) )

      RETURN
      END

