      DOUBLE PRECISION FUNCTION GETESTAR(ESFIXTYP,ESFPAR1,ESFPAR2)
      DOUBLE PRECISION ESFPAR1,ESFPAR2
      INTEGER ESFIXTYP
      DOUBLE PRECISION R1,R2,DT_RNDM
      EXTERNAL DT_RNDM
C...
C... Mark D. Baker 05-Jun-2023
C...
C... Roll desired E* value from either gaussian (ESFIXTYP=1) or
C... flat (ESFIXTYP=2) distribution. 
C... ESFPAR1 = center of gaussian or left edge of flat
C... ESFPAR2 = sigma of gaussian or right edge of flat

      IF (ESFIXTYP.EQ.1) THEN
 100     CALL DT_RANNOR(R1,R2)
         GETESTAR = ESFPAR1 + ESFPAR2*R1
         IF (GETESTAR.LT.0.0) GETESTAR = ESFPAR1 + ESFPAR2*R2
         IF (GETESTAR.LT.0.0) GOTO 100
      ELSEIF (ESFIXTYP.EQ.2) THEN
         GETESTAR = ESFPAR1 + (ESFPAR2-ESFPAR1)*DT_RNDM(R1)
      ELSE
         STOP "FATAL ERROR in getestar.f. Illegal ESFIXTYP."
      ENDIF

      IF (GETESTAR.NE.GETESTAR) THEN
         WRITE (*,*) 'NaN in GETESTAR!'
         WRITE (*,*) 'GETESTAR,ESFPAR1,ESFPAR2',GETESTAR,ESFPAR1,ESFPAR2
      ENDIF

      RETURN
      END
