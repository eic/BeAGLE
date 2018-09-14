      FUNCTION ULANGL(X,Y)  

C... Purpose: to reconstruct an angle from given x and y coordinates.    
C... Taken from sPHENIX GitHub 2016-08-20. Changed LUDAT1 to PYDAT1.
C... Note: PARU(1) is of course still pi

      DOUBLE PRECISION X,Y,R,ULANGL
      INTEGER MSTU,MSTJ
      DOUBLE PRECISION PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /PYDAT1/ 
    
      ULANGL=0.D0 
      R=SQRT(X**2+Y**2) 
      IF(R.LT.1E-20) RETURN 
      IF(ABS(X)/R.LT.0.8) THEN  
        ULANGL=SIGN(ACOS(X/R),Y)    
      ELSE  
        ULANGL=ASIN(Y/R)    
        IF(X.LT.0..AND.ULANGL.GE.0.) THEN   
          ULANGL=PARU(1)-ULANGL 
        ELSEIF(X.LT.0.) THEN    
          ULANGL=-PARU(1)-ULANGL    
        ENDIF   
      ENDIF 
    
      RETURN    
      END 
