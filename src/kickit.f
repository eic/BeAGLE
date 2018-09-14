      SUBROUTINE KICKIT(PXPRT,PYPRT,PZPRT,P0PRT,MPRT,
     &                  PXNUC,PYNUC,PZNUC,P0NUC,MNUC)
C
C... Mark Baker 2016-08-24
C... Roll an intrinsic kt and apply it to the main struck parton
C... Apply negative kt (recoil) to the nucleon
C... Modify pz, p0 to maintain 4-momentum conservation (elastic scattering)
      DOUBLE PRECISION PXPRT, PYPRT, PZPRT, P0PRT, MPRT         
      DOUBLE PRECISION PXNUC, PYNUC, PZNUC, P0NUC, MNUC
      DOUBLE PRECISION PTK1, PTK2, PZSUM, P0SUM, PXSUM, PYSUM, MSUM
      DOUBLE PRECISION MT2PRT, MT2NUC, MU2         
      DOUBLE PRECISION DETERM, DETEMP
      DOUBLE PRECISION MPRTCH, MNUCCH

      PXSUM = PXPRT + PXNUC
      PYSUM = PYPRT + PYNUC
      PZSUM = PZPRT + PZNUC
      P0SUM = P0PRT + P0NUC 
      MSUM  = DSQRT(P0SUM*P0SUM-PXSUM*PXSUM-PYSUM*PYSUM-PZSUM*PZSUM)
      MPRTCH = DSQRT(P0PRT*P0PRT-PXPRT*PXPRT-PYPRT*PYPRT-PZPRT*PZPRT)
      MNUCCH = DSQRT(P0NUC*P0NUC-PXNUC*PXNUC-PYNUC*PYNUC-PZNUC*PZNUC)
C      WRITE(*,*) 'KICKIT: Before picture:'
C      WRITE(*,*) '          px       py      pz      E      M      Mchk'
C      WRITE(*,*) 'parton:  ',PXPRT,' ',PYPRT,' ',PZPRT,' ',P0PRT,' ',
C     &                       MPRT, ' ', MPRTCH
C      WRITE(*,*) 'nucleon: ',PXNUC,' ',PYNUC,' ',PZNUC,' ',P0NUC,' ',
C     &                       MNUC, ' ', MNUCCH
C      WRITE(*,*) 'sum:     ',PXPRT+PXNUC,' ',PYPRT+PYNUC, PZSUM, ' ',
C     &                       P0SUM,' ',MSUM
C      WRITE(*,*) ' '
      CALL RMKICK(PTK1,PTK2)
      PXPRT = PXPRT + PTK1
      PYPRT = PYPRT + PTK2
      PXNUC = PXNUC - PTK1
      PYNUC = PYNUC - PTK2
      MT2PRT = MPRT*MPRT + PXPRT*PXPRT + PYPRT*PYPRT
      MT2NUC = MNUC*MNUC + PXNUC*PXNUC + PYNUC*PYNUC
      MU2 = 0.5D0 * (P0SUM*P0SUM - PZSUM*PZSUM + MT2NUC - MT2PRT)
      DETEMP = 1.0D0 - MT2NUC*(P0SUM*P0SUM-PZSUM*PZSUM)/MU2/MU2
      IF (DETEMP .GE. 0.0D0) THEN
         DETERM = MU2 * P0SUM * DSQRT(DETEMP)
C Note: The + DETERM root corresponds to the parton stopping and the
C       giving most of its momentum to the nucleon. We don't want that.
         PZNUC  = (PZSUM*MU2 - DETERM ) / (P0SUM*P0SUM-PZSUM*PZSUM)
         P0NUC  = DSQRT(PZNUC*PZNUC   + MT2NUC)
         P0PRT  = P0SUM - P0NUC
         PZPRT  = PZSUM - PZNUC
      ELSE
C FAILED! Restore inputs
         PXPRT = PXPRT - PTK1
         PYPRT = PYPRT - PTK2
         PXNUC = PXNUC + PTK1
         PYNUC = PYNUC + PTK2
         WRITE(*,*)"KICKIT ERROR: Kick pT violates energy conservation"
         WRITE(*,*) '          px      py     pz     E     M'
         WRITE(*,*) 'parton:  ',PXPRT,PYPRT,PZPRT,P0PRT,MPRT 
         WRITE(*,*) 'nucleon: ',PXNUC,PYNUC,PZNUC,P0NUC,MNUC 
         WRITE(*,*) 'PT Kick: ',PTK1, PTK2
      ENDIF         
C         MPRTCH = DSQRT(P0PRT*P0PRT-PXPRT*PXPRT-PYPRT*PYPRT-PZPRT*PZPRT)
C         MNUCCH = DSQRT(P0NUC*P0NUC-PXNUC*PXNUC-PYNUC*PYNUC-PZNUC*PZNUC)
C         MSUM  = DSQRT(P0SUM*P0SUM-PXSUM*PXSUM-PYSUM*PYSUM-PZSUM*PZSUM)
C      WRITE(*,*) 'KICKIT: After picture:'
C      WRITE(*,*) '          px       py      pz      E      M      Mchk'
C      WRITE(*,*) 'parton:  ',PXPRT,' ',PYPRT,' ',PZPRT,' ',P0PRT,' ',
C     &                       MPRT, ' ', MPRTCH
C      WRITE(*,*) 'nucleon: ',PXNUC,' ',PYNUC,' ',PZNUC,' ',P0NUC,' ',
C     &                       MNUC, ' ', MNUCCH
C      WRITE(*,*) 'sum:     ',PXPRT+PXNUC,' ',PYPRT+PYNUC, PZPRT+PZNUC,
C     &                       ' ',P0PRT+P0NUC,' ', MSUM
C      WRITE(*,*) ' '
   
      RETURN
      END
