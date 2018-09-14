C*********************************************************************
 
C... RMKICK - Mark 2016-08-18 
C... Stripped down version of PYREMN to apply primordial kt recoil
C... to extra interactions in multinucleon eA scattering
C... Note: Includes liang's changes to allow double gaussian.
C...
      SUBROUTINE RMKICK(PT1,PT2)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/
C...Local arrays.
      DIMENSION KFLCH(2),KFLSP(2),CHI(2),PMS(0:6),IS(2),ISN(2),ROBO(5),
     &PSYS(0:2,5),PMIN(0:2),QOLD(4),QNEW(4),DBE(3),PSUM(4)
 
      IF(MSTP(91).LE.0) THEN
         PT=0D0
      ELSEIF(MSTP(91).EQ.1) THEN
         PT=PARP(91)*SQRT(-LOG(PYR(0)))
         !added by liang to implement double gauss kt++
      ELSEIF(MSTP(91).EQ.5) THEN
         !sigma^2 for two rayleigh
         WDTH1=PARP(91)**2
         WDTH2=PARP(195)**2
         !relative factor for second rayleigh
         RF=PARP(196)
         RATIO=RF*WDTH2/(WDTH1+RF*WDTH2)
         IF(PYR(0).lt.RATIO) THEN
            PT=PARP(195)*SQRT(-LOG(PYR(0)))
         ELSE
            PT=PARP(91)*SQRT(-LOG(PYR(0)))
         ENDIF
         !added by liang to implement double gauss kt--
      ELSE
         RPT1=PYR(0)
         RPT2=PYR(0)
         PT=-PARP(92)*LOG(RPT1*RPT2)
      ENDIF
      PHI=PARU(2)*PYR(0)
      PT1=PT*COS(PHI)
      PT2=PT*SIN(PHI)
 
      RETURN
      END
