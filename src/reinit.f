      SUBROUTINE REINIT(NEWID)

C     Mark D. Baker  2018-02-24 Initial Version
C
C     Reinitializes Pythia and DPMJET for ep vs. en hard collisions
C
      include "pythia.inc"
      include "mc_set.inc"
      include "beagle.inc"

      INTEGER NEWID

      INTEGER NCALLS
      DATA NCALLS /0/
      SAVE NCALLS

      idNucPY=NEWID
      Mnucl=PYMASS(idNucPY)
      massp=real(PYMASS(idNucPY))
      P(1,1) = 0.0
      P(1,2) = 0.0
      P(1,3) = PZLEP
      P(2,1) = 0.0
      P(2,2) = 0.0

      IF (NEWID.EQ.2212) THEN
         idNucBAM=1
         pbeamdbl = pbeamP
         PZNUCL   = pbeamP
         P(2,3)   = pbeamP
         pbeam    = real(pbeamP)
         call pyinit('3MOM', lName, 'p+', WIN)
      ELSEIF (NEWID.EQ.2112) THEN
         idNucBAM=8
         pbeamdbl = pbeamN
         PZNUCL   = pbeamN
         P(2,3)   = pbeamN
         pbeam    = real(pbeamN)
         call pyinit('3MOM', lName, 'n0', WIN)
      ELSE
         WRITE(*,*) "ERROR: Can't initialize with ID: ",NEWID
      ENDIF
 
C Note: Radgen not tested in BeAGLE
C Does it need reinitializing?
      if (iModel.eq.1 .and. qedrad.eq.1) then
         MSTP(199)=1
         call radgen_init(.true.,.false.)
      endif

      NCALLS = NCALLS + 1
      IF (NCALLS.EQ.4) MSTP(122)=0
      RETURN
      END
