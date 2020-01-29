      SUBROUTINE RGTOPY
C
C... Mark Baker 2019-01-27
C... 
C... Modify RG event to better match Pythia so that we can use 
C... BeAGLE routines which are designed for Pythia.
C...
C... Swap Pythia particles 3 and 4
C... Add some missing mother-daughter relationships.
C... In particular, BeAGLE relies on K(3,3)=1 to identify the e'
C
      include 'pythia.inc'
      INTEGER KTMP(5)
      DOUBLE PRECISION PTMP(5),VTMP(5)
      WRITE(*,*) "RGTOPY: PYLIST before:"
      CALL PYLIST(2)
      DO I=1,5
         KTMP(I)=K(3,I)
         PTMP(I)=P(3,I)
         VTMP(I)=V(3,I)
         K(3,I)=K(4,I)
         P(3,I)=P(4,I)
         V(3,I)=V(4,I)
         K(4,I)=KTMP(I)
         P(4,I)=PTMP(I)
         V(4,I)=VTMP(I)
      ENDDO
      K(1,4) = 3
      K(1,5) = 4
      K(2,4) = 5
       WRITE(*,*) "RGTOPY: PYLIST after:"
      CALL PYLIST(2)
      RETURN
      END
