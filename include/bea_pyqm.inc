C     MDB 2017-10-06 BeAGLE-PyQM common blocks 
C     PYQM recoil for BeAGLE
      double precision PYQREC
      COMMON /CORECT/ PYQREC(4)

C     PYQM control from BeAGLE
      double precision PQRECF, PYQ_SUPF, PYQ_IET
      integer PYQ_IPTF, PYQ_IEG, PYQ_HQ, PYQ_VERB
      COMMON /PQCTRL/ PQRECF, PYQ_SUPF, PYQ_IPTF, PYQ_IEG, PYQ_HQ,
     &                PYQ_IET, PYQ_VERB

* position of interacted nucleon from BeAGLE for PyQM
      DOUBLE PRECISION PosNuc
      COMMON /NPARINT/ PosNuc(4)

