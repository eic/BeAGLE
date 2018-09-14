       DOUBLE PRECISION FUNCTION pyth_xsec(dy,dx,dQ2,dF1,dF2)

       implicit none
        include "mc_set.inc"

*-----DEKLARATION----------------------

       DOUBLE PRECISION DX, DQ2,dy,DEBEA    ! kinematic variables
       DOUBLE PRECISION DBRACK   ! intermed. factor
       double precision emass2,alpha2,pi
       double precision df1,df2
! Load up some convenient variables
        pi=3.1415926d0
        alpha2=0.729735d-2**2
        emass2=masse**2
        DEBEA = dble(mcSet_EneBeam)
! Check that the kinematics are reasonable
        if (dx.eq.0D0) then
          goto 999
        endif
       IF(DX.LE.0.D0.OR.DX.GE.1.D0.OR.
     1    Dy.LE.0.D0.OR.Dy.GE.1.D0.OR.dq2.le.0) then
          goto 999
       endif
! Determine the UNPOLARIZED xsec

       DBRACK = (4.D0*pi*alpha2)/dq2**2
! ........ that's dSigma/dlogQ2/dlogy

       pyth_xsec = DBRACK*(dy**2*(1-2D0*emass2/dq2)*dF1+
     1                 (1-dy-dq2/(4D0*DEBEA**2))*dF2/dx)*dx*dQ2

! Leave by the front door

       RETURN

! Leave by the back door

999    continue
       pyth_xsec = 0.0D0
       return
       end

