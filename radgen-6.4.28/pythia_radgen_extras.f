C A collection of routines needed by RADGEN
 
!---------------------------------------------------------------
! Calculate the F2 structure function

      subroutine mkf2 (DQ2,DX,A,Z,DF2,DF1)

      implicit none
      include "py6strf.inc"
      include "mc_set.inc"

        double precision dQ2, dx, dF2, dF1, DF2NF2P
        DOUBLE PRECISION F2ALLM, DNP(5), DR, XPQ
        DOUBLE PRECISION gamma2, dnu, w2, pmass2
        integer A, Z
        DIMENSION XPQ(-25:25)

C parameters for ratio F2(n)/F2(p)
C measured at NMC, (NMC Amaudruz et al. CERN-PPE/91-167))
        data DNP   / 0.976D0,    -1.34D0,      1.319D0,
     &              -2.133D0,     1.533D0/
 
        pmass2=massp**2
        w2=pmass2+(dq2*(1/dx-1))
        dnu=(w2-pmass2+dq2)/(2.*massp)
        gamma2=dq2/(dnu**2)
    
        if ((A.eq.1).and.(Z.eq.1)) then
*:        ALLM:
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
            Call MKR(DQ2,DX,DR)
            DF2=F2ALLM(dx,dq2)
            DF1=(1.D0+gamma2)/(2.D0*dx)/(1.D0+DR)*DF2
         ELSEIF(genSet_FStruct(1:4).EQ.'F2PY') THEN
            call F2PYTH(dx,dq2,df1,df2,z)
         ELSEIF(genSet_FStruct(1:4).EQ.'PDF') THEN
            call PYPDFU(2212,dX,dQ2,XPQ)
            DF2=1D0/9D0*(XPQ(1)+XPQ(-1)+XPQ(3)+XPQ(-3))+
     &          4D0/9D0*(XPQ(2)+XPQ(-2))
            DF1=(1.D0+gamma2)/(2.D0*dx)*DF2
         ELSE
*:           error:
            write(*,*)('invalid parametrisation choice in mkf2')
         ENDIF
      ELSEIF(A.EQ.2.and.z.eq.1)THEN
*:        ALLM:
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
            Call MKR(DQ2,DX,DR)
            DF2=F2ALLM(dx,dq2)
            DF2NF2P=DNP(1)+dx*(DNP(2)+dx*(DNP(3)+dx*(DNP(4)+dx*DNP(5))))
            DF2=DF2*0.5*(df2nf2p+1.)
            DF1=(1.D0+gamma2)/(2.D0*dx)/(1.D0+DR)*DF2
         ELSE
*:           error
            write(*,*)('MKF2: invalid parametrisation choice FStruct')
         ENDIF
*
! ... neutron = 2*(deuterium_per_nucleon) - proton:
      ELSEIF(A.EQ.1.and.z.eq.0)THEN
*:        ALLM:
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
            DF2=F2ALLM(dx,dq2)
            Call MKR(DQ2,DX,DR)
            DF2NF2P=DNP(1)+dx*(DNP(2)+dx*(DNP(3)+dx*(DNP(4)+dx*DNP(5))))
            DF2=DF2*df2nf2p
            DF1=(1.D0+gamma2)/(2.D0*dx)/(1.D0+DR)*DF2
         ELSEIF(genSet_FStruct(1:4).EQ.'F2PY') THEN
            call F2PYTH(dx,dq2,df1,df2,z)
         ELSEIF(genSet_FStruct(1:4).EQ.'PDF') THEN
            call PYPDFU(2112,dX,dQ2,XPQ)
            DF2=1D0/9D0*(XPQ(1)+XPQ(-1)+XPQ(3)+XPQ(-3))+
     &          4D0/9D0*(XPQ(2)+XPQ(-2))
            DF1=(1.D0+gamma2)/(2.D0*dx)*DF2
         ELSE
*:           error
            write(*,*)('MKF2: invalid parametrisation choice FStruct')
         ENDIF
      ELSEIF((A.gt.1).and.(Z.gt.1)) then
         IF(genSet_FStruct(1:4).EQ.'PDF') THEN
            call PYPDFU(2212,dX,dQ2,XPQ)
            DF2=1D0/9D0*(XPQ(1)+XPQ(-1)+XPQ(3)+XPQ(-3))+
     &          4D0/9D0*(XPQ(2)+XPQ(-2))
            DF1=(1.D0+gamma2)/(2.D0*dx)*DF2
         ELSE
*:           error:
            write(*,*)('invalid parametrisation choice in mkf2')
         ENDIF
      ELSE
*:        error:
         write(*,*)('MKF2: invalid target type'), A, Z
      ENDIF

       return
       end

!---------------------------------------------------------------
! Calculate R = sigma_L/sigma_T

      SUBROUTINE MKR(DQ2,DX,DR)
      IMPLICIT NONE
      include "mc_set.inc"
      include "py6strf.inc"

      DOUBLE PRECISION DQ2, DX
      DOUBLE PRECISION DR,DELTAR

      IF ( genSet_R .EQ. '1990' ) THEN
*        Whitlow et al.,  Phys.Lett.B 250(1990),193
         CALL R1990(DQ2,DX,DR,DELTAR)
         py6R=DR
      ELSE IF ( genSet_R .EQ. '1998' ) THEN
*        E143, hep-ex/9808028
         CALL R1998(DQ2,DX,DR,DELTAR)
         py6R=DR
      ELSE IF ( genSet_R .eq. '0' ) THEN
* pure transverse (sigma_L=0)
         DR = 0.d0
         py6R=0.d0
      ELSE
         write(*,*)( 'MKR: invalid choice for R parametrization' )
      ENDIF
     
      RETURN
      END

C------------------------------------------------------------------

      SUBROUTINE R1990(DQ2,DX,DR,DELTAR)

      IMPLICIT NONE

      DOUBLE PRECISION DQ2, DX
      DOUBLE PRECISION DR, DELTAR

      REAL R
      REAL QQ35, XX
      REAL FAC, RLOG, Q2THR
      REAL R_A, R_B, R_C
C
C Data-Definition of R-Calculation, see
C            L.W.WHITLOW, SLAC-REPORT-357,
C            PH.D. THESIS, STANFORD UNIVERSITY,
C            MARCH 1990.
      REAL AR1990(3), BR1990(3), CR1990(3)
      DATA AR1990  / .06723, .46714, 1.89794 /
      DATA BR1990  / .06347, .57468, -.35342 /
      DATA CR1990  / .05992, .50885, 2.10807 /

      DELTAR = 0.

      XX=real(DX)
      IF (DQ2.LT.0.35) THEN
        QQ35=0.35
      ELSE
        QQ35=real(DQ2)
      ENDIF
C
C *** If Q2 < 0.35 then variable "R" is calculated at the fixed Q2 of 0.35
C
      FAC   = 1+12.*(QQ35/(1.+QQ35))*(.125**2/(XX**2+.125**2))
      RLOG  = FAC/LOG(QQ35/.04)
      Q2THR = 5.*(1.-XX)**5

      R_A   = AR1990(1)*RLOG +
     &        AR1990(2)/SQRT(SQRT(QQ35**4+AR1990(3)**4))
      R_B   = BR1990(1)*RLOG +
     &        BR1990(2)/QQ35 + BR1990(3)/(QQ35**2+.3**2)
      R_C   = CR1990(1)*RLOG +
     &        CR1990(2)/SQRT((QQ35-Q2THR)**2+CR1990(3)**2)
      R     = (R_A+R_B+R_C)/3.

      IF (DQ2.GE.0.35) THEN
        DR=dble(R)
      ELSE
        DR=dble(R)*DQ2/0.35
      ENDIF

c      print*,'R:',R
      
      END

C-----------------------------------------------------------------------
      SUBROUTINE R1998(DQ2,DX,DR,DELTAR)

C new fit to R  hep-ex/9808028 E143 Collab.
C it is based on the old 3 paramter forms
C 0.005<x<0.86, 0.5<Q2<130 GeV2
C E143 x-section measurement 0.03<x<0.4
C with  overall norm uncertainty 2.5%
C
C U. Stoesslein, October 1998
C

      IMPLICIT NONE

      DOUBLE PRECISION DQ2,DX,DR,DELTAR
      DOUBLE PRECISION Q2,Q2MAX,FAC,RLOG,Q2THR
      DOUBLE PRECISION R_A_NEW,R_A,R_B_NEW,R_B,R_C

      DOUBLE PRECISION A(6),B(6),C(6)

      DATA A/ .0485, 0.5470, 2.0621, -.3804,  0.5090, -.0285/
      DATA B/ .0481, 0.6114, -.3509, -.4611,  0.7172, -.0317/
      DATA C/ .0577, 0.4644, 1.8288,12.3708,-43.1043,41.7415/ 

      DOUBLE PRECISION DR1998
      EXTERNAL DR1998

* use R(0.35) if q2 is below 0.35
      Q2=DQ2
      Q2MAX=0.35
      IF(Q2.LT.Q2MAX) Q2=Q2MAX

      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(DX**2+.125**2))
      RLOG  = FAC/LOG(Q2/.04)
      Q2thr = C(4)*DX+C(5)*DX*DX+C(6)*DX*DX*DX

* new additional terms
      R_A_NEW = (1.+A(4)*DX+A(5)*DX*DX)*DX**(A(6))
      R_A   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))*R_A_NEW
      R_B_NEW = (1.+B(4)*DX+B(5)*DX*DX)*DX**(B(6))
      R_B   = B(1)*RLOG + (B(2)/Q2 + B(3)/(Q2**2+0.3**2))*R_B_NEW
      R_C   = C(1)*RLOG + C(2)/SQRT((Q2-Q2thr)**2+C(3)**2)
      DR    = (R_A+R_B+R_C)/3.

* straight line fit extrapolation to R(Q2=0)=0
      if (dq2.lt.q2max) DR = DR*DQ2/Q2MAX

* I assume the fit uncertainty only for measured Q2 range
      if (Q2 .GT. 0.5) then
         DELTAR = DR1998(DX,Q2)
      else
         DELTAR=DR
      endif

      RETURN
      END

C--------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DR1998(DX,DQ2)

* Parameterize uncertainty in R1998 
* associated to the fitting procedure only

      IMPLICIT NONE
      DOUBLE PRECISION DX,DQ2

      DR1998 = 0.0078-0.013*DX+(0.070-0.39*DX+0.7*DX*DX)/(1.7+DQ2)

      RETURN
      END

!---------------------------------------------------------------
! Calculate the asymmetries A1 and A2. Routine is empty because Pythia
! is unpolarised, but radgen expects it

       subroutine mkasym (dQ2, dX, A, Z, dA1, dA2)

       implicit none

       double precision dQ2, dx, dA1, dA2
       integer A, Z

        da1 = 0.D0
        da2 = 0.D0

       return
       end

!---------------------------------------------------------------
! Calculate the dilution factor

       double precision function fdilut (dQ2, dx, A)

       implicit none

       DOUBLE PRECISION DQ2, DX, DF
       DOUBLE PRECISION DNP, DFN, DFP
       DOUBLE PRECISION DZ, DF2NF2P
       INTEGER A
       dimension dnp(7)
*       
C ... fit to NMC  F2n/F2p data (86/87+89 T1,T14)
       data dnp/
     +   0.67225D+00,
     +   0.16254D+01,
     +  -0.15436D+00,
     +   0.48301D-01,
     +        0.41979D+00,
     +   0.47331D-01,
     +  -0.17816D+00/

! Definitions of 'intrinsic' dilutions for neutron and proton (GNOME confused)

       dfn=1.
       dfp=1.

! Only 3He has a dilution, and we determine it as F2n/(2*F2p + F2n)

       if (A.ne.3) then
         df = 1.
       else
         dZ = 1./2.*DLOG(1.+DEXP(2.0-1000.*dX))
         df2nf2p = dnp(1)*(1.0-dx)**dnp(2)+dnp(3)*dx**dnp(4)
     1            +(dnp(5)+dnp(6)*dz+dnp(7)*dz**2)
         df = dfn*(1./((2./df2nf2p)+1))
       endif

       FDILUT = DF
       return

       END

C MDB 2018-06-17 Remove for BeAGLE. Identical routine in RAPGAP library
C                taken from Heracles466.
C!---------------------------------------------------------------
C! Function to calculate F2 from ALLM parametrisation
C
C      double precision FUNCTION f2allm(x,q2)
C
C      implicit double precision (a-h,o-z)
C
C      REAL M02,M12,LAM2,M22
C      COMMON/ALLM/SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
CC  POMERON
C      PARAMETER (
C     , S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
C     , A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
C     , B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
C     , M12   =  49.457 )
C 
CC  REGGEON
C     PARAMETER (
C     , S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
C     , A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
C     , B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
C     , M22   =   0.15052 )
CC
C      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )
C      PARAMETER ( ALFA=112.2, XMP2=0.8802)
CC                                                                             C  
C      W2=q2*(1./x -1.)+xmp2
C      W=dsqrt(w2)
CC
C      IF(Q2.EQ.0.) THEN
C       S=0.
C       Z=1.
CC                                                                             C  
CC   POMERON                                                                   C  
CC                                                                             C  
C       XP=1./(1.+(W2-XMP2)/(Q2+M12))
C       AP=A11
C       BP=B11 
C       SP=S11 
C       F2P=SP*XP**AP 
CC                                                                             C  
CC   REGGEON                                                                   C  
CC                                                                             C  
C       XR=1./(1.+(W2-XMP2)/(Q2+M22)) 
C       AR=A21 
C       BR=B21 
C       SR=S21
C       F2R=SR*XR**AR 
CC                                                                             C  
C      ELSE
C       S=DLOG(DLOG((Q2+Q02)/LAM2)/DLOG(Q02/LAM2)) 
C       Z=1.-X 
CC                                                                             C  
CC   POMERON                                                                   C  
CC                                                                             C  
C       XP=1./(1.+(W2-XMP2)/(Q2+M12)) 
C       AP=A11+(A11-A12)*(1./(1.+S**A13)-1.) 
C       BP=B11+B12*S**B13 
C       SP=S11+(S11-S12)*(1./(1.+S**S13)-1.) 
C       F2P=SP*XP**AP*Z**BP 
CC                                                                             C  
CC   REGGEON                                                                   C  
CC                                                                             C  
C       XR=1./(1.+(W2-XMP2)/(Q2+M22)) 
C       AR=A21+A22*S**A23 
C       BR=B21+B22*S**B23  
C       SR=S21+S22*S**S23 
C       F2R=SR*XR**AR*Z**BR 
C    
CC                                                                             C  
C      ENDIF
C
Cc      CIN=ALFA/(Q2+M02)*(1.+4.*XMP2*Q2/(Q2+W2-XMP2)**2)/Z 
Cc      SIGal=CIN*(F2P+F2R) 
Cc      f2allm=sigal/alfa*(q2**2*(1.-x))/(q2+4.*xmp2*x**2)
C      f2allm = q2/(q2+m02)*(F2P+F2R)
C
C      RETURN  
C      END 
CC--------------------------------------------------------------------

C.....Total inclusive Pythia cross section according to the model in  
C.....C.Friberg, T.Sjoestrand, J. High Energy Phys. JHEP 0009, 010 (2000)
C.....(neglecting GVMD cross sections)
C.....P.Liebing, 08/16/2003

      SUBROUTINE F2PYTH(x,q2,f1,f2,z)

      IMPLICIT NONE
      COMMON/PYINT1/MINT(400),VINT(400)
      INTEGER MINT
      DOUBLE PRECISION VINT
      SAVE/PYINT1/
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      INTEGER MSTP,MSTI
      DOUBLE PRECISION PARP,PARI
      SAVE/PYPARS/
      DOUBLE PRECISION f2dis,f1dis,sdis,svmd1,svmd2,sigvm,x,q2,dipol
      DOUBLE PRECISION mrho2,alpha,eps,eta,pmass,pmass2,gevmb,rvmd
      DOUBLE PRECISION w2,gamma2,nu,convf2,convf1,conv,sigh,pi
      DOUBLE PRECISION f1,f2,df2allm,df2allml,df1allm,df1allml,r,rl
      DOUBLE PRECISION w2l,gamma2l,nul,xl,sw2
      INTEGER z
C...Local arrays.
      DOUBLE PRECISION XPQ,XPAR,YPAR
      INTEGER i
      DIMENSION XPQ(-25:25),XPAR(4),YPAR(4)
     
C...X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
      DATA XPAR/2*13.63D0,10.01D0,0.970D0/
      DATA YPAR/2*31.79D0,-1.51D0,-0.146D0/
      DATA EPS/0.0808D0/,ETA/-0.4525D0/
      DATA MRHO2/0.591822D0/,ALPHA/7.297352533D-3/,PMASS/0.93827D0/,
     &     GEVMB/0.389379292D0/,pi/3.14159265358979324D0/

      EXTERNAL PYPDFU 
      f2dis=0D0
      f1dis=0D0
      sdis=0D0
      svmd1=0D0
      svmd2=0D0
      sigvm=0D0
      rvmd=0D0
      dipol=0D0
      w2=0D0
      gamma2=0D0
      nu=0D0
      convf2=0D0
      convf1=0D0
      conv=0D0
      sigh=0D0
      f1=0D0
      f2=0D0
      df1allm=0D0
      df2allm=0D0
      df1allml=0D0
      df2allml=0D0
      r=0D0
      rl=0D0
      sw2=0D0
c....Kinematic factors and constants
      pmass2=pmass**2
      nul=-1D0
      gamma2l=0D0
      w2l=-1D0
      xl=0D0
      
      if ((x.gt.0D0).and.(x.le.1D0)) then
       w2=pmass2+(q2*(1./x-1.))
      else
       f1=0D0
       f2=0D0
       return
      endif
      
      if (w2.lt.4D0) then
       w2l=w2
       w2=4D0
       xl=x
       nul=(w2l-pmass2+q2)/(2D0*pmass)
      endif
      
      nu=(w2-pmass2+q2)/(2D0*pmass)

      if (nu.gt.0D0) then
       gamma2=q2/(nu**2)
       if (nul.gt.0D0) then
        gamma2l=q2/(nul**2)
       endif
      else
       f1=0D0
       f2=0D0
       return
      endif

c....For W2<4, we don't trust the PYTHIA F2, so we calculate 
C....F2-ALLM(W2,Q2). The real kinematics have and "l" at the end, 
C....the kinematics without the "l" are the ones we get by setting W2 to 4
c....Output: f2allm=F2-ALLM(W2=4,Q2),f2allml=F2-ALLM(W2=w2l,Q2)       
      if (w2l.gt.0D0) then
       sw2=((w2l-pmass2)/(4D0-pmass2))**10
      endif 
c....This factor is needed to convert the Pythia virtual gamma cross 
c....section for VMD to the same level as F2
c....The kinematic factors making the (ep) cross section out of F2 are
c....provided by RADGEN      
      conv=q2*(1D0-x)/(4D0*pi**2*alpha)/gevmb
c....Pythia PDF call, sum PDFs to F2
      if (z.eq.1) then
      call PYPDFU(2212,X,Q2,XPQ)
      elseif (z.eq.0) then
      call PYPDFU(2112,X,Q2,XPQ)
      endif
       f2dis=1D0/9D0*(XPQ(1)+XPQ(-1)+XPQ(3)+XPQ(-3))+
     &       4D0/9D0*(XPQ(2)+XPQ(-2))
c....Suppression factor for DIS
      if (MSTP(19).eq.0) then
       sdis=1.
      else  
       sdis=q2/(q2+mrho2)
       if (MSTP(19).gt.1) then
        sdis=sdis**2
       endif
      endif
C....Sum of Hadronic (Vector Meson) cross sections * Photon couplings
C....const.
      sigh=0.
      do 10 i=1,4
       sigh=sigh+alpha/PARP(160+i)*(XPAR(i)*w2**eps+YPAR(i)*w2**eta) 
   10 continue
C....W2/Q2 suppression of VMD and (1+epsilon R_VMD)
      svmd1=(w2/(w2+q2))**MSTP(20)
      if (MSTP(20).eq.0) then
       dipol=2.575D0
      else
       dipol=2D0
      endif
      if (MSTP(17).eq.6) then
       rvmd=PARP(165)*(q2/mrho2)**PARP(166)
      else
c    ...Attention: This is only good for MSTP(17)=4, i.e., the Pythia
c    ...default       
       rvmd=PARP(165)*(4.*mrho2*q2)/(mrho2+q2)**2
      endif
C  .... Dipole factor for VMD      
      svmd2=(mrho2/(mrho2+q2))**dipol
C.....virtual photon xsec for VMD
      sigvm=svmd1*svmd2*sigh
      convf2=(1D0+rvmd)/(1D0+gamma2)
      convf1=1D0/(2D0*x)
c.....Total "F2"
      f2=sdis*f2dis+conv*convf2*sigvm
      f1dis=(1.D0+gamma2)/(2.D0*x)*f2dis
      f1=sdis*f1dis+conv*convf1*sigvm
      if (w2l.gt.0D0) then
C.....Here we scale F2-ALLM(W2=w2l,Q2) by the factor 
C.....F2-PYTH(W2=4,Q2)/F2-ALLM(W2=4,Q2) (normalize ALLM to PYTHIA model)
        f2=sw2*f2
        f1=sw2*f1
      endif
      RETURN
      END
