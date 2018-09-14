************************************************************************
*
* II. Nuclear density (fm^-3) v.2.0
      function nucdens(r,Z,A,idist,irho)
*     programmer: Alberto Accardi
*     date: 17/03/02 & 14 Nov 03 (new version)
*
*  A. COMMENTARY
*
*     Returns the nuclear density normalized to 1, according to the
*     value of "thickfn". After each call, in the common block /nukrad/
*     one can find the "hard-sphere equivalent" radius
*       RReq = sqrt(5/3) * sqrt(<r^2>)
*     and when a Wood-Saxon distribution is selected also the radius
*     parameter RRws of the WS parametrization.
*
*     r     = (fm) distance from the nuclear center
*     Z     = Atomic number
*     A     = Atomic mass
*     idist = 0 --> Hard Sphere (HS)
*             1 --> Woods-saxon (WS)
*     irho  = user defined number of the nuclear density: it allows
*             to use several densities at the same time (important
*             when using WS distributions).
*
*     NOTES ON USAGE:
*     --------------
*     + if Z=1,A=2 (Deuterium) the Reid's soft core potential [5]
*       is always used
*     + Woods-Saxon (WS) distributions
*       - if Z>0 uses experimental values for 2- or 3-parameters WS when
*         available, otherwise a 2-parameter WS is used with a fitted
*         parameters as a function of A, see below. NOTE that only the
*         most common isotope at a given Z is considered.
*       - if Z=0 then uses the fitted parameters always.
*
*     HARD-SPHERE DISTRIBUTIONS:
*     --------------------------
*     A uniform distribution in a sphere of radius RA = 1.12 * A^(1/3)
*     This is the value stored in the common block /nukrad/
*
*     WOODS-SAXON DISTRIBUTIONS:
*     -------------------------
*     The 3-parameters Woods-Saxon distribution rho is:
*                       1 + c0*(r^2/RA^2)
*       rho(r) = rho0 -------------------
*                      1 + exp((r-RA)/a0)
*     which reduces to the standard 2-parameter WS when c0=0. When
*     experimental data for the parameters RA,a0,c0 are available from
*     Refs.[1,2] we use them, otherwise we use the following values [4]:
*       RA = (0.978d0+0.0206d0*(A**1/3))*(A**1/3)
*       a0 = 0.523
*       c0 = 0
*     This we prefer to the one used in [3], namely
*       RR = 1.19d0*(A**1/3) - 1.61d0*(A**(-1/3))
*       a0 = 0.54
*       c0 = 0
*     as we checked that it reproduces more correctly the results for nuclei
*     with known experimental values. If wanted, one can uncomment the latter
*     in the code below and comment the Bialas parametrization.
*
*     Here are the policies we followed in choosing the experimental values:
*     1) We considered for each Z only the most common isotope.
*     2) Ref.[2] preferred if values available in both Ref.[1] and [2]
*        (Actually, ref.[1] used only for 4He and 32S).
*     3) We chose 3-params WS instead of 2-params WS if both were available.
*     4) If more than one set of values available, we picked up that with
*        the highest momentum transfer (q) coverage.
*
*     During the first call, the routine computes the normalization
*     for the WS distributions and stores it in the Nws(irho) array.
*     It is then used in the following calls. 
*     ******* This is to be improved in the future: just compute 
*     ******* it once for all and store it as a parameter.
*     ******* The usage would be simplified quite a lot
*     ******* because the 'irho' variable would no more be necessary
*
*     3-parameters distributions may become negative at large r, because
*     often c0<0. When this happens the routine sets rho=0.
*
*     The nuclear radius stored in the common blok /nukrad/ is
*     computed as RA^2 = (3/5) <r^2>. This is the "Hard-sphere equivalent 
*     radius" (the definition reproduces identically RA of a hard-spere 
*     distribution).  
*
*     REFERENCES:
*     -----------
*     [1] Landolt-Boernstein, "Numerical data and functional relationships 
*         in science and technology. Vol.2, nuclear radii", Springer-Verlag 
*         ~1969
*     [2] H.De Vries,C.W.De Jager, C.De Vries, Atomic Data and Nuclear Data 
*         Tables, 36(87)495-536
*     [3] K.J.Eskola, X.N.Wang, Int.J.Mod.Phys.A10(95)3087-3090
*     [4] A.Bialas, Phys.Lett.B133(83)241
*     [5] R.V.Reid Jr., Ann.Phys.50(68)411-448
*
*  B. DECLARATIONS
*
      implicit none

*     MDB Allow BeAGLE to override the Woods-Saxon parameters
      include 'beagle.inc'

*    *** variables

      integer Z,A,idist,irho

      double precision nucdens,r

*   *** functions

      double precision dgauss11,DeuteronDensity,WoodsSaxon3d
     &     , r2DD, r2WS3d
      external WoodsSaxon3d, r2DD, r2WS3d

*    *** common blocks

      double precision RR,a0,c0
      common/ woodssaxon/ RR,a0,c0

      double precision RAeq,RAws
      common/nukrad/RAeq,RAws

*    *** parameters and initial values

      double precision onethird, srft
      parameter (onethird=0.333333333333333d0,srft=1.29099)

      double precision twopi, threefourthoverpi
      parameter (twopi=6.2831530718d0
     &     , threefourthoverpi=0.2387336394417d0)

      integer Asave(10),Zsave(10)
      logical firsttime(10)
      double precision Rws(92),aws(92),cws(92) 
      double precision NWS(10),RA(10),RAequiv(10),aa(10),cc(10)
      save firsttime,NWS,RA,RAequiv,aa,cc,Asave,Zsave

*    ... first time flag fr normalization of WS
      data firsttime /.true.,.true.,.true.,.true.,.true.
     &     ,.true.,.true.,.true.,.true.,.true./
      
*    ... Woods-Saxon params (as a function of the atomic no. Z)
      data Rws(2)  , aws(2) , cws(2)   / 
     :     1.01d0  , 0.327d0, 0.445d0  / 
      data Rws(3)  , aws(3) , cws(3)   / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(4)  , aws(4) , cws(4)   / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(5)  , aws(5) , cws(5)   / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(6)  , aws(6) , cws(6)   / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(7)  , aws(7) , cws(7)   / 
     :     2.570d0 , 0.505d0, -0.180d0 / 
      data Rws(8)  , aws(8) , cws(8)   / 
     :     2.608   , 0.513  , -0.051d0 / 
      data Rws(9)  , aws(9) , cws(9)   / 
     :     2.58d0  , 0.567d0, 0d0      /
      data Rws(10) , aws(10), cws(10)  / 
     :     2.791d0 , 0.698d0, -0.168d0 /
      data Rws(11) , aws(11), cws(11)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(12) , aws(12), cws(12)  / 
     :     3.192d0 , 0.604d0, -0.249   /
      data Rws(13) , aws(13), cws(13)  / 
     :     3.07d0  , 0.519d0, 0d0      /
      data Rws(14) , aws(14), cws(14)  / 
     :     3.340d0 , 0.580d0, -0.233d0 /
      data Rws(15) , aws(15), cws(15)  / 
     :     3.369d0 , 0.582d0, -0.173d0 /
      data Rws(16) , aws(16), cws(16)  / 
     :     3.20d0  , 0.59d0 , 0d0      /
      data Rws(17) , aws(17), cws(17)  / 
     :     3.476d0 , 0.599d0, -0.10d0  /
      data Rws(18) , aws(18), cws(18)  / 
     :     3.73d0  , 0.63d0 , -0.19d0  /
      data Rws(19) , aws(19), cws(19)  / 
     :     3.743d0 , 0.585d0, -0.201d0 /
      data Rws(20) , aws(20), cws(20)  / 
     :     3.766d0 , 0.586d0, -0.161d0 /
      data Rws(21) , aws(21), cws(21)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(22) , aws(22), cws(22)  / 
     :     3.843d0 , 0.588d0, 0d0      /
      data Rws(23) , aws(23), cws(23)  / 
     :     3.94d0  , 0.505d0, 0d0      /
      data Rws(24) , aws(24), cws(24)  / 
     :     3.98d0  , 0.542d0, 0d0      /
      data Rws(25) , aws(25), cws(25)  / 
     :     3.89d0  , 0.567d0, 0d0      /
      data Rws(26) , aws(26), cws(26)  / 
     :     4.111d0 , 0.558d0, 0d0      /
      data Rws(27) , aws(27), cws(27)  / 
     :     4.158d0 , 0.575d0, 0d0      /
      data Rws(28) , aws(28), cws(28)  / 
     :     4.309d0 , 0.517d0, -0.131d0 /
      data Rws(29) , aws(29), cws(29)  / 
     :     4.218d0 , 0.596d0, 0d0      /
      data Rws(30) , aws(30), cws(30)  / 
     :     4.285   , 0.584  , 0d0      /
      data Rws(31) , aws(31), cws(31)  / 
     :     0d0     , 0d0    , 0d0      / 
      data Rws(32) , aws(32), cws(32)  / 
     :     4.45d0  , 0.573d0, 0d0      /
      data Rws(33) , aws(33), cws(33)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(34) , aws(34), cws(34)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(35) , aws(35), cws(35)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(36) , aws(36), cws(36)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(37) , aws(37), cws(37)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(38) , aws(38), cws(38)  / 
     :     4.83d0  , 0.496d0, 0d0      /
      data Rws(39) , aws(39), cws(39)  / 
     :     4.76d0  , 0.571d0, 0d0      /
      data Rws(40) , aws(40), cws(40)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(41) , aws(41), cws(41)  / 
     :     4.87d0  , 0.573d0, 0d0      /
      data Rws(42) , aws(42), cws(42)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(43) , aws(43), cws(43)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(44) , aws(44), cws(44)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(45) , aws(45), cws(45)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(46) , aws(46), cws(46)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(47) , aws(47), cws(47)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(48) , aws(48), cws(48)  / 
     :     5.38d0  , 0.532d0, 0d0      /
      data Rws(49) , aws(49), cws(49)  / 
     :     5.357d0 , 0.563d0, 0d0      /
      data Rws(50) , aws(50), cws(50)  / 
     :     5.442d0 , 0.543d0, 0d0      /
      data Rws(51) , aws(51), cws(51)  / 
     :     5.32d0  , 0.57d0 , 0d0      /
      data Rws(52) , aws(52), cws(52)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(53) , aws(53), cws(53)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(54) , aws(54), cws(54)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(55) , aws(55), cws(55)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(56) , aws(56), cws(56)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(57) , aws(57), cws(57)  / 
     :     5.71d0  , 0.535d0, 0d0      /
      data Rws(58) , aws(58), cws(58)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(59) , aws(59), cws(59)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(60) , aws(60), cws(60)  / 
     :     5.626d0 , 0.618d0, 0d0      /
      data Rws(61) , aws(61), cws(61)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(62) , aws(62), cws(62)  / 
     :     5.804d0 , 0.581d0, 0d0      /
      data Rws(63) , aws(63), cws(63)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(64) , aws(64), cws(64)  / 
     :     5.930d0 , 0.576d0, 0d0      /
      data Rws(65) , aws(65), cws(65)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(66) , aws(66), cws(66)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(67) , aws(67), cws(67)  / 
     :     6.18d0  , 0.57d0 , 0d0      /
      data Rws(68) , aws(68), cws(68)  / 
     :     5.98d0  , 0.446  , 0.19d0   /
      data Rws(69) , aws(69), cws(69)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(70) , aws(70), cws(70)  / 
     :     6.127d0 , 0.363d0, 0d0      /
      data Rws(71) , aws(71), cws(71)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(72) , aws(72), cws(72)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(73) , aws(73), cws(73)  / 
     :     6.38d0  , 0.64d0 , 0d0      /
      data Rws(74) , aws(74), cws(74)  / 
     :     6.51d0  , 0.535d0, 0d0      /
      data Rws(75) , aws(75), cws(75)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(76) , aws(76), cws(76)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(77) , aws(77), cws(77)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(78) , aws(78), cws(78)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(79) , aws(79), cws(79)  / 
     :     6.38d0  , 0.535d0, 0d0      /
      data Rws(80) , aws(80), cws(80)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(81) , aws(81), cws(81)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(82) , aws(82), cws(82)  / 
     :     6.62d0  , 0.546d0, 0d0      /
      data Rws(83) , aws(83), cws(83)  / 
     :     6.75d0  , 0.468d0, 0d0      /
      data Rws(84) , aws(84), cws(84)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(85) , aws(85), cws(85)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(86) , aws(86), cws(86)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(87) , aws(87), cws(87)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(88) , aws(88), cws(88)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(89) , aws(89), cws(89)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(90) , aws(90), cws(90)  / 
     :     6.792d0 , 0.571d0, 0d0      /
      data Rws(91) , aws(91), cws(91)  / 
     :     0d0     , 0d0    , 0d0      /
      data Rws(92) , aws(92), cws(92)  / 
     :     6.86d0  , 0.42d0 , 0d0      /
C MDB kludge - replace original with Tang U-HN above
C     :     6.874d0 , 0.556d0, 0d0      /


*
*  C. ACTION   
*
      if(.not. firsttime(irho)) then
        if ((Asave(irho).ne.A).or.(Zsave(irho).ne.Z)) then
          firsttime(irho)=.true.
        end if
      end if
      Asave(irho)=A
      Zsave(irho)=Z
      
      if ((Z.eq.1).and.(A.eq.2).and.(idist.ne.1)) then
*       *** REID's SOFT-CORE 
         if (firsttime(irho)) then
            RAequiv(irho) = srft 
     &           * sqrt(dgauss11(r2DD,0d0,20d0,1d-5))
            firsttime(irho) = .false.
         end if
         nucdens = DeuteronDensity(r)
         RAeq = RAequiv(irho)
         RAws = 1.4111d0
      else if (idist.eq.0) then
*       *** HARD SPHERE
         RAeq = 1.12d0*(A**onethird) 
         RAws = RAeq
         if (r.le.RAeq) then
            nucdens = threefourthoverpi / (RAeq**3) 
         else
            nucdens = 0
         end if
      else if (idist.eq.1) then
*       *** WOODS - SAXON  
         if (firsttime(irho)) then
            if (USER3D) then
*              BeAGLE override of Woods-Saxon (Fermi) parameters
               RR = RG3D
               a0 = AG3D
               c0 = WG3D
               WRITE(*,*) "NUCDENS: BeAGLE override for RR,a0,c0:",
     &              RR, ", ", a0, ", ", c0
            elseif ((Z.eq.0).or.(Z.gt.92).or.(Rws(Z).lt.0.1d0)) then
*             ... parametrization from Eskola,Wang, Ref.[3]
*               RR = 1.19d0*(A**onethird) - 1.61d0*(A**(-onethird))
*               a0 = 0.54d0
*               c0 = 0d0
*             ... parametrization from Bialas, Ref.[4]
               RR = (0.978d0+0.0206d0*(A**onethird))*(A**onethird)
               a0 = 0.523d0
               c0 = 0d0
            else
               RR = Rws(Z)
               a0 = aws(Z)
               c0 = cws(Z)
            end if
            NWS(irho) = 1d0/dgauss11(WoodsSaxon3d,0d0,10d0*RR,1d-5)
            RAequiv(irho) = srft * dsqrt( NWS(irho) 
     &           * dgauss11(r2WS3d,0d0,10d0*RR,1d-5))
            RA(irho) = RR
            aa(irho) = a0
            cc(irho) = c0
            firsttime(irho) = .false.
         end if
         RAeq = RAequiv(irho)
         RAws = RA(irho)
         if (r-RA(irho).lt.700d0) then
            nucdens = NWS(irho) 
     &           * (1d0+cc(irho)*(r*r/(RA(irho)*RA(irho))))
     &           / (1d0+dexp((r-RA(irho))/aa(irho)))

         else 
            nucdens = 0d0
         end if
         if (nucdens.lt.0d0) nucdens = 0d0
      else 
         print*, 'ERROR (nucdens): called out of range:'
     :        , 'idist,Z,A=', idist,Z,A
         stop
      end if
      
      return
      end


************************************************************************
*
*VI.  3D (unnormalized) Wodds-Saxon distribution 
      function WoodsSaxon3d(r)
*     programmer: Alberto Accardi
*     date: 10/05/00
*
*  A. COMMENTARY
*
*     It returns the value of the Nuclear thickness function (GeV^-2)
*     times the Jacobian for d^3r = 4*pi * r^2 * dr
*     where:
*     bx, by = transverse coordinates (origin in the centre of the
*              nucleus
*     A = atomic number
*
*  B. DECLARATIONS
*
      implicit none
 
*    *** variables

      double precision WoodsSaxon3d, r

      double precision pi
      parameter(pi=3.141592653d0)

*   *** common blocks

      double precision RA,aa,cc
      common/ woodssaxon/ RA,aa,cc

*
*  C. ACTION
*

      WoodsSaxon3d = 4*pi*r**2 *(1+cc*(r*r/(RA*RA)))/(1+dexp((r-RA)/aa))
      if (WoodsSaxon3d.lt.0d0) WoodsSaxon3d = 0d0 

      return
      end

************************************************************************
*     integrand for <r^2>
      double precision function r2WS3d(r)
*     programmer: Alberto Accardi
*     date: 20/07/04
*
*  A. COMMENTARY
*
*     It returns the value of the r^2 * Wood-Saxon nuclear thickness function 
*     (GeV^-2) at a radius r [fm]
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables

      double precision r

*    *** functions

      double precision WoodsSaxon3d

*
*  C. ACTION
*

      r2WS3D = r*r * WoodsSaxon3d(r)
*      r2WS3D = 4*3.14* r**4  

      return
      end


************************************************************************
      function DeuteronDensity(r)
*     programmer: Daniel Gruenewald
*     date: 23/02/04
*
*  A. Commentary
*     Returns the value of the deuteron thickness function.
*     r       =(fm) distance from the origin
*     Based on the Reid's soft core Deuteron wave function
*
*  B. Declarations
      double precision DeuteronDensity, reidsc, r, twopi
      parameter (twopi=6.2831530718d0)
      
*  C. Action
      DeuteronDensity = reidsc(2*r)/(twopi*r*r)
      return 
      end
      

************************************************************************
*     integrand for <r^2>
      double precision function r2DD(r)
*     programmer: Alberto Accardi
*     date: 20 July 2004
*
*  A. Commentary
*     Returns r^2 times J=4*pi*r^2 times
*     the value of the deuteron thickness function.
*     r       =(fm) distance from the origin
*     Based on the Reid's soft core Deuteron wave function
*
*  B. Declarations
      double precision reidsc, r, twopi
      parameter (twopi=6.2831530718d0)
      
*  C. Action
      r2DD = 2*r*r*reidsc(2*r)
      return 
      end
      

************************************************************************
*
*VI.  Reid's soft core deuteron wave function [fm^-3]
      function reidsc(r)
*     programmer: Alberto Accardi
*     date: 26/03/02
*
*  A. COMMENTARY
*
*     Returns the Reid's soft core wave-function of the deuteron squared,
*     i.e. the deuteron electron density normalized to 1.
*
*     r       = (fm) distance from the origin
*
*     REF: R.V.Reid, Ann.Phys.(NY)50(68)411-448, appendix and fig.9
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables

      double precision reidsc, r, x, h, p, A1, A2, A3, A4, alphax
     :     , u2, w2
      integer j

*    *** data, parameters and initial values

      double precision fourpi, mu, mu2, AS2, ADS2, alpha, twoalpha
      parameter (fourpi=2d0*6.2831530718d0, mu=0.7d0, mu2=mu**2 
     :     , AS2=0.7749985d0, ADS2=6.7081d-4
     :     , alpha=0.33088d0, twoalpha=0.66176d0)


      integer N
      parameter (N=33)

      double precision xx(N), u(N), w(N), du(N), dw(N)
*     tabulated u, w, du, dw
      DATA xx/ 1.000d-2, 4.125d-2, 7.250d-2, 1.350d-1, 1.975d-1,
     :         2.600d-1, 3.225d-1, 3.850d-1, 4.475d-1, 5.100d-1,
     :         5.725d-1, 6.350d-1, 6.975d-1, 7.600d-1, 8.850d-1,
     :         1.010d0 , 1.135d0 , 1.260d0 , 1.385d0 , 1.510d0 , 
     :         1.760d0 , 2.010d0 , 2.510d0 , 3.010d0 , 3.510d0 ,
     :         4.010d0 , 4.510d0 , 5.010d0 , 5.510d0 , 6.010d0 ,
     :         7.010d0 , 8.010d0 , 9.010d0 /
      DATA u / 0.0000d0, 3.3373d-5, 2.3901d-4, 2.7621d-3, 1.2737d-2,
     :         3.6062d-2, 7.5359d-2, 1.2847d-1, 1.8993d-1, 2.5349d-1,
     :         3.1390d-1, 3.6770d-1, 4.1317d-1, 4.4992d-1, 4.9953d-1, 
     :         5.2406d-1, 5.3166d-1, 5.2846d-1, 5.1926d-1, 5.0621d-1, 
     :         4.7505d-1, 4.4200d-1, 3.7864d-1, 3.2249d-1, 2.7399d-1, 
     :         2.3251d-1, 1.9719d-1, 1.6718d-1, 1.4172d-1, 1.2012d-1, 
     :         8.6290d-2, 6.1983d-2, 4.4523d-2 /
      DATA w / 0.0000d0 , 1.0850d-5, 8.4073d-5, 1.0369d-3, 4.9642d-3,
     :         1.4446d-2, 3.0795d-2, 5.3157d-2, 7.8995d-2, 1.0525d-1,
     :         1.2933d-1, 1.4958d-1, 1.6529d-1, 1.7645d-1, 1.8710d-1, 
     :         1.8654d-1, 1.7946d-1, 1.6910d-1, 1.5742d-1, 1.4553d-1, 
     :         1.2314d-1, 1.0373d-1, 7.3859d-2, 5.3293d-2, 3.9077d-2, 
     :         2.9115d-2, 2.2016d-2, 1.6871d-2, 1.3079d-2, 1.0243d-2, 
     :         6.4412d-3, 4.1575d-3, 2.7363d-3 /
      DATA du/ 2.9751d-4, 2.6127d-3, 1.2335d-2, 8.2951d-2, 2.5326d-1,
     :         5.0052d-1, 7.5072d-1, 9.3349d-1, 1.0162d0 , 1.0034d0 ,
     :         9.2042d-1, 7.9674d-1, 6.5744d-1, 5.1985d-1, 2.8494d-1, 
     :         1.1865d-1, 1.1397d-2,-5.4090d-2,-9.2523d-2,-1.1418d-1,
     :        -1.3097d-1,-1.3193d-1,-1.2004d-1,-1.0453d-1,-8.9709d-2, 
     :        -7.6510d-2,-6.5054d-2,-5.5229d-2,-4.6850d-2,-3.9726d-2,
     :        -2.8547d-2,-2.0508d-2,-1.4731d-2 /
      DATA dw/ 7.4202d-5, 8.9321d-4, 4.4755d-3, 3.1892d-2, 1.0121d-1, 
     :         2.0596d-1, 3.1477d-1, 3.9371d-1, 4.2466d-1, 4.0842d-1, 
     :         3.5772d-1, 2.8848d-1, 2.1428d-1, 1.4427d-1, 3.3294d-2, 
     :        -3.5899d-2,-7.3151d-2,-9.0091d-2,-9.5299d-2,-9.4158d-2,
     :        -8.3973d-2,-7.1282d-2,-4.9327d-2,-3.3940d-2,-2.3612d-2,
     :        -1.6691d-2,-1.2002d-2,-8.7768d-3,-6.5201d-3,-4.9136d-3,
     :        -2.8956d-3,-1.7758d-3,-1.1227d-3 /

*
*  C. ACTION
*

      x = mu * r

      if (x.lt.xx(2)) then
         reidsc = (u(2)*u(2)+w(2)*w(2))
      else if (x.gt.xx(N)) then
         alphax = alpha*x
         reidsc = AS2 * dexp(-twoalpha*x) 
     :        * ((1d0 + ADS2*(1d0+(3d0/alphax)+(3d0/(alphax**2))))**2)   
     
      else
         call locatetable(xx,N,x,j)
         h = xx(j+1) - xx(j)
         p = (x - xx(j))/h
         A1 = u(j)
         A2 = h*du(j)
         A3 = 3d0*(u(j+1)-u(j)) - (2d0*du(j)+du(j+1))*h
         A4 = 2d0*(u(j)-u(j+1)) + (du(j)+du(j+1))*h
         u2 = (A1 + A2*p + A3*p*p + A4*p*p*p)**2
         A1 = w(j)
         A2 = h*dw(j)
         A3 = 3d0*(w(j+1)-w(j)) - (2d0*dw(j)+dw(j+1))*h
         A4 = 2d0*(w(j)-w(j+1)) + (dw(j)+dw(j+1))*h
         w2 = (A1 + A2*p + A3*p*p + A4*p*p*p)**2
         reidsc = (u2+w2)
      end if

      
      return
      end


************************************************************************
* I.  Search an ordered table
      subroutine locatetable(xx,n,x,j)
*
*  A. COMMENTARY
*
*     Returns the Reid's soft core wave-function of the deuteron squared,
*     taken from: Numerical Recipes in Fortran, latest edition available 
*                 in 2002
*
*     Given an array xx(1:n), and given a value x, returns a value j
*     such that x is between xx(j) and xx(j+1). xx(1:n) must be
*     monotonic, either increasing or decreasing. j=0 or j=n is
*     returnedto indicate that x is out of range.
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables

      integer j, n, jl, jm, ju
      double precision x, xx(n)
*
*  C. ACTION
*
      jl = 0
      ju = n+1
 10   if (ju-jl.gt.1) then
         jm = (ju+jl)/2
         if ((xx(n).ge.xx(1).eqv.(x.ge.xx(jm)))) then
            jl = jm
         else
            ju = jm
         end if
         goto 10
      end if
      if (x.eq.xx(1)) then
         j = 1
      else if (x.eq.xx(n)) then
         j = n-1
      else
         j=jl
      end if

      return
      end

***********************************************************************
*
      DOUBLE PRECISION FUNCTION DGAUSS11(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS11 IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS11=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS11=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS11=DGAUSS11+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS11 ... TOO HIGH ACCURACY REQUIRED') 
      END      


