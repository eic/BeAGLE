      DOUBLE PRECISION FUNCTION DCALC(Z0)
C... MDB 2017-03-09 Initial Version
C... Use PyQM routines to calculate d = integral_z0^infinity rho(b,z) dz
C... Where rho is normalized to be 1 at r=0 for Pb.
C...
      implicit none
      include 'beagle.inc'

      double precision RR,a0,c0
      common/ woodssaxon/ RR,a0,c0
      double precision RHOOFZ, z0, zinfty, dgauss11
      external dgauss11, RHOOFZ

      zinfty = RR*RATMAX3D + 10.0d0*a0
      DCALC = DGAUSS11(RHOOFZ,z0,zinfty,1d-5)

      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RHOOFZ(Z)
C... MDB 2017-03-09 rho(b,z) for fixed b to integrate
C... MDB 2018-05-03 Generalize to 3D
C... Note: rho is normalized to be 1 at r=0.
C... Note: This won't work quite right for Ca. Have to think about it.
      implicit none
      include 'common.f'
      include 'beagle.inc'

C     Common from nucdens.f in PyQM
      double precision RR,a0,c0
      common/ woodssaxon/ RR,a0,c0

      double precision SCLFAC, DT_DENSIT, x, y, z, r, TINY6
      double precision CANGL2, Y20, Y40, RADFAC
      external SCLFAC, DT_DENSIT
      parameter (TINY6=1.0d-6)
      integer IDUM

* properties of interacting particles
c...target/proj mass, charge and projectile internal ID
      INTEGER IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,ITMMOD,
     &     MODHYP,NHYPER,IDHYP 
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5)

      r = dsqrt(z*z+bbea*bbea)
      if (USER3D) then
         x = bbea*cos(phib)
         y = bbea*sin(phib)
         CANGL2 = (x*XORIENT+y*YORIENT+z*ZORIENT)/r
         CANGL2 = CANGL2*CANGL2
         Y20 = 0.31591565D0*(3.0D0*CANGL2-1.0D0)
         Y40 = 0.10578555D0*
     &        (35.0D0*CANGL2*CANGL2-30.0D0*CANGL2+3.0D0)
         RADFAC=1.0D0+B2G3D*Y20+B4G3D*Y40
         IF (RADFAC.LE.TINY6) THEN
            WRITE(*,*)'DCALC ERROR: 3D RADIUS < 0:',RADFAC
            RADFAC=TINY6
         ENDIF
C        Effective r for a naive spherical radius (r<0 treated as 0)
         r = r+RR*(1.0-RADFAC)
      endif

C      rhoofz = DBLE(IT)*nucdens(r, ITZ, IT, 1, 1)/SCLFAC(IDUM)
      rhoofz = DENSFAC*DBLE(IT)*DT_DENSIT(ITZ,r)/SCLFAC(IDUM)

      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SCLFAC(IDUM)
      implicit none
      include 'common.f'
      include 'beagle.inc'
     
* properties of interacting particles
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5) 
c...target/proj mass, charge and projectile internal ID
      integer IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,ITMMOD
      integer MODHYP,NHYPER,IDHYP 

      double precision DT_DENSIT
      external DT_DENSIT

      double precision RR,a0,c0
      common/ woodssaxon/ RR,a0,c0

      integer idum
      logical first 
      data first /.TRUE./

      IF (first) THEN
         first = .FALSE.
         if (USER3D) then
            if ((RR-RG3D).gt.1d-3) STOP 'DCALC ERROR RR.NE.RG3D'
            if ((a0-AG3D).gt.1d-3) STOP 'DCALC ERROR a0.NE.AG3D'
            if ((c0-WG3D).gt.1d-3) STOP 'DCALC ERROR c0.NE.WG3D'
         endif
         SCLFAC = DT_DENSIT(ITZ,0.0D0)*DBLE(IT)
         WRITE(*,*) 'SCLFAC: 3A/(4piR^3)=',0.23873241D0*DBLE(IT)/RR**3,
     &        'nucleons/fm^3'
         if (USEB3D) THEN
            WRITE(*,*) 'SCLFAC: 3A/R^3(4pi+3beta2^2+3beta4^2)=',
     &           DENSFAC*0.23873241D0*DBLE(IT)/RR**3,'nucleons/fm^3'
            WRITE(*,*) 'SCLFAC: A=',IT,' Nuclear density at origin: ', 
     &           DENSFAC*SCLFAC, 'nucleons/fm^3'
         else
            WRITE(*,*) 'SCLFAC: A=',IT,' Nuclear density at origin: ', 
     &           SCLFAC, 'nucleons/fm^3'
         endif
         IF (USEB3D .OR. abs(sclfac-0.1603902d0).lt.1.0d-9)  WRITE(*,*) 
     &        'SCLFAC: Using standardized Pb density for consistency: ',
     &        '0.1603902 nucleons/fm^3'
      ENDIF
      SCLFAC = 0.1603902d0

      RETURN
      END
