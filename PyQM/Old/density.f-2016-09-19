ccccc Generate the interaction position
      subroutine InterPos
      implicit none
      include 'common.f'

* position of interacted nucleon
      double precision PosNuc
      COMMON /NPARINT/ PosNuc(4)

      double precision r
      integer i

ccccc Randomize the radius and angles
C      pos_radius = ranf(0)*quantity_table(2000)
C      pos_theta = acos(2*ranf(0)-1)
C      pos_phi = 2*pi*ranf(0)

C      r = init_dens
C      i = 1

C      do while (pos_radius.ge.quantity_table(i))
C        r = r + step_size_dens
C        i = i + 1
C      enddo
C
C      pos_radius = r + ranf(0)*step_size_dens

cccc Calculate position on cartesian axis
C      x_inter = pos_radius * sin(pos_theta) * cos(pos_phi)
C      y_inter = pos_radius * sin(pos_theta) * sin(pos_phi)
C      z_inter = pos_radius * cos(pos_theta)

c...read the interaction position(coodinate for interacted nucleon)
      x_inter = PosNuc(1)*1E12
      y_inter = PosNuc(2)*1E12
      z_inter = PosNuc(3)*1E12
      

      end

ccccc Generate the table of density in function of r
      subroutine GenNucDens(iZ,iA)
      implicit none

      include 'common.f'

      integer i,idist,irho
      integer iZ,iA
      double precision nucdens,r
      double precision integral

ccccc Some init
      QW_nb = 0
      QW_qhat = 0.

      integral = 0.
      idist = 1
      irho = 1
      step_size_dens = 0.01
      init_dens = 0.005
      r = init_dens

      do i=1,2000
        density_table(i) = nucdens(r,iZ,iA,idist,irho)
ccccc calculate the quantity of mater in function of the radius
        integral = integral + 4*pi*density_table(i)*r**2*step_size_dens
        quantity_table(i) = integral
        r = r + step_size_dens
      enddo

      write(*,*) 'Density profile of the target: DONE'
      write(*,*) 'Charge number=',iZ,' Mass number=',iA

      end
